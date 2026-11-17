"""Core of llomics: plain tables in, payload units out. Standard library only.

A table is a dict of equal-length columns. Required columns: ``feature``, ``effect``, ``p``,
``p_adj``. Optional: ``symbol``, ``stat``, ``ave_expr``, ``pct_in``, ``pct_out``. Adapters for
scanpy, Seurat, limma and others only map their columns onto these names and collect provenance;
everything that decides what goes into a payload lives here, and the R package mirrors it.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .contract import CONTRACT, DEFAULTS, VERSION

Table = Mapping[str, Sequence[Any]]
OPTIONAL_ROW_COLUMNS = ("symbol", "ave_expr", "stat", "pct_in", "pct_out")
_ORDER = ["pass", "warn", "fail"]


# ---------------------------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------------------------

def _missing(x: Any) -> bool:
    return x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x)))


def sig(x: Any, digits: int | None = None) -> float | None:
    """Round to significant digits. Missing and non-finite values become None (null, never 0)."""
    if x is None:
        return None
    x = float(x)
    if math.isnan(x) or math.isinf(x):
        return None
    return float(f"{x:.{digits or DEFAULTS['digits']}g}")


def quantile(values: Sequence[float], q: float) -> float | None:
    """Linear-interpolation quantile (R type 7, the numpy default)."""
    v = sorted(values)
    if not v:
        return None
    h = (len(v) - 1) * q
    lo = math.floor(h)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])


def median(values: Sequence[float]) -> float | None:
    return quantile([v for v in values if not _missing(v)], 0.5)


def jsonable(obj: Any) -> Any:
    """Convert numpy and pandas scalars and containers to plain JSON types."""
    if isinstance(obj, Mapping):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v) for v in obj]
    if hasattr(obj, "tolist") and not isinstance(obj, (str, bytes)):
        return jsonable(obj.tolist())
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def _fmt(x: float) -> str:
    return f"{x:g}"


# ---------------------------------------------------------------------------------------------
# Trust checks
# ---------------------------------------------------------------------------------------------

def check(name: str, status: str, note: str | None = None, value: Any = None) -> dict[str, Any]:
    out: dict[str, Any] = {"check": name, "status": status}
    if value is not None:
        out["value"] = value
    if note:
        out["note"] = note
    return out


def worst(statuses: Iterable[str]) -> str:
    return _ORDER[max(_ORDER.index(s) for s in statuses)]


def apply_trust_overrides(checks: list[dict[str, Any]], overrides: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Replace computed checks with the analyst's declarations: {"name": ("pass", "note")} or a dict."""
    if not overrides:
        return checks
    by_name = {c["check"]: c for c in checks}
    for name, decl in overrides.items():
        if isinstance(decl, Mapping):
            status, note = decl["status"], decl.get("note")
        else:
            status, note = decl
        by_name[name] = check(name, status, f"declared by the analyst: {note}" if note else "declared by the analyst")
    return list(by_name.values())


def p_histogram(pvals: Iterable[Any]) -> list[int]:
    """Counts of raw p-values in ten equal bins: [0, 0.1), [0.1, 0.2), ... [0.9, 1]."""
    counts = [0] * 10
    for p in pvals:
        if not _missing(p):
            counts[min(int(p * 10), 9)] += 1
    return counts


def p_shape(h: Sequence[int]) -> tuple[str, str]:
    """Heuristic health check of a p-value histogram. The thresholds come from the contract."""
    cfg = CONTRACT["p_shape"]
    if sum(h) == 0:
        return "warn", "no p-values"
    mid = sorted(h[4:9])[2]
    later = max(h[1:])
    if later > cfg["later_bin_ratio"] * h[0]:
        return "warn", (
            f"bin {list(h[1:]).index(later) + 2} holds {later} p-values, more than "
            f"{_fmt(cfg['later_bin_ratio'])}x the first bin ({h[0]}): p-values are not concentrated near 0, "
            "which points to a misspecified test or variance"
        )
    if h[9] > cfg["last_bin_ratio"] * mid:
        return "warn", f"excess of p-values near 1: last bin {h[9]} versus a mid-range level of {mid}"
    return "pass", "non-increasing from the first bin, flat on the right"


def storey_pi0(pvals: Sequence[float]) -> float | None:
    """Storey's estimate of the share of true nulls at lambda = 0.5."""
    p = [v for v in pvals if not _missing(v)]
    if not p:
        return None
    return min(1.0, sum(v > 0.5 for v in p) / len(p) / 0.5)


def replicates_check(group_sizes: Mapping[str, int] | None, min_replicates: int | None = None,
                     units_are_replicates: bool = True, what: str = "samples") -> dict[str, Any]:
    """Enough independent replicates per group? ``units_are_replicates=False`` marks cell-level tests."""
    min_replicates = min_replicates or DEFAULTS["min_replicates"]
    if not group_sizes:
        note = f"number of {what} per group not supplied"
        if not units_are_replicates:
            note += "; the test units were not independent replicates, so p-values do not reflect variation between samples"
        return check("replicates", "not_run", note)
    smallest = min(group_sizes.values())
    listing = ", ".join(f"{k}={v}" for k, v in group_sizes.items())
    if not units_are_replicates:
        status = "fail" if smallest < 2 else "warn"
        note = (f"{what} per group: {listing}. The test units were not these replicates, so p-values do not "
                "reflect variation between samples; aggregate per sample (pseudobulk) for inference")
    else:
        status = "pass" if smallest >= min_replicates else ("warn" if smallest >= 2 else "fail")
        note = f"{what} per group: {listing}; pass needs at least {min_replicates}"
    return check("replicates", status, note, smallest)


def batch_confounding_check(crosstab: Mapping[str, Mapping[str, int]] | None, batch_in_design: bool | None = None,
                            adjustable: bool = True, single_batch: bool = False) -> dict[str, Any]:
    """Is the group effect separable from batch? ``crosstab[group][batch]`` holds observation counts."""
    name = "batch_confounding"
    if single_batch:
        return check(name, "not_applicable", "declared by the analyst as a single batch")
    if not crosstab:
        return check(name, "not_run", "no batch variable was supplied")
    batches = sorted({b for row in crosstab.values() for b, n in row.items() if n > 0})
    if len(batches) <= 1:
        return check(name, "not_applicable", "all observations come from one batch")
    groups_per_batch = [sum(1 for row in crosstab.values() if row.get(b, 0) > 0) for b in batches]
    if all(n == 1 for n in groups_per_batch):
        return check(name, "fail", "every batch contains a single group: the group effect cannot be separated from batch", crosstab)
    if not adjustable:
        return check(name, "warn", "the test has no covariates, so batch is not adjusted", crosstab)
    if batch_in_design is not True:
        note = ("not declared whether batch is in the design (set batch_in_design)" if batch_in_design is None
                else "batch is not in the design")
        return check(name, "warn", note, crosstab)
    if any(row.get(b, 0) == 0 for row in crosstab.values() for b in batches):
        return check(name, "warn", "batch is in the design, but some group x batch cells are empty", crosstab)
    return check(name, "pass", "batch is in the design and every group occurs in every batch", crosstab)


def multiple_testing_check(n_contrasts: int, method: str = "Benjamini-Hochberg") -> dict[str, Any]:
    if n_contrasts == 1:
        return check("multiple_testing", "pass", f"{method} within the one contrast reported")
    return check("multiple_testing", "warn",
                 f"{method} within each of {n_contrasts} contrasts; not adjusted across contrasts", n_contrasts)


def mean_variance_check(declared: str | None, detected: str | None = None) -> dict[str, Any]:
    name = "mean_variance"
    if declared == "none":
        return check(name, "warn", "declared: no mean-variance modelling")
    if declared == "not_applicable":
        return check(name, "not_applicable", "declared not applicable to this test or data type")
    if declared:
        return check(name, "pass", f"declared: {declared}")
    if detected:
        return check(name, "pass", f"detected: {detected}")
    return check(name, "not_run", "not declared; set mean_variance to 'voom', 'trend', 'none' or 'not_applicable'")


# ---------------------------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------------------------

def _rows(table: Table) -> list[dict[str, Any]]:
    for col in CONTRACT["feature_columns"]["required"]:
        if col not in table:
            raise KeyError(f"table is missing the required column '{col}'")
    cols = [c for c in ("feature", "effect", "p", "p_adj", *OPTIONAL_ROW_COLUMNS) if c in table]
    n = len(table["feature"])
    rows = [{c: table[c][i] for c in cols} for i in range(n)]
    rows = [r for r in rows if not _missing(r["p"]) and not _missing(r["effect"])]
    for r in rows:
        # An adjusted p-value can be missing where the raw one is not (DESeq2 independent filtering).
        # Such a feature is never eligible, and its p_adj stays null.
        if _missing(r["p_adj"]):
            r["p_adj"] = None
        if "stat" in r and _missing(r["stat"]):
            r["stat"] = 0.0
    return rows


def _below(p_adj: Any, cutoff: float) -> bool:
    return p_adj is not None and p_adj < cutoff


def _floored(p_adj: Any, p_floor: float) -> float | None:
    return None if p_adj is None else sig(max(p_adj, p_floor))


def _label(row: Mapping[str, Any]) -> str:
    return str(row["symbol"]) if row.get("symbol") not in (None, "") else str(row["feature"])


def _feature_row(row: Mapping[str, Any], p_floor: float) -> dict[str, Any]:
    out: dict[str, Any] = {"feature": str(row["feature"])}
    if "symbol" in row:
        out["symbol"] = None if row["symbol"] in (None, "") else str(row["symbol"])
    out["effect"] = sig(row["effect"])
    if "ave_expr" in row:
        out["ave_expr"] = sig(row["ave_expr"])
    if "stat" in row:
        out["stat"] = sig(row["stat"])
    out["p_adj"] = _floored(row["p_adj"], p_floor)
    for c in ("pct_in", "pct_out"):
        if c in row:
            out[c] = sig(row[c], 3)
    return out


def _rank(rows: list[dict[str, Any]], has_stat: bool, signed: int = 0) -> list[dict[str, Any]]:
    """Rank by test statistic where there is one, else by p-value then effect size. Stable for ties."""
    if has_stat:
        key = (lambda r: -abs(r["stat"])) if signed == 0 else (lambda r: -signed * r["stat"])
        return sorted(rows, key=key)
    return sorted(rows, key=lambda r: (r["p"], -abs(r["effect"])))


def count_floored(table: Table, p_floor: float | None = None) -> int:
    p_floor = p_floor or DEFAULTS["p_floor"]
    return sum(1 for v in table["p_adj"] if not _missing(v) and v < p_floor)


def contrast_unit(table: Table, *, unit_id: str, definition: str, n_features: int | None = None,
                  fdr: float | None = None, min_effect: float = 0.0, effect_cutoff: float | None = None,
                  exclude: str | None = None, panel: Sequence[str] = (), size: Mapping[str, Any] | None = None,
                  extra_facts: Mapping[str, Any] | None = None, contrast_weights: Mapping[str, float] | None = None,
                  p_shape_skip: str | None = None, p_floor: float | None = None) -> dict[str, Any]:
    """One contrast: hit counts, p-value shape, top up and down features, and the panel.

    ``p_shape_skip`` gives the reason the p-value histogram cannot be judged (for example the tool
    pre-filtered features before testing); the check is then recorded as not_run with that reason.
    """
    n_features = n_features or DEFAULTS["n_features"]
    fdr = fdr or DEFAULTS["fdr"]
    effect_cutoff = DEFAULTS["effect_cutoff"] if effect_cutoff is None else effect_cutoff
    p_floor = p_floor or DEFAULTS["p_floor"]
    pattern = re.compile(exclude) if exclude else None

    rows = _rows(table)
    has_stat = "stat" in table
    informative = [r for r in rows if not has_stat or r["stat"] != 0]
    eligible = [r for r in rows if _below(r["p_adj"], fdr) and abs(r["effect"]) > min_effect
                and not (pattern and pattern.search(_label(r)))]
    up = _rank([r for r in eligible if r["effect"] > 0], has_stat)
    down = _rank([r for r in eligible if r["effect"] < 0], has_stat)
    top_up, top_down = up[:n_features], down[:n_features]

    h = p_histogram(r["p"] for r in informative)
    status, note = ("not_run", p_shape_skip) if p_shape_skip else p_shape(h)
    effects = [r["effect"] for r in informative]
    facts: dict[str, Any] = {
        "fdr": fdr,
        "effect_cutoff": effect_cutoff,
        "n_hits": len(eligible),
        "n_up": len(up),
        "n_down": len(down),
        "n_up_effect": sum(1 for r in up if r["effect"] >= effect_cutoff),
        "n_down_effect": sum(1 for r in down if r["effect"] <= -effect_cutoff),
        "n_near_miss": sum(1 for r in rows if not _below(r["p_adj"], fdr) and _below(r["p_adj"], 2 * fdr)),
        "n_p_adj_missing": sum(1 for r in rows if r["p_adj"] is None),
        "n_uninformative": len(rows) - len(informative),
        "pi0": sig(storey_pi0([r["p"] for r in informative]), 3),
        "p_histogram": h,
        "effect_quantiles": {f"q{int(q * 100):02d}": sig(quantile(effects, q), 3) for q in (0.01, 0.25, 0.5, 0.75, 0.99)},
    }
    facts.update(extra_facts or {})

    by_feature = {str(r["feature"]): r for r in rows}
    by_symbol = {str(r["symbol"]): r for r in rows if r.get("symbol") not in (None, "")} if "symbol" in table else {}
    panel_state: dict[str, Any] = {}
    for name in panel:
        r = by_feature.get(name) or by_symbol.get(name)
        if r is None:
            panel_state[name] = None  # not measured or filtered out: null, never 0
        else:
            state = {"effect": sig(r["effect"])}
            if "ave_expr" in r:
                state["ave_expr"] = sig(r["ave_expr"])
            if "stat" in r:
                state["stat"] = sig(r["stat"])
            state["p_adj"] = _floored(r["p_adj"], p_floor)
            panel_state[name] = state

    rule = f"p_adj < {_fmt(fdr)}"
    if min_effect > 0:
        rule += f" and |effect| > {_fmt(min_effect)}"
    rule += f", ranked by {'|stat|' if has_stat else 'p then |effect|'}, up to {n_features} per direction"

    unit: dict[str, Any] = {"unit_id": str(unit_id), "unit_type": "contrast", "definition": definition}
    if contrast_weights:
        unit["contrast_weights"] = dict(contrast_weights)
    unit["size"] = {"n_tested": len(rows), **(size or {})}
    unit["facts"] = facts
    unit["top_features"] = {"up": [_feature_row(r, p_floor) for r in top_up],
                            "down": [_feature_row(r, p_floor) for r in top_down]}
    unit["panel"] = panel_state
    unit["tail"] = {"n_shown": len(top_up) + len(top_down), "n_eligible": len(eligible), "rule": rule}
    unit["trust"] = [check("pvalue_distribution", status, note)]
    return unit


def cluster_unit(table: Table, *, unit_id: str, n_features: int | None = None, fdr: float | None = None,
                 min_effect: float = 0.0, exclude: str | None = None, panel_values: Mapping[str, Any] | None = None,
                 size: Mapping[str, Any] | None = None, facts: Mapping[str, Any] | None = None,
                 p_floor: float | None = None) -> dict[str, Any]:
    """One cluster: its up-regulated markers against the rest, plus facts and panel supplied by the adapter."""
    n_features = n_features or DEFAULTS["n_features"]
    fdr = fdr or DEFAULTS["fdr"]
    p_floor = p_floor or DEFAULTS["p_floor"]
    pattern = re.compile(exclude) if exclude else None

    rows = _rows(table)
    has_stat = "stat" in table
    up = [r for r in rows if _below(r["p_adj"], fdr) and r["effect"] > min_effect
          and not (pattern and pattern.search(_label(r)))]
    up = _rank(up, has_stat, signed=1)
    top = up[:n_features]
    rule = (f"p_adj < {_fmt(fdr)} and effect > {_fmt(min_effect)}, "
            f"ranked by {'stat' if has_stat else 'p then effect'}, top {n_features}")
    return {
        "unit_id": str(unit_id),
        "unit_type": "cluster",
        "size": dict(size or {}),
        "facts": dict(facts or {}),
        "top_features": [_feature_row(r, p_floor) for r in top],
        "panel": dict(panel_values or {}),
        "tail": {"n_shown": len(top), "n_eligible": len(up), "rule": rule},
    }


def cluster_facts(labels: Sequence[str], *, qc: Mapping[str, Sequence[float]] | None = None,
                  doublet: Sequence[bool] | None = None, batch: Sequence[str] | None = None,
                  sample: Sequence[str] | None = None) -> tuple[dict[str, dict[str, Any]], dict[str, Any], list[str]]:
    """Per-cluster facts from per-observation metadata, each compared with the whole dataset.

    Returns (facts per cluster, dataset reference, clusters dominated by one batch).
    """
    labels = [str(v) for v in labels]
    n = len(labels)
    clusters = list(dict.fromkeys(labels))
    index: dict[str, list[int]] = {c: [] for c in clusters}
    for i, lab in enumerate(labels):
        index[lab].append(i)

    reference: dict[str, Any] = {"definition": "medians and fractions across all observations in the dataset"}
    for key, values in (qc or {}).items():
        reference[f"median_{key}"] = sig(median(values))
    if doublet is not None:
        reference["pct_predicted_doublet"] = sig(100 * sum(bool(d) for d in doublet) / n, 3)

    cfg = CONTRACT["batch_dominance"]
    batch_share: dict[str, float] = {}
    if batch is not None:
        batch = [str(b) for b in batch]
        for b in batch:
            batch_share[b] = batch_share.get(b, 0) + 1 / n

    facts: dict[str, dict[str, Any]] = {}
    flagged: list[str] = []
    for c in clusters:
        idx = index[c]
        f: dict[str, Any] = {}
        for key, values in (qc or {}).items():
            med = median([values[i] for i in idx])
            ref = reference[f"median_{key}"]
            f[f"median_{key}"] = sig(med)
            f[f"{key}_vs_reference"] = sig(med / ref, 3) if ref and med is not None else None
        if doublet is not None:
            pct = 100 * sum(bool(doublet[i]) for i in idx) / len(idx)
            ref = reference["pct_predicted_doublet"]
            f["pct_predicted_doublet"] = sig(pct, 3)
            f["doublet_fraction_vs_reference"] = sig(pct / ref, 3) if ref else None
        if sample is not None:
            f["n_samples_contributing"] = len({str(sample[i]) for i in idx})
        if batch is not None:
            counts: dict[str, int] = {}
            for i in idx:
                counts[batch[i]] = counts.get(batch[i], 0) + 1
            top_batch = max(sorted(counts), key=lambda b: counts[b])
            share = counts[top_batch] / len(idx)
            f["dominant_batch"] = top_batch
            f["dominant_batch_fraction"] = sig(share, 3)
            if share >= cfg["cluster_fraction"] and batch_share[top_batch] < cfg["dataset_share_below"]:
                flagged.append(f"{c} ({share:.0%} {top_batch})")
        facts[c] = f
    return facts, reference, flagged


def cluster_trust(*, test: str, n_samples: int | None, doublet_note: str | None, n_batches: int | None,
                  batch_flagged: Sequence[str] = ()) -> list[dict[str, Any]]:
    """The four required checks for a clusters payload built from a cell-level marker test."""
    cfg = CONTRACT["batch_dominance"]
    samples = (f"{n_samples} biological sample(s)" if n_samples is not None
               else "the number of biological samples is unknown (no sample column supplied)")
    checks = [
        check("pseudoreplication", "warn",
              f"{test} with cells as replicates; {samples}. P-values rank markers and are not evidence of a population-level difference."),
        check("doublet_detection", "pass", doublet_note) if doublet_note
        else check("doublet_detection", "not_run", "no predicted-doublet column was supplied"),
        check("cluster_stability", "not_run", "cannot be computed from one clustering; declare it with trust="),
    ]
    if n_batches is None:
        checks.append(check("batch_integration", "not_run", "no batch column supplied; declare a single batch with trust="))
    elif n_batches == 1:
        checks.append(check("batch_integration", "not_applicable", "all observations come from one batch"))
    elif batch_flagged:
        checks.append(check("batch_integration", "warn", "clusters dominated by one batch: " + ", ".join(batch_flagged)))
    else:
        checks.append(check("batch_integration", "pass",
                            f"no cluster is {cfg['cluster_fraction']:.0%} or more from a batch that holds under "
                            f"{cfg['dataset_share_below']:.0%} of observations"))
    return checks


# ---------------------------------------------------------------------------------------------
# Payloads
# ---------------------------------------------------------------------------------------------

def envelope(method: str, run: Mapping[str, Any], dataset: Mapping[str, Any], trust: Sequence[Mapping[str, Any]],
             units: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if method not in CONTRACT["methods"]:
        raise ValueError(f"unknown method '{method}'; the contract defines {sorted(CONTRACT['methods'])}")
    return jsonable({"llomics": VERSION, "method": method, "run": dict(run), "dataset": dict(dataset),
                     "trust": list(trust), "units": list(units)})


def standard_filters(fdr: float, n_features: int, min_effect: float, exclude: str | None, ranked_by: str,
                     extra: Sequence[Mapping[str, str]] = ()) -> list[dict[str, str]]:
    filters = [dict(f) for f in extra]
    if exclude:
        filters.append({"rule": f"top features exclude names matching {exclude}",
                        "why": "technical signal; applied to the feature list only, not to the test"})
    filters.append({"rule": f"p_adj < {_fmt(fdr)}", "why": "FDR control"})
    if min_effect > 0:
        filters.append({"rule": f"|effect| > {_fmt(min_effect)}", "why": "minimum effect size"})
    filters.append({"rule": f"top {n_features} by {ranked_by}; never padded", "why": "payload size"})
    return filters


def dge_payload(tables: Mapping[str, Table], definitions: Mapping[str, str], *, analysis_id: str,
                question: str | None, tool: str, steps: Sequence[str], n_observations: int,
                observation_type: str = "sample", group_sizes: Mapping[str, int] | None = None,
                units_are_replicates: bool = True, batch_crosstab: Mapping[str, Mapping[str, int]] | None = None,
                batch_in_design: bool | None = None, batch_adjustable: bool = True, single_batch: bool = False,
                mean_variance: str | None = None, mean_variance_detected: str | None = None,
                correction: str = "Benjamini-Hochberg", min_replicates: int | None = None,
                n_features: int | None = None, fdr: float | None = None, min_effect: float = 0.0,
                effect_cutoff: float | None = None, exclude: str | None = None, panel: Sequence[str] = (),
                sizes: Mapping[str, Mapping[str, Any]] | None = None,
                contrast_weights: Mapping[str, Mapping[str, float]] | None = None,
                filters: Sequence[Mapping[str, str]] = (), source_tables: Mapping[str, str] | None = None,
                run_extra: Mapping[str, Any] | None = None, dataset_reference: Mapping[str, Any] | None = None,
                p_shape_skip: str | None = None, trust: Mapping[str, Any] | None = None,
                p_floor: float | None = None) -> dict[str, Any]:
    """A differential-analysis payload from one standard-column table per contrast."""
    n_features = n_features or DEFAULTS["n_features"]
    fdr = fdr or DEFAULTS["fdr"]
    p_floor = p_floor or DEFAULTS["p_floor"]
    missing = [u for u in tables if u not in definitions]
    if missing:
        raise ValueError(f"definitions needs an entry in words, including the direction, for: {missing}")

    units = [
        contrast_unit(t, unit_id=uid, definition=definitions[uid], n_features=n_features, fdr=fdr,
                      min_effect=min_effect, effect_cutoff=effect_cutoff, exclude=exclude, panel=panel,
                      size=(sizes or {}).get(uid), contrast_weights=(contrast_weights or {}).get(uid),
                      p_shape_skip=p_shape_skip, p_floor=p_floor)
        for uid, t in tables.items()
    ]
    shapes = [(u["unit_id"], u["trust"][0]) for u in units]
    has_stat = all("stat" in t for t in tables.values())
    checks = [
        replicates_check(group_sizes, min_replicates, units_are_replicates),
        check("pvalue_distribution", "not_run", p_shape_skip) if p_shape_skip else
        check("pvalue_distribution", worst(c["status"] for _, c in shapes),
              " | ".join(f"{uid}: {c['note']}" for uid, c in shapes)),
        batch_confounding_check(batch_crosstab, batch_in_design, batch_adjustable, single_batch),
        multiple_testing_check(len(units), correction),
        mean_variance_check(mean_variance, mean_variance_detected),
    ]
    run = {
        "analysis_id": analysis_id,
        "question": question,
        "tool": tool,
        "steps": list(steps),
        "filters": standard_filters(fdr, n_features, min_effect, exclude,
                                    "|stat| per direction" if has_stat else "p then |effect| per direction", filters),
        "p_floor": p_floor,
        "n_p_values_floored": sum(count_floored(t, p_floor) for t in tables.values()),
        **(run_extra or {}),
        "source_tables": dict(source_tables or {}),
    }
    dataset = {
        "n_observations": n_observations,
        "observation_type": observation_type,
        "n_features": max(len(t["feature"]) for t in tables.values()),
        "reference": {"n_contrasts": len(units), **(dataset_reference or {})},
    }
    return envelope("dge", run, dataset, apply_trust_overrides(checks, trust), units)


def clusters_payload(tables: Mapping[str, Table], *, analysis_id: str, question: str | None, tool: str,
                     steps: Sequence[str], test: str, labels: Sequence[str],
                     qc: Mapping[str, Sequence[float]] | None = None, doublet: Sequence[bool] | None = None,
                     batch: Sequence[str] | None = None, sample: Sequence[str] | None = None,
                     panel_values: Mapping[str, Mapping[str, Any]] | None = None,
                     extra_facts: Mapping[str, Mapping[str, Any]] | None = None, n_dataset_features: int | None = None,
                     observation_type: str = "cell", n_features: int | None = None, fdr: float | None = None,
                     min_effect: float = 0.0, exclude: str | None = None, filters: Sequence[Mapping[str, str]] = (),
                     source_tables: Mapping[str, str] | None = None, run_extra: Mapping[str, Any] | None = None,
                     trust: Mapping[str, Any] | None = None, p_floor: float | None = None) -> dict[str, Any]:
    """A clusters payload from one standard-column marker table per cluster plus per-observation metadata."""
    n_features = n_features or DEFAULTS["n_features"]
    fdr = fdr or DEFAULTS["fdr"]
    p_floor = p_floor or DEFAULTS["p_floor"]
    labels = [str(v) for v in labels]
    n = len(labels)
    facts, reference, flagged = cluster_facts(labels, qc=qc, doublet=doublet, batch=batch, sample=sample)
    counts: dict[str, int] = {}
    for lab in labels:
        counts[lab] = counts.get(lab, 0) + 1

    units = []
    for uid, t in tables.items():
        uid = str(uid)
        f = {**facts.get(uid, {}), **(extra_facts or {}).get(uid, {})}
        units.append(cluster_unit(
            t, unit_id=uid, n_features=n_features, fdr=fdr, min_effect=min_effect, exclude=exclude,
            panel_values=(panel_values or {}).get(uid, {}), facts=f, p_floor=p_floor,
            size={"n_cells": counts.get(uid, 0), "pct_of_total": sig(100 * counts.get(uid, 0) / n, 3)}))

    doublet_note = None
    if doublet is not None:
        doublet_note = f"{_fmt(reference['pct_predicted_doublet'])}% of observations predicted doublets"
    checks = cluster_trust(
        test=test, n_samples=len(set(map(str, sample))) if sample is not None else None, doublet_note=doublet_note,
        n_batches=len(set(map(str, batch))) if batch is not None else None, batch_flagged=flagged)
    has_stat = all("stat" in t for t in tables.values())
    run = {
        "analysis_id": analysis_id,
        "question": question,
        "tool": tool,
        "steps": list(steps),
        "filters": standard_filters(fdr, n_features, min_effect, exclude, "stat" if has_stat else "p then effect", filters),
        "p_floor": p_floor,
        "n_p_values_floored": sum(count_floored(t, p_floor) for t in tables.values()),
        "panel_units": "mean expression across the unit's observations, in the units of the tested matrix",
        **(run_extra or {}),
        "source_tables": dict(source_tables or {}),
    }
    dataset = {
        "n_observations": n,
        "observation_type": observation_type,
        "n_biological_samples": len(set(map(str, sample))) if sample is not None else None,
        "n_features": n_dataset_features or max(len(t["feature"]) for t in tables.values()),
        "reference": reference,
    }
    return envelope("clusters", run, dataset, apply_trust_overrides(checks, trust), units)


def write_payload(payload: Mapping[str, Any], out: str | Path) -> Path:
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return path
