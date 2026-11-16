"""Shared helpers for building llomics payloads (contract v0.1)."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

CONTRACT_VERSION = "0.1"
P_FLOOR = 1e-300
DEFAULT_EXCLUDE = r"^(MT-|mt-|RP[LS]\d|Rp[ls]\d|MALAT1$|Malat1$)"
_ORDER = ["pass", "warn", "fail"]


def sig(x: Any, digits: int = 4) -> float | None:
    """Round to significant digits; None for missing or non-finite values (null, never 0)."""
    if x is None:
        return None
    x = float(x)
    if math.isnan(x) or math.isinf(x):
        return None
    return float(f"{x:.{digits}g}")


def jsonable(obj: Any) -> Any:
    """Convert numpy and pandas scalars and containers to plain JSON types."""
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [jsonable(v) for v in obj]
    if hasattr(obj, "tolist") and not isinstance(obj, (str, bytes)):
        return jsonable(obj.tolist())
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def check(name: str, status: str, note: str | None = None, value: Any = None) -> dict[str, Any]:
    out: dict[str, Any] = {"check": name, "status": status}
    if value is not None:
        out["value"] = value
    if note:
        out["note"] = note
    return out


def worst(statuses: list[str]) -> str:
    return _ORDER[max(_ORDER.index(s) for s in statuses)]


def apply_trust_overrides(checks: list[dict[str, Any]], overrides: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Replace computed checks with the analyst's declarations: {"name": ("pass", "note")} or a dict."""
    if not overrides:
        return checks
    by_name = {c["check"]: c for c in checks}
    for name, decl in overrides.items():
        if isinstance(decl, dict):
            status, note = decl["status"], decl.get("note")
        else:
            status, note = decl
        by_name[name] = check(name, status, f"declared by the analyst: {note}" if note else "declared by the analyst")
    return list(by_name.values())


def p_histogram(pvals: Any) -> list[int]:
    """Counts of raw p-values in ten equal bins."""
    counts = [0] * 10
    for p in pvals:
        if p is None or p != p:
            continue
        counts[min(int(p * 10), 9)] += 1
    return counts


def p_shape(h: list[int]) -> tuple[str, str]:
    """Heuristic health check of a p-value histogram. Thresholds are stated in the note."""
    if sum(h) == 0:
        return "warn", "no p-values"
    mid = sorted(h[4:9])[2]
    later = max(h[1:])
    if later > 1.25 * h[0]:
        return "warn", (
            f"bin {h[1:].index(later) + 2} holds {later} p-values, more than 1.25x the first bin ({h[0]}): "
            "p-values are not concentrated near 0, which points to a misspecified test or variance"
        )
    if h[9] > 1.5 * mid:
        return "warn", f"excess of p-values near 1: last bin {h[9]} versus a mid-range level of {mid}"
    return "pass", "non-increasing from the first bin, flat on the right"


def write_payload(payload: dict[str, Any], out: str | Path) -> Path:
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return path
