#!/usr/bin/env python3
"""Validate an llomics payload against the contract. Standard library only.

Usage:
    python validate_payload.py payload.json [--unit-budget 8000]

Reads llomics-contract.json from the directory this file is in. Exit code 0 means no errors
(warnings may still be printed); 1 means errors.

This file is the same in the Python package (llomics/validate.py) and in the llomics-core skill
(scripts/validate_payload.py); tools/sync.py keeps the copies identical.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

CONTRACT = json.loads(Path(__file__).with_name("llomics-contract.json").read_text(encoding="utf-8"))
P_KEYS = {k.lower() for k in CONTRACT["p_value_keys"]}


def walk(node: Any, path: str, errors: list[str]) -> None:
    """Recursively check numbers: no NaN/inf anywhere, p-values in (0, 1]."""
    if isinstance(node, dict):
        for key, value in node.items():
            child = f"{path}.{key}"
            if key.lower() in P_KEYS and isinstance(value, (int, float)) and not isinstance(value, bool):
                if value == 0:
                    errors.append(f"{child}: p-value is exactly 0 (underflow); clamp to run.p_floor")
                elif not (0 < value <= 1):
                    errors.append(f"{child}: p-value {value} is outside (0, 1]")
            walk(value, child, errors)
    elif isinstance(node, list):
        for i, value in enumerate(node):
            walk(value, f"{path}[{i}]", errors)
    elif isinstance(node, float) and (math.isnan(node) or math.isinf(node)):
        errors.append(f"{path}: NaN or infinity; use null and say why")


def feature_rows(top_features: Any) -> list[Any] | None:
    """Return all rows from a list, or from an object of lists; None if the shape is wrong."""
    if isinstance(top_features, list):
        return top_features
    if isinstance(top_features, dict) and all(isinstance(v, list) for v in top_features.values()):
        return [row for rows in top_features.values() for row in rows]
    return None


def validate(payload: Any, unit_budget: int | None = None) -> tuple[list[str], list[str]]:
    """Return (errors, warnings) for a parsed payload."""
    unit_budget = unit_budget or CONTRACT["defaults"]["unit_budget_chars"]
    errors: list[str] = []
    warnings: list[str] = []

    if not isinstance(payload, dict):
        return ["payload must be a JSON object"], warnings
    for key in CONTRACT["top_level_keys"]:
        if key not in payload:
            errors.append(f"missing top-level key: {key}")
    if errors:
        return errors, warnings

    method = payload["method"]
    run, dataset, trust, units = payload["run"], payload["dataset"], payload["trust"], payload["units"]
    if str(payload["llomics"]) != CONTRACT["version"]:
        warnings.append(f"payload is contract {payload['llomics']}; this validator is {CONTRACT['version']}")

    # L0 run card
    for key in CONTRACT["run_keys"]:
        if key not in run:
            errors.append(f"run.{key} is missing")
    if isinstance(run.get("filters"), list):
        for i, flt in enumerate(run["filters"]):
            if not isinstance(flt, dict) or "rule" not in flt:
                errors.append(f"run.filters[{i}] needs a 'rule'")
    if not run.get("source_tables"):
        warnings.append("run.source_tables is empty: no L2 drill-down path for dropped rows")
    if not run.get("question"):
        warnings.append("run.question is empty: the model is not told what the analysis is for")

    for key in CONTRACT["dataset_keys"]:
        if key not in dataset:
            errors.append(f"dataset.{key} is missing")
    if "reference" in dataset and not isinstance(dataset["reference"], dict):
        errors.append("dataset.reference must be an object (may be empty)")

    # Trust layer
    statuses = set(CONTRACT["statuses"])
    seen: dict[str, str] = {}
    if not isinstance(trust, list):
        errors.append("trust must be a list of checks")
        trust = []
    for i, chk in enumerate(trust):
        if not isinstance(chk, dict) or "check" not in chk or "status" not in chk:
            errors.append(f"trust[{i}] needs 'check' and 'status'")
            continue
        status = chk["status"]
        if status not in statuses:
            errors.append(f"trust[{i}] ({chk['check']}): status '{status}' not in {sorted(statuses)}")
        if status in CONTRACT["note_required_for"] and not chk.get("note"):
            errors.append(f"trust[{i}] ({chk['check']}): status '{status}' needs a 'note'")
        seen[chk["check"]] = status

    spec = CONTRACT["methods"].get(method)
    if spec:
        for name in spec["required_checks"]:
            if name not in seen:
                errors.append(f"trust: required check '{name}' for method '{method}' is absent (use not_run)")
    else:
        errors.append(f"method '{method}' is not in the contract; known methods: {sorted(CONTRACT['methods'])}")

    not_run = sorted(name for name, status in seen.items() if status == "not_run")
    failed = sorted(name for name, status in seen.items() if status == "fail")
    if not_run:
        warnings.append(f"checks not run: {', '.join(not_run)}")
    if failed:
        warnings.append(f"checks failed: {', '.join(failed)}")

    # L1 unit states
    if not isinstance(units, list) or not units:
        errors.append("units must be a non-empty list")
        units = []
    ids: set[str] = set()
    unknown_tail: list[str] = []
    for i, unit in enumerate(units):
        where = f"units[{i}]"
        if not isinstance(unit, dict):
            errors.append(f"{where} must be an object")
            continue
        for key in CONTRACT["unit_keys"]:
            if key not in unit:
                errors.append(f"{where}.{key} is missing")
        uid = str(unit.get("unit_id"))
        if uid in ids:
            errors.append(f"{where}: duplicate unit_id '{uid}'")
        ids.add(uid)
        if spec and "unit_type" in unit and unit["unit_type"] not in spec["unit_types"]:
            errors.append(f"{where} ({uid}): unit_type '{unit['unit_type']}' is not valid for method '{method}' "
                          f"(expected one of {spec['unit_types']})")

        rows = feature_rows(unit.get("top_features"))
        if rows is None:
            if "top_features" in unit:
                errors.append(f"{where}.top_features must be a list, or an object of lists")
            rows = []
        for j, row in enumerate(rows):
            if not isinstance(row, dict) or "feature" not in row:
                errors.append(f"{where}.top_features row {j} needs a 'feature'")

        tail = unit.get("tail")
        if isinstance(tail, dict):
            if any(k not in tail for k in CONTRACT["tail_keys"]):
                errors.append(f"{where}.tail needs 'n_shown' and 'n_eligible' (n_eligible may be null)")
            else:
                if tail["n_shown"] != len(rows):
                    errors.append(f"{where}.tail.n_shown={tail['n_shown']} but {len(rows)} rows are present")
                if tail["n_eligible"] is None:
                    unknown_tail.append(uid)
                elif tail["n_eligible"] < len(rows):
                    errors.append(f"{where}.tail.n_eligible is smaller than the number of rows shown")
        elif "tail" in unit:
            errors.append(f"{where}.tail must be an object")

        if "panel" not in unit:
            warnings.append(f"{where} ({uid}): no panel; add expected features or an empty object")

        size = len(json.dumps(unit, separators=(",", ":")))
        if size > unit_budget:
            warnings.append(f"{where} ({uid}): {size} chars exceeds the unit budget of {unit_budget}")

    if unknown_tail:
        warnings.append(
            f"tail size unknown (n_eligible is null) in {len(unknown_tail)} of {len(units)} units: "
            + ", ".join(unknown_tail)
        )

    walk(payload, "payload", errors)
    return errors, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate an llomics payload against the contract.")
    parser.add_argument("payload", type=Path)
    parser.add_argument("--unit-budget", type=int, default=None, help="max compact-JSON chars per unit")
    args = parser.parse_args()

    try:
        payload = json.loads(args.payload.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"ERROR  cannot read payload: {exc}")
        return 1

    errors, warnings = validate(payload, args.unit_budget)
    for msg in warnings:
        print(f"WARN   {msg}")
    for msg in errors:
        print(f"ERROR  {msg}")

    n_units = len(payload.get("units", [])) if isinstance(payload, dict) else 0
    verdict = "INVALID" if errors else "VALID"
    print(f"{verdict}: {n_units} units, {len(errors)} errors, {len(warnings)} warnings")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
