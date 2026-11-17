#!/usr/bin/env python3
"""Copy the single-source files to the places that need their own copy, or check the copies.

    python tools/sync.py           # copy
    python tools/sync.py --check   # exit 1 if any copy differs (used by the tests)

The contract is authored in contract/llomics-contract.json. The validator is authored in
python/llomics/validate.py. Skills are installed on their own, and the R and Python packages are
installed on their own, so each needs a copy; none of the copies is edited by hand.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COPIES = {
    "contract/llomics-contract.json": [
        "python/llomics/llomics-contract.json",
        "r/llomics/inst/llomics-contract.json",
        "skills/llomics-core/scripts/llomics-contract.json",
    ],
    "python/llomics/validate.py": [
        "skills/llomics-core/scripts/validate_payload.py",
        "r/llomics/inst/validate_payload.py",
    ],
}


def main() -> int:
    check = "--check" in sys.argv
    stale = []
    for source, targets in COPIES.items():
        content = (ROOT / source).read_bytes().replace(b"\r\n", b"\n")
        for target in targets:
            path = ROOT / target
            current = path.read_bytes().replace(b"\r\n", b"\n") if path.exists() else None
            if current == content:
                continue
            stale.append(target)
            if not check:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
                print(f"synced {target}")
    if check and stale:
        print("out of sync (run python tools/sync.py): " + ", ".join(stale))
        return 1
    if not stale:
        print("all copies are in sync")
    return 0


if __name__ == "__main__":
    sys.exit(main())
