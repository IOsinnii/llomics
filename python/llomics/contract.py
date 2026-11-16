"""The llomics contract, loaded from the JSON file that R, Python and the validator all share."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CONTRACT_PATH = Path(__file__).with_name("llomics-contract.json")
CONTRACT: dict[str, Any] = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))

VERSION: str = CONTRACT["version"]
DEFAULTS: dict[str, Any] = CONTRACT["defaults"]
METHODS: dict[str, Any] = CONTRACT["methods"]
