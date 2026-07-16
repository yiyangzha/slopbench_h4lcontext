"""Persist the ggH response write/close and raw-content validation result."""

from __future__ import annotations

import json
import logging
from datetime import UTC, datetime
from pathlib import Path

from validate_raw_pseudodata_regeneration_v1 import (
    REPAIRS,
    check_response_sidecars,
    compare_raw_pair,
)


log = logging.getLogger(__name__)

SAMPLE = "GluGluToHToZZ"
OUTPUT = Path("phase1_strategy/outputs/raw_pseudodata_regeneration_ggh_raw_v1.json")


def main() -> None:
    """Validate the terminated ggH response task without replacing evidence."""
    if OUTPUT.exists():
        raise FileExistsError("ggH raw validation evidence already exists; create a new version for a rerun")

    spec = REPAIRS[SAMPLE]
    response_check = check_response_sidecars(SAMPLE, spec)
    raw_check = compare_raw_pair(SAMPLE, spec)
    result = {
        "schema_version": 1,
        "generated_utc": datetime.now(UTC).isoformat(),
        "status": "pass",
        "scope": "ggH-only regenerated modified-MC raw validation; not collision-data provenance.",
        "response_validation": response_check,
        "regenerated_raw_validation": raw_check,
        "next_required_action": "Re-ntuplize the versioned ggH raw file before the combined 652-row validation.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    log.info("Wrote %s.", OUTPUT)


if __name__ == "__main__":
    main()
