"""Validate the DY repair stage before any ggH response submission."""

from __future__ import annotations

import json
import logging
from datetime import UTC, datetime
from pathlib import Path

from validate_raw_pseudodata_regeneration_v1 import (
    REPAIRS,
    compare_h4l_pair,
    h4l_sources_with_repairs,
    resolve_fake_sources,
)


log = logging.getLogger(__name__)

SAMPLE = "DYJetsToLL"
OUTPUT = Path("phase1_strategy/outputs/raw_pseudodata_regeneration_dy_gate_v1.json")
MAPPING_OUTPUT = Path(
    "phase1_strategy/outputs/raw_pseudodata_regeneration_dy_gate_mapping_v1.json"
)
RAW_VALIDATION = Path("phase1_strategy/outputs/raw_pseudodata_regeneration_dy_raw_v1.json")


def main() -> None:
    """Persist the DY raw/H4l/pseudo-data gate without overwriting evidence."""
    if OUTPUT.exists() or MAPPING_OUTPUT.exists():
        raise FileExistsError("DY gate evidence already exists; create a new version for a rerun")

    if not RAW_VALIDATION.is_file():
        raise FileNotFoundError(f"DY raw validation evidence is missing: {RAW_VALIDATION}")
    raw_validation = json.loads(RAW_VALIDATION.read_text())
    if raw_validation.get("status") != "pass":
        raise RuntimeError("DY raw validation did not pass")

    spec = REPAIRS[SAMPLE]
    h4l_check = compare_h4l_pair(SAMPLE, spec)

    h4l_sources, replacements = h4l_sources_with_repairs({SAMPLE})
    fake_arrays, mapping, h4l_source_scan = resolve_fake_sources(h4l_sources, replacements)
    affected = set(spec["affected_fake_identities"])
    affected_rows = [
        row
        for row in mapping
        if (row["run"], row["lumi"], row["event"]) in affected
    ]
    if len(affected_rows) != len(affected):
        raise RuntimeError("The three DY pseudo-data rows were not all resolved")
    repaired_h4l = str(spec["repaired_h4l"])
    if any(row["source_h4l_file"] != repaired_h4l for row in affected_rows):
        raise RuntimeError("A DY pseudo-data row did not resolve to the regenerated H4l source")

    mapping_payload = {
        "schema_version": 1,
        "purpose": "DY-only gate before the ggH response repair; each listed row is an exact 111-field pseudo-data/H4l payload match.",
        "fake_rows_examined_for_full_payload_resolution": len(mapping),
        "h4l_payload_fields_compared": len(fake_arrays),
        "mixer_h4l_source_scan": h4l_source_scan,
        "affected_dy_rows": affected_rows,
    }
    result = {
        "schema_version": 1,
        "generated_utc": datetime.now(UTC).isoformat(),
        "status": "pass",
        "scope": "DY-only synthetic modified-MC pseudo-data repair gate; not collision-data provenance.",
        "dy_raw_validation_evidence": str(RAW_VALIDATION),
        "regenerated_h4l_validation": h4l_check,
        "affected_dy_pseudodata_agreement": {
            "affected_rows": len(affected_rows),
            "exact_full_payload_matches": len(affected_rows),
            "payload_fields_compared": len(fake_arrays),
            "source_token_manifest": str(MAPPING_OUTPUT),
            "mixer_h4l_source_scan": h4l_source_scan,
        },
        "next_authorized_action": "Submit exactly one versioned ggH response task only after this pass record exists.",
    }
    MAPPING_OUTPUT.write_text(json.dumps(mapping_payload, indent=2) + "\n")
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    log.info("Wrote %s and %s.", OUTPUT, MAPPING_OUTPUT)


if __name__ == "__main__":
    main()
