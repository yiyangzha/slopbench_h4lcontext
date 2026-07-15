"""Inventory current-flat MC yields for the Phase 1 background hierarchy."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import uproot
from rich.logging import RichHandler


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(rich_tracebacks=True)],
)
log = logging.getLogger(__name__)

ROOT = Path(
    "/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/"
    "cms_opendata_2017_full_production/h4l_mc_nominal"
)
OUTPUT = Path("phase1_strategy/outputs/mc_preselection_inventory_v1.json")
LUMINOSITY_PB_INV = 10_000.0
XSECS_PB = {
    "GluGluToHToZZ": 6.024e-3,
    "VBF_HToZZ": 4.8794e-4,
    "ZHToZZ": 9.8394e-5,
    "WPHToZZ": 1.072352e-4,
    "WMHToZZ": 6.706e-5,
    "ZZTo4L": 1.325,
    "DYJetsToLL": 5.396e3,
    "TTBar": 52.70,
    "GGZZ2E2Mu": 3.185e-3,
    "GGZZ4Mu": 1.575e-3,
    "GGZZ4E": 1.619e-3,
}
FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}


def main() -> None:
    """Sum Metadata.nEvents and current-flat candidate counts without overwriting."""
    if OUTPUT.exists():
        raise FileExistsError(
            "Inventory already exists; preserve it and create a new version for a rerun."
        )

    samples: dict[str, dict[str, object]] = {}
    for sample, cross_section_pb in XSECS_PB.items():
        files = sorted((ROOT / sample).rglob("*.root"))
        if not files:
            raise FileNotFoundError(f"No current-flat ROOT files found for {sample}")
        n_events = 0
        n_selected = 0
        counts = {state: 0 for state in FINAL_STATES.values()}
        for path in files:
            with uproot.open(path) as source:
                n_events += int(source["Metadata"]["nEvents"].array(library="np").sum())
                arrays = source["h4lTree"].arrays(["m4l", "finalState"], library="np")
                n_selected += len(arrays["m4l"])
                in_window = (arrays["m4l"] > 105.0) & (arrays["m4l"] < 140.0)
                for code, state in FINAL_STATES.items():
                    counts[state] += int(
                        np.count_nonzero(in_window & (arrays["finalState"] == code))
                    )
        weight = cross_section_pb * LUMINOSITY_PB_INV / n_events
        samples[sample] = {
            "files": len(files),
            "cross_section_pb": cross_section_pb,
            "metadata_nEvents": n_events,
            "flat_selected_entries": n_selected,
            "entries_105_to_140_GeV": counts,
            "weight_at_10fb": weight,
            "expected_105_to_140_GeV": {
                state: count * weight for state, count in counts.items()
            },
            "expected_total_105_to_140_GeV": sum(counts.values()) * weight,
        }

    output = {
        "purpose": "Current-flat preselection inventory, not a final category or fit yield table.",
        "input_root": str(ROOT),
        "luminosity_pb_inverse": LUMINOSITY_PB_INV,
        "normalization": "cross_section_pb * luminosity_pb_inverse / sum(Metadata.nEvents)",
        "mass_window_GeV": [105.0, 140.0],
        "category_limit": "The current flat tree has no jets or FSR; VBF/non-VBF and NN-score splits require the raw-object re-ntupling route.",
        "samples": samples,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    log.info("Wrote current-flat preselection inventory for %d samples to %s.", len(samples), OUTPUT)


if __name__ == "__main__":
    main()
