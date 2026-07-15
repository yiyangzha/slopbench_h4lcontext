"""Calculate the data-statistical Z-scale floor for the Phase 1 contract."""

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

INPUT = (
    "/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/"
    "cms_opendata_2017_full_production/fake_data/fake_data_dilep_10fb.root:"
    "dilepTree"
)
OUTPUT = Path("phase1_strategy/outputs/z_scale_sanity_preflight_v1.json")
MZ_REFERENCE_GEV = 91.1876
MH_REFERENCE_GEV = 125.26
HIG_16_041_MASS_TOTAL_SYSTEMATIC_GEV = 0.08
HIG_16_041_MU_UPPER_UNCERTAINTY = 0.19
REFERENCE_LUMINOSITY_FB_INV = 35.9
TARGET_LUMINOSITY_FB_INV = 10.0
MUMU_TRIGGER_MASK = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 8)
EE_TRIGGER_MASK = (1 << 4) | (1 << 5)


def add_moments(accumulator: list[float], values: np.ndarray) -> None:
    """Accumulate count, sum, and sum of squares in chunked columnar I/O."""
    accumulator[0] += len(values)
    accumulator[1] += float(values.sum())
    accumulator[2] += float(np.square(values).sum())


def summarize_moments(accumulator: list[float]) -> dict[str, float | int]:
    """Return the conservative broad-window statistical precision proxy."""
    count, total, total_square = accumulator
    mean = total / count
    rms = float(np.sqrt(total_square / count - mean * mean))
    mean_error = rms / np.sqrt(count)
    return {
        "selected_candidates": int(count),
        "broad_window_mean_GeV": mean,
        "broad_window_rms_GeV": rms,
        "rms_over_sqrtN_GeV": mean_error,
        "relative_scale_statistical_floor": mean_error / MZ_REFERENCE_GEV,
    }


def main() -> None:
    """Write a non-result numerical scale-sanity input record."""
    if OUTPUT.exists():
        raise FileExistsError(
            "Preflight already exists; preserve it and create a new version for a rerun."
        )

    expressions = [
        "mll",
        "finalState",
        "trigBits",
        "nPV",
        "l1pfRelIso03",
        "l2pfRelIso03",
        "l1muMedium",
        "l2muMedium",
        "l1muPF",
        "l2muPF",
        "l1muGlobal",
        "l2muGlobal",
        "l1elMvaWP90",
        "l2elMvaWP90",
        "l1dz",
        "l2dz",
    ]
    moments = {"ee": [0.0, 0.0, 0.0], "mumu": [0.0, 0.0, 0.0]}
    for arrays in uproot.iterate(
        INPUT, expressions=expressions, step_size="100 MB", library="np"
    ):
        common = (
            (arrays["mll"] > 70.0)
            & (arrays["mll"] < 110.0)
            & (arrays["nPV"] > 0)
            & (np.abs(arrays["l1dz"] - arrays["l2dz"]) < 0.1)
            & (arrays["l1pfRelIso03"] < 0.15)
            & (arrays["l2pfRelIso03"] < 0.15)
        )
        mumu = (
            common
            & (arrays["finalState"] == 0)
            & ((arrays["trigBits"] & MUMU_TRIGGER_MASK) != 0)
            & (arrays["l1muMedium"] > 0)
            & (arrays["l2muMedium"] > 0)
            & ((arrays["l1muPF"] > 0) | (arrays["l1muGlobal"] > 0))
            & ((arrays["l2muPF"] > 0) | (arrays["l2muGlobal"] > 0))
        )
        ee = (
            common
            & (arrays["finalState"] == 1)
            & ((arrays["trigBits"] & EE_TRIGGER_MASK) != 0)
            & (arrays["l1elMvaWP90"] > 0)
            & (arrays["l2elMvaWP90"] > 0)
        )
        add_moments(moments["mumu"], arrays["mll"][mumu].astype(float))
        add_moments(moments["ee"], arrays["mll"][ee].astype(float))

    controls = {flavour: summarize_moments(values) for flavour, values in moments.items()}
    electron_floor = controls["ee"]["relative_scale_statistical_floor"]
    muon_floor = controls["mumu"]["relative_scale_statistical_floor"]
    scale_factor = np.sqrt(REFERENCE_LUMINOSITY_FB_INV / TARGET_LUMINOSITY_FB_INV)
    output = {
        "purpose": "Pre-fit numerical scale-reach input; this is not a fitted Z peak, calibration correction, mH, or mu result.",
        "input": INPUT,
        "control_selection": {
            "common": "OSSF pair supplied by dilep ntuplizer; 70<mll<110 GeV, nPV>0, |dz1-dz2|<0.1 cm, both relIso03<0.15.",
            "mumu": "finalState=0, one of the single/double/triple-muon trigger bits 0,1,2,3,8, both medium-ID and PF-or-global muons.",
            "ee": "finalState=1, one of electron trigger bits 4,5, both MVA-WP90 electrons.",
        },
        "method": "RMS/sqrt(N) over the broad fit domain is a conservative statistical floor proxy; the Phase 3 peak-fit covariance supersedes it.",
        "controls": controls,
        "mass_floor_from_control_statistics_GeV": {
            "4e": MH_REFERENCE_GEV * electron_floor,
            "4mu": MH_REFERENCE_GEV * muon_floor,
            "2e2mu": float(
                np.hypot(0.5 * MH_REFERENCE_GEV * electron_floor, 0.5 * MH_REFERENCE_GEV * muon_floor)
            ),
        },
        "reference_scaled_investigation_thresholds": {
            "mH_scale_contribution_GeV": HIG_16_041_MASS_TOTAL_SYSTEMATIC_GEV * scale_factor,
            "absolute_mu_scale_contribution": HIG_16_041_MU_UPPER_UNCERTAINTY * scale_factor,
            "meaning": "Cross-check thresholds only, scaled from HIG-16-041 total uncertainty; they are not assigned nuisance sizes.",
        },
        "sources": {
            "MZ_reference": "PDG 2025 Z listing, STRATEGY.md [R4]",
            "mH_and_HIG_16_041_uncertainties": "CMS HIG-16-041, STRATEGY.md [R1]",
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    log.info("Wrote numerical Z-scale preflight to %s.", OUTPUT)


if __name__ == "__main__":
    main()
