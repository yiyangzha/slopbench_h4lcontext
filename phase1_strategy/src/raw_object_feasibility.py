"""Write a matched raw-object VBF/FSR pilot from the supplied NanoAOD MC."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import awkward as ak
import numpy as np
import uproot
from rich.logging import RichHandler


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(rich_tracebacks=True)],
)
log = logging.getLogger(__name__)

RAW_INPUT = Path(
    "/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/"
    "cms_opendata_2017_full_production/nanoaod_modified/VBF_HToZZ/230000/"
    "3602787B-E4E1-784B-95EE-663C48DB1CE6_scale_smear.root"
)
FLAT_INPUT = Path(
    "/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/"
    "cms_opendata_2017_full_production/h4l_mc_modified/VBF_HToZZ/230000/"
    "3602787B-E4E1-784B-95EE-663C48DB1CE6_scale_smear_h4l.root"
)
OUTPUT_ROOT = Path("phase1_strategy/outputs/raw_object_pilot_vbf_v1.root")
OUTPUT_JSON = Path("phase1_strategy/outputs/raw_object_pilot_vbf_v1.json")

# These are the HIG-16-041 jet-object baseline requirements used only to
# demonstrate raw-object availability. The final VBF category is the
# reference-like D_2jet category documented in STRATEGY_v2.md.
JET_PT_MIN_GEV = 30.0
JET_ABS_ETA_MAX = 4.7
LEPTON_JET_DR2_MIN = 0.16


def event_keys(run: np.ndarray, lumi: np.ndarray, event: np.ndarray) -> np.ndarray:
    """Return sortable, exact run/lumi/event structured keys."""
    return np.rec.fromarrays(
        [run, lumi, event],
        names=("run", "lumi", "event"),
    )


def raw_to_flat_matches(
    raw_keys: np.ndarray, flat_keys: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Match each selected raw event to exactly one existing flat candidate."""
    if len(np.unique(flat_keys)) != len(flat_keys):
        raise RuntimeError("Flat pilot input has duplicate run/lumi/event candidates")

    raw_order = np.argsort(raw_keys)
    flat_order = np.argsort(flat_keys)
    sorted_raw = raw_keys[raw_order]
    sorted_flat = flat_keys[flat_order]
    positions = np.searchsorted(sorted_flat, sorted_raw)
    valid = positions < len(sorted_flat)
    comparison_positions = np.minimum(positions, len(sorted_flat) - 1)
    valid &= sorted_flat[comparison_positions] == sorted_raw
    return raw_order[valid], flat_order[positions[valid]]


def candidate_lepton_coordinates(
    flat: dict[str, np.ndarray], indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Stack the four preserved flat-candidate lepton directions."""
    eta = np.column_stack([flat[f"l{index}eta"][indices] for index in range(1, 5)])
    phi = np.column_stack([flat[f"l{index}phi"][indices] for index in range(1, 5)])
    return eta, phi


def clean_selected_jets(
    raw: ak.Array, lepton_eta: np.ndarray, lepton_phi: np.ndarray
) -> tuple[ak.Array, ak.Array, ak.Array, ak.Array, ak.Array]:
    """Select and clean jets while retaining their b-tag discriminator."""
    jet_pt = raw["Jet_pt"]
    jet_eta = raw["Jet_eta"]
    jet_phi = raw["Jet_phi"]
    baseline = (
        (jet_pt > JET_PT_MIN_GEV)
        & (abs(jet_eta) < JET_ABS_ETA_MAX)
        & (raw["Jet_jetId"] >= 2)
    )
    deta = jet_eta[:, :, np.newaxis] - lepton_eta[:, np.newaxis, :]
    dphi = np.arctan2(
        np.sin(jet_phi[:, :, np.newaxis] - lepton_phi[:, np.newaxis, :]),
        np.cos(jet_phi[:, :, np.newaxis] - lepton_phi[:, np.newaxis, :]),
    )
    clean = ~ak.any(deta * deta + dphi * dphi < LEPTON_JET_DR2_MIN, axis=2)
    mask = baseline & clean
    return (
        jet_pt[mask],
        jet_eta[mask],
        jet_phi[mask],
        raw["Jet_mass"][mask],
        raw["Jet_btagDeepFlavB"][mask],
    )


def leading_dijet_mass(pt: ak.Array, eta: ak.Array, phi: ak.Array, mass: ak.Array) -> ak.Array:
    """Calculate m_jj for the two leading retained jets, or NaN if absent."""
    order = ak.argsort(pt, axis=1, ascending=False)
    pt = pt[order]
    eta = eta[order]
    phi = phi[order]
    mass = mass[order]
    padded_pt = ak.pad_none(pt, 2, clip=True)
    padded_eta = ak.pad_none(eta, 2, clip=True)
    padded_phi = ak.pad_none(phi, 2, clip=True)
    padded_mass = ak.pad_none(mass, 2, clip=True)
    px = padded_pt * np.cos(padded_phi)
    py = padded_pt * np.sin(padded_phi)
    pz = padded_pt * np.sinh(padded_eta)
    energy = np.sqrt((padded_pt * np.cosh(padded_eta)) ** 2 + padded_mass**2)
    mass_squared = (
        (energy[:, 0] + energy[:, 1]) ** 2
        - (px[:, 0] + px[:, 1]) ** 2
        - (py[:, 0] + py[:, 1]) ** 2
        - (pz[:, 0] + pz[:, 1]) ** 2
    )
    return ak.fill_none(np.sqrt(ak.where(mass_squared > 0, mass_squared, 0)), np.nan)


def main() -> None:
    """Create the non-overwriting raw-object feasibility artifact."""
    if OUTPUT_ROOT.exists() or OUTPUT_JSON.exists():
        raise FileExistsError(
            "Pilot outputs already exist; preserve them and create a new version "
            "for a subsequent rerun."
        )

    required = {
        "event_identity": ["run", "luminosityBlock", "event"],
        "jets_and_btag": [
            "Jet_pt",
            "Jet_eta",
            "Jet_phi",
            "Jet_mass",
            "Jet_jetId",
            "Jet_rawFactor",
            "Jet_btagDeepFlavB",
        ],
        "fsr": [
            "FsrPhoton_pt",
            "FsrPhoton_eta",
            "FsrPhoton_phi",
            "FsrPhoton_relIso03",
            "FsrPhoton_muonIdx",
            "Photon_pt",
            "Photon_eta",
            "Photon_phi",
            "Photon_pfRelIso03_all",
        ],
    }
    raw_file = uproot.open(RAW_INPUT)
    flat_file = uproot.open(FLAT_INPUT)
    raw_tree = raw_file["Events"]
    flat_tree = flat_file["h4lTree"]
    raw_branches = set(raw_tree.keys())
    missing = {
        group: [branch for branch in branches if branch not in raw_branches]
        for group, branches in required.items()
    }
    if any(missing.values()):
        raise RuntimeError(f"Raw NanoAOD is missing required pilot branches: {missing}")

    raw_ids = raw_tree.arrays(required["event_identity"], library="np")
    flat_fields = ["run", "lumi", "event", "m4l"]
    flat_fields += [f"l{index}{quantity}" for index in range(1, 5) for quantity in ("eta", "phi")]
    flat = flat_tree.arrays(flat_fields, library="np")
    raw_keys = event_keys(
        raw_ids["run"], raw_ids["luminosityBlock"], raw_ids["event"]
    )
    flat_keys = event_keys(flat["run"], flat["lumi"], flat["event"])
    raw_indices, flat_indices = raw_to_flat_matches(raw_keys, flat_keys)
    if len(raw_indices) != len(flat_keys):
        raise RuntimeError(
            f"Only {len(raw_indices)} of {len(flat_keys)} flat candidates match raw events"
        )

    raw_fields = sorted({branch for branches in required.values() for branch in branches})
    raw = raw_tree.arrays(raw_fields, library="ak")[raw_indices]
    lepton_eta, lepton_phi = candidate_lepton_coordinates(flat, flat_indices)
    selected_pt, selected_eta, selected_phi, selected_mass, selected_btag = clean_selected_jets(
        raw, lepton_eta, lepton_phi
    )
    n_selected_jets = ak.num(selected_pt, axis=1)
    m_jj = leading_dijet_mass(selected_pt, selected_eta, selected_phi, selected_mass)
    fsr_muon_idx = raw["FsrPhoton_muonIdx"]
    n_fsr = ak.num(raw["FsrPhoton_pt"], axis=1)

    pilot = {
        "run": raw["run"],
        "lumi": raw["luminosityBlock"],
        "event": raw["event"],
        "m4l_flat": flat["m4l"][flat_indices],
        "nSelectedJet": n_selected_jets,
        "selectedJet_pt": selected_pt,
        "selectedJet_eta": selected_eta,
        "selectedJet_phi": selected_phi,
        "selectedJet_mass": selected_mass,
        "selectedJet_btagDeepFlavB": selected_btag,
        "selectedJet_mjj": m_jj,
        "nFsrPhoton": n_fsr,
        "fsrPhoton_pt": raw["FsrPhoton_pt"],
        "fsrPhoton_eta": raw["FsrPhoton_eta"],
        "fsrPhoton_phi": raw["FsrPhoton_phi"],
        "fsrPhoton_relIso03": raw["FsrPhoton_relIso03"],
        "fsrPhoton_muonIdx": fsr_muon_idx,
        "nPhoton": ak.num(raw["Photon_pt"], axis=1),
        "photon_pt": raw["Photon_pt"],
        "photon_eta": raw["Photon_eta"],
        "photon_phi": raw["Photon_phi"],
        "photon_relIso03": raw["Photon_pfRelIso03_all"],
    }
    with uproot.recreate(OUTPUT_ROOT) as output:
        output["rawObjectPilot"] = pilot

    in_mass_window = (flat["m4l"][flat_indices] > 105.0) & (
        flat["m4l"][flat_indices] < 140.0
    )
    summary = {
        "purpose": "MC raw-object re-ntupling feasibility pilot; not a VBF yield result",
        "raw_input": str(RAW_INPUT),
        "flat_input": str(FLAT_INPUT),
        "output_root": str(OUTPUT_ROOT),
        "raw_events": int(len(raw_keys)),
        "flat_candidates": int(len(flat_keys)),
        "matched_candidates": int(len(raw_indices)),
        "match_fraction": float(len(raw_indices) / len(flat_keys)),
        "candidate_events_in_105_140_GeV": int(np.count_nonzero(in_mass_window)),
        "candidate_events_with_at_least_two_clean_jets": int(
            np.count_nonzero(ak.to_numpy(n_selected_jets >= 2))
        ),
        "mass_window_events_with_at_least_two_clean_jets": int(
            np.count_nonzero(ak.to_numpy(n_selected_jets[in_mass_window] >= 2))
        ),
        "candidate_events_with_fsr_photon": int(
            np.count_nonzero(ak.to_numpy(n_fsr > 0))
        ),
        "candidate_events_with_muon_fsr_association": int(
            np.count_nonzero(ak.to_numpy(ak.any(fsr_muon_idx >= 0, axis=1)))
        ),
        "required_branch_groups": required,
        "missing_required_branches": missing,
        "jet_object_probe": {
            "pt_min_GeV": JET_PT_MIN_GEV,
            "abs_eta_max": JET_ABS_ETA_MAX,
            "jetId_min": 2,
            "lepton_jet_deltaR_min": float(np.sqrt(LEPTON_JET_DR2_MIN)),
        },
        "limits": [
            "This verifies a matched MC raw-object route only.",
            "No raw-data NanoAOD provenance is supplied by the production root.",
            "The final VBF category requires the separate reference-like D_2jet implementation and expected-yield study.",
            "The final FSR algorithm must validate muon-indexed FsrPhoton and electron-associated Photon recovery.",
        ],
    }
    OUTPUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    log.info(
        "Wrote %s matched raw-object candidates to %s (%.1f%% flat-to-raw match).",
        len(raw_indices),
        OUTPUT_ROOT,
        100.0 * summary["match_fraction"],
    )


if __name__ == "__main__":
    main()
