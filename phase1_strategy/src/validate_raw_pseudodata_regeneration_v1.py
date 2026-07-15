"""Validate the versioned repair of two modified NanoAOD pseudo-data sources."""

from __future__ import annotations

import hashlib
import json
import logging
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

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

PRODUCTION_ROOT = Path(
    "/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/"
    "cms_opendata_2017_full_production"
)
REPAIR_ROOT = Path(f"{PRODUCTION_ROOT}_regenerated_raw_v1")
FAKE_DATA = PRODUCTION_ROOT / "fake_data/fake_data_10fb.root"
MIXER_MANIFEST = Path(
    "/eos/home-y/yiyangz/codex/slopbench_code/results/"
    "cms_opendata_2017_full_production/modified_h4l_mixer_manifest.json"
)
OUTPUT = Path("phase1_strategy/outputs/raw_pseudodata_regeneration_v1.json")
MAPPING_OUTPUT = Path(
    "phase1_strategy/outputs/raw_pseudodata_regeneration_mapping_v1.json"
)
CONDOR_ROOT = Path("phase1_strategy/outputs/raw_pseudodata_regeneration_v1")

KEY_DTYPE = np.dtype([("run", "<u8"), ("lumi", "<u8"), ("event", "<u8")])
IDENTITY_FIELDS = ["run", "luminosityBlock", "event"]
FAKE_IDENTITY_FIELDS = ["run", "lumi", "event"]
EXPECTED_MUTATED_FIELDS = {
    "Muon_pt",
    "Muon_eta",
    "Muon_phi",
    "Electron_pt",
    "Electron_eta",
    "Electron_phi",
}
REQUIRED_SCHEMA = {
    "identity": ["run", "luminosityBlock", "event"],
    "jets_raw_jec_jer_deepjet_inputs": [
        "Jet_pt",
        "Jet_eta",
        "Jet_phi",
        "Jet_mass",
        "Jet_jetId",
        "Jet_rawFactor",
        "Jet_area",
        "Jet_muonSubtrFactor",
        "Jet_genJetIdx",
        "Jet_btagDeepFlavB",
    ],
    "fsr_photons": [
        "FsrPhoton_pt",
        "FsrPhoton_eta",
        "FsrPhoton_phi",
        "FsrPhoton_relIso03",
        "FsrPhoton_muonIdx",
    ],
    "generic_photons": [
        "Photon_pt",
        "Photon_eta",
        "Photon_phi",
        "Photon_pfRelIso03_all",
    ],
    "muons": ["nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"],
    "electrons": [
        "nElectron",
        "Electron_pt",
        "Electron_eta",
        "Electron_phi",
        "Electron_mass",
    ],
}

REPAIRS: dict[str, dict[str, Any]] = {
    "DYJetsToLL": {
        "entries": 1_283_487,
        "response_config": Path("phase1_strategy/outputs/response_repair_dy_v1.json"),
        "original_raw": Path(
            "/eos/opendata/cms/mc/RunIISummer20UL17NanoAODv9/"
            "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/"
            "106X_mc2017_realistic_v9-v1/280000/"
            "3183559F-C042-5944-84D4-63F2E7CF0D0F.root"
        ),
        "repaired_raw": REPAIR_ROOT
        / "nanoaod_modified/DYJetsToLL/280000/"
        "3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear_regenerated_v1.root",
        "original_h4l": PRODUCTION_ROOT
        / "h4l_mc_modified/DYJetsToLL/280000/"
        "3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear_h4l.root",
        "repaired_h4l": REPAIR_ROOT
        / "h4l_mc_modified/DYJetsToLL/280000/"
        "3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear_regenerated_v1_h4l.root",
        "affected_fake_identities": [
            (1, 390_917, 97_729_202),
            (1, 927_325, 231_831_246),
            (1, 672_725, 168_181_225),
        ],
    },
    "GluGluToHToZZ": {
        "entries": 827_970,
        "response_config": Path("phase1_strategy/outputs/response_repair_ggh_v1.json"),
        "original_raw": Path(
            "/eos/opendata/cms/mc/RunIISummer20UL17NanoAODv9/"
            "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_minloHJJ_"
            "JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/"
            "2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A.root"
        ),
        "repaired_raw": REPAIR_ROOT
        / "nanoaod_modified/GluGluToHToZZ/2430000/"
        "78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear_regenerated_v1.root",
        "original_h4l": PRODUCTION_ROOT
        / "h4l_mc_modified/GluGluToHToZZ/2430000/"
        "78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear_h4l.root",
        "repaired_h4l": REPAIR_ROOT
        / "h4l_mc_modified/GluGluToHToZZ/2430000/"
        "78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear_regenerated_v1_h4l.root",
        "affected_fake_identities": [(1, 559, 215_650), (1, 1_739, 670_997)],
    },
}


def event_keys(run: np.ndarray, lumi: np.ndarray, event: np.ndarray) -> np.ndarray:
    """Return sortable exact run/lumi/event keys with one fixed dtype."""
    keys = np.empty(len(run), dtype=KEY_DTYPE)
    keys["run"] = np.asarray(run, dtype=np.uint64)
    keys["lumi"] = np.asarray(lumi, dtype=np.uint64)
    keys["event"] = np.asarray(event, dtype=np.uint64)
    return keys


def key_tuple(key: np.void) -> tuple[int, int, int]:
    """Convert a structured identity key to JSON-safe integers."""
    return int(key["run"]), int(key["lumi"]), int(key["event"])


def equal_numpy(left: np.ndarray, right: np.ndarray) -> bool:
    """Compare NumPy values exactly while treating matching NaNs as equal."""
    if left.dtype != right.dtype or left.shape != right.shape:
        return False
    if np.issubdtype(left.dtype, np.inexact):
        return bool(np.array_equal(left, right, equal_nan=True))
    return bool(np.array_equal(left, right))


def equal_awkward(left: ak.Array, right: ak.Array) -> bool:
    """Compare Awkward layouts exactly, including jagged structure and values."""
    left_form, left_length, left_buffers = ak.to_buffers(left)
    right_form, right_length, right_buffers = ak.to_buffers(right)
    if left_length != right_length or left_form.to_json() != right_form.to_json():
        return False
    if set(left_buffers) != set(right_buffers):
        return False
    return all(
        equal_numpy(np.asarray(left_buffers[name]), np.asarray(right_buffers[name]))
        for name in left_buffers
    )


def flat_changed_values(left: ak.Array, right: ak.Array) -> tuple[int, int]:
    """Return changed and total payload values for one jagged numeric branch."""
    left_values = np.asarray(ak.to_numpy(ak.flatten(left, axis=None)))
    right_values = np.asarray(ak.to_numpy(ak.flatten(right, axis=None)))
    if left_values.shape != right_values.shape:
        raise RuntimeError("Configured lepton collection structure changed unexpectedly")
    equal = left_values == right_values
    if np.issubdtype(left_values.dtype, np.inexact):
        equal |= np.isnan(left_values) & np.isnan(right_values)
    return int(np.count_nonzero(~equal)), int(len(left_values))


def require_existing(path: Path, label: str) -> None:
    """Fail clearly rather than accepting an absent prerequisite as a result."""
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")


def schema_missing(tree: uproot.behaviors.TTree.TTree) -> dict[str, list[str]]:
    """Report missing raw-object schema fields by binding group."""
    fields = set(tree.keys())
    return {
        group: [field for field in required if field not in fields]
        for group, required in REQUIRED_SCHEMA.items()
    }


def source_raw_from_h4l(source_h4l: Path) -> Path:
    """Map a versioned modified-H4l source path to its paired raw NanoAOD."""
    source = str(source_h4l)
    if "/h4l_mc_modified/" not in source or not source.endswith("_h4l.root"):
        raise RuntimeError(f"Cannot derive raw NanoAOD source from {source_h4l}")
    return Path(source.replace("/h4l_mc_modified/", "/nanoaod_modified/").removesuffix("_h4l.root") + ".root")


def h4l_sources_with_repairs() -> tuple[list[Path], dict[str, str]]:
    """Load mixer inputs and substitute only the two verified repair outputs."""
    require_existing(MIXER_MANIFEST, "Modified-H4l mixer manifest")
    manifest = json.loads(MIXER_MANIFEST.read_text())
    original_to_repaired = {
        str(spec["original_h4l"]): str(spec["repaired_h4l"])
        for spec in REPAIRS.values()
    }
    source_paths: list[Path] = []
    for sample in manifest["samples"]:
        source_paths.extend(Path(original_to_repaired.get(path, path)) for path in sample["files"])
    if len(source_paths) != 400 or len(set(source_paths)) != len(source_paths):
        raise RuntimeError("Mixer manifest does not provide 400 unique H4l source files")
    return source_paths, original_to_repaired


def first_matching_positions(
    sorted_targets: np.ndarray, tested: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return tested-row and sorted-target indices for exact structured-key matches."""
    positions = np.searchsorted(sorted_targets, tested)
    valid = positions < len(sorted_targets)
    comparison = np.minimum(positions, len(sorted_targets) - 1)
    valid &= sorted_targets[comparison] == tested
    rows = np.nonzero(valid)[0]
    return rows, positions[rows]


def compare_raw_pair(sample: str, spec: dict[str, Any]) -> dict[str, Any]:
    """Check output integrity, schema, identities, and allowed response changes."""
    original = spec["original_raw"]
    repaired = spec["repaired_raw"]
    require_existing(original, f"{sample} original NanoAOD")
    require_existing(repaired, f"{sample} regenerated NanoAOD")
    with uproot.open(original) as source_file, uproot.open(repaired) as repaired_file:
        if "Events" not in source_file or "Events" not in repaired_file:
            raise RuntimeError(f"{sample} is missing an Events tree after regeneration")
        source_tree = source_file["Events"]
        repaired_tree = repaired_file["Events"]
        source_branches = list(source_tree.keys())
        repaired_branches = list(repaired_tree.keys())
        if source_branches != repaired_branches:
            raise RuntimeError(f"{sample} regenerated Events branch schema differs from original")
        if source_tree.typenames() != repaired_tree.typenames():
            raise RuntimeError(f"{sample} regenerated Events branch types differ from original")
        if source_tree.num_entries != spec["entries"] or repaired_tree.num_entries != spec["entries"]:
            raise RuntimeError(
                f"{sample} event count is {source_tree.num_entries}/{repaired_tree.num_entries}; "
                f"expected {spec['entries']}"
            )
        missing = schema_missing(repaired_tree)
        if any(missing.values()):
            raise RuntimeError(f"{sample} regenerated NanoAOD is missing required schema: {missing}")

        source_ids = source_tree.arrays(IDENTITY_FIELDS, library="np")
        repaired_ids = repaired_tree.arrays(IDENTITY_FIELDS, library="np")
        source_keys = event_keys(
            source_ids["run"], source_ids["luminosityBlock"], source_ids["event"]
        )
        repaired_keys = event_keys(
            repaired_ids["run"], repaired_ids["luminosityBlock"], repaired_ids["event"]
        )
        source_unique = len(np.unique(source_keys))
        repaired_unique = len(np.unique(repaired_keys))
        if source_unique != len(source_keys) or repaired_unique != len(repaired_keys):
            raise RuntimeError(f"{sample} original-to-repaired identity mapping is not unique")
        if not equal_numpy(source_keys, repaired_keys):
            raise RuntimeError(f"{sample} original-to-repaired identity ordering/content differs")

        source_classes = source_file.classnames(recursive=True)
        repaired_classes = repaired_file.classnames(recursive=True)
        if source_classes != repaired_classes:
            raise RuntimeError(f"{sample} non-Events object schema differs after regeneration")

        nonlepton_fields = [
            field
            for field in source_branches
            if field not in EXPECTED_MUTATED_FIELDS
            and field not in {"nMuon", "nElectron"}
            and not field.startswith("Muon_")
            and not field.startswith("Electron_")
        ]
        unchanged_lepton_fields = [
            field
            for field in source_branches
            if field not in EXPECTED_MUTATED_FIELDS and field not in nonlepton_fields
        ]
        nonlepton_exact = True
        unchanged_lepton_exact = True
        mismatch_fields: set[str] = set()
        changed_values = {field: 0 for field in sorted(EXPECTED_MUTATED_FIELDS)}
        total_values = {field: 0 for field in sorted(EXPECTED_MUTATED_FIELDS)}
        chunk_size = 50_000
        for start in range(0, source_tree.num_entries, chunk_size):
            stop = min(start + chunk_size, source_tree.num_entries)
            source_arrays = source_tree.arrays(
                source_branches, entry_start=start, entry_stop=stop, library="ak"
            )
            repaired_arrays = repaired_tree.arrays(
                source_branches, entry_start=start, entry_stop=stop, library="ak"
            )
            if not equal_awkward(source_arrays[nonlepton_fields], repaired_arrays[nonlepton_fields]):
                nonlepton_exact = False
                mismatch_fields.update(
                    field
                    for field in nonlepton_fields
                    if not equal_awkward(source_arrays[field], repaired_arrays[field])
                )
            if not equal_awkward(
                source_arrays[unchanged_lepton_fields], repaired_arrays[unchanged_lepton_fields]
            ):
                unchanged_lepton_exact = False
                mismatch_fields.update(
                    field
                    for field in unchanged_lepton_fields
                    if not equal_awkward(source_arrays[field], repaired_arrays[field])
                )
            for field in EXPECTED_MUTATED_FIELDS:
                changed, total = flat_changed_values(source_arrays[field], repaired_arrays[field])
                changed_values[field] += changed
                total_values[field] += total
        if not nonlepton_exact or not unchanged_lepton_exact or mismatch_fields:
            raise RuntimeError(
                f"{sample} has unexpected unchanged-content differences in "
                f"{sorted(mismatch_fields)}"
            )
        if any(changed_values[field] == 0 for field in EXPECTED_MUTATED_FIELDS):
            raise RuntimeError(f"{sample} did not change every configured pt/eta/phi response field")

    return {
        "original_raw": str(original),
        "repaired_raw": str(repaired),
        "events_entries": int(spec["entries"]),
        "events_tree_readable": True,
        "original_to_repaired_identity_mapping": {
            "complete": True,
            "unique": True,
            "source_unique": int(source_unique),
            "repaired_unique": int(repaired_unique),
        },
        "required_schema_missing": missing,
        "events_schema_and_types_exact": True,
        "non_events_object_schema_exact": True,
        "nonlepton_event_content_exact": True,
        "unchanged_lepton_content_exact": True,
        "configured_response_fields": {
            field: {"changed_values": changed_values[field], "total_values": total_values[field]}
            for field in sorted(EXPECTED_MUTATED_FIELDS)
        },
        "all_jet_branches_present": sorted(
            field for field in repaired_branches if field.startswith("Jet_")
        ),
    }


def check_response_sidecars(sample: str, spec: dict[str, Any]) -> dict[str, Any]:
    """Require the response summary, histogram sidecar, and successful task record."""
    config_path = spec["response_config"]
    require_existing(config_path, f"{sample} repair config")
    config = json.loads(config_path.read_text())
    summary_dir = Path(config["validation"]["directory"])
    summary_path = summary_dir / "lepton_response_summary.json"
    require_existing(summary_path, f"{sample} response validation summary")
    summary = json.loads(summary_path.read_text())
    if int(summary["n_events"]) != spec["entries"]:
        raise RuntimeError(f"{sample} response summary has an unexpected event count")
    if summary["outputs"] != [
        {
            "sample": sample,
            "input": str(spec["original_raw"]),
            "output": str(spec["repaired_raw"]),
        }
    ]:
        raise RuntimeError(f"{sample} response summary does not bind the expected input/output")
    histogram_path = Path(summary["validation_histograms"]["path"])
    require_existing(histogram_path, f"{sample} response validation histogram sidecar")
    task_records = sorted(summary_dir.rglob("task_execution_*.txt"))
    if not task_records:
        raise RuntimeError(f"{sample} has no response worker execution record")
    records = [record.read_text() for record in task_records]
    if not any("exit_code=0" in record for record in records):
        raise RuntimeError(f"{sample} response worker did not record successful completion")
    failure_text = "\n".join(
        path.read_text(errors="replace")
        for pattern in ("*.out", "*.err")
        for path in CONDOR_ROOT.rglob(pattern)
    )
    forbidden = ["TFile::Flush", "No such device", "Unable to create output ROOT file"]
    found = [message for message in forbidden if message in failure_text]
    if found:
        raise RuntimeError(f"{sample} Condor logs contain ROOT write failures: {found}")
    return {
        "config": str(config_path),
        "validation_summary": str(summary_path),
        "validation_histogram_sidecar": str(histogram_path),
        "worker_execution_records": [str(record) for record in task_records],
        "write_close_failure_strings_absent": True,
    }


def equal_h4l_rows(
    fake_arrays: dict[str, np.ndarray], fake_index: int, source_arrays: dict[str, np.ndarray]
) -> bool:
    """Apply the ticket's all-111-copied-branch exact payload rule."""
    for field, fake_values in fake_arrays.items():
        if not equal_numpy(
            np.asarray(fake_values[fake_index : fake_index + 1]),
            np.asarray(source_arrays[field]),
        ):
            return False
    return True


def compare_h4l_pair(sample: str, spec: dict[str, Any]) -> dict[str, Any]:
    """Require every regenerated H4l source entry to match the current payload."""
    original = spec["original_h4l"]
    repaired = spec["repaired_h4l"]
    require_existing(original, f"{sample} current modified H4l")
    require_existing(repaired, f"{sample} regenerated modified H4l")
    with uproot.open(original) as source_file, uproot.open(repaired) as repaired_file:
        if "h4lTree" not in source_file or "h4lTree" not in repaired_file:
            raise RuntimeError(f"{sample} H4l output lacks h4lTree")
        source_tree = source_file["h4lTree"]
        repaired_tree = repaired_file["h4lTree"]
        fields = list(source_tree.keys())
        if fields != list(repaired_tree.keys()) or len(fields) != 111:
            raise RuntimeError(f"{sample} regenerated H4l schema differs from the 111 copied fields")
        if source_tree.num_entries != repaired_tree.num_entries:
            raise RuntimeError(f"{sample} regenerated H4l entry count differs from current payload")
        source = source_tree.arrays(fields, library="np")
        repaired_arrays = repaired_tree.arrays(fields, library="np")
        mismatches = [field for field in fields if not equal_numpy(source[field], repaired_arrays[field])]
        if mismatches:
            raise RuntimeError(f"{sample} regenerated H4l payload differs in {mismatches}")
    return {
        "current_h4l": str(original),
        "regenerated_h4l": str(repaired),
        "h4l_entries": int(source_tree.num_entries),
        "h4l_fields": len(fields),
        "all_source_entries_exact": True,
    }


def resolve_fake_sources(
    h4l_sources: list[Path], original_to_repaired: dict[str, str]
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    """Resolve every fake row by identity plus exact payload, preserving source tokens."""
    require_existing(FAKE_DATA, "Synthetic pseudo-data H4l file")
    with uproot.open(FAKE_DATA) as fake_file:
        fake_tree = fake_file["h4lTree"]
        fields = list(fake_tree.keys())
        if len(fields) != 111:
            raise RuntimeError("Fake H4l tree no longer has the expected 111 copied fields")
        fake_arrays = fake_tree.arrays(fields, library="np")
    fake_keys = event_keys(fake_arrays["run"], fake_arrays["lumi"], fake_arrays["event"])
    if len(fake_keys) != 652 or len(np.unique(fake_keys)) != len(fake_keys):
        raise RuntimeError("Fake H4l identity set is not the expected 652 unique rows")
    sorted_indices = np.argsort(fake_keys)
    sorted_fake_keys = fake_keys[sorted_indices]
    candidates: dict[int, list[tuple[Path, int]]] = defaultdict(list)
    for source_number, source_path in enumerate(h4l_sources, start=1):
        require_existing(source_path, "Modified-H4l mixer source")
        with uproot.open(source_path) as source_file:
            tree = source_file["h4lTree"]
            if list(tree.keys()) != fields:
                raise RuntimeError(f"H4l source schema differs from fake schema: {source_path}")
            entry_start = 0
            for ids in tree.iterate(FAKE_IDENTITY_FIELDS, library="np", step_size="50 MB"):
                source_keys = event_keys(ids["run"], ids["lumi"], ids["event"])
                rows, positions = first_matching_positions(sorted_fake_keys, source_keys)
                for row, position in zip(rows, positions, strict=True):
                    candidates[int(sorted_indices[position])].append(
                        (source_path, entry_start + int(row))
                    )
                entry_start += len(source_keys)
        if source_number % 50 == 0:
            log.info("Resolved identity candidates from %d/%d H4l mixer files.", source_number, len(h4l_sources))

    mapping: list[dict[str, Any]] = []
    for fake_index in range(len(fake_keys)):
        exact: list[tuple[Path, int]] = []
        for source_path, source_entry in candidates[fake_index]:
            with uproot.open(source_path) as source_file:
                source_values = source_file["h4lTree"].arrays(
                    fields, entry_start=source_entry, entry_stop=source_entry + 1, library="np"
                )
            if equal_h4l_rows(fake_arrays, fake_index, source_values):
                exact.append((source_path, source_entry))
        if len(exact) != 1:
            raise RuntimeError(
                f"Fake row {key_tuple(fake_keys[fake_index])} has {len(exact)} exact full-payload source matches"
            )
        source_path, source_entry = exact[0]
        prior_source = next(
            (original for original, repaired in original_to_repaired.items() if repaired == str(source_path)),
            str(source_path),
        )
        mapping.append(
            {
                "fake_row": fake_index,
                "run": int(fake_keys[fake_index]["run"]),
                "lumi": int(fake_keys[fake_index]["lumi"]),
                "event": int(fake_keys[fake_index]["event"]),
                "source_h4l_file": str(source_path),
                "prior_source_h4l_file": prior_source,
                "source_h4l_entry": int(source_entry),
                "raw_file": str(source_raw_from_h4l(source_path)),
                "raw_entry": None,
            }
        )
    return fake_arrays, mapping


def repeat_complete_raw_scan(mapping: list[dict[str, Any]]) -> dict[str, Any]:
    """Repeat the complete 202-file raw scan with the two repaired counterparts."""
    by_raw: dict[Path, list[int]] = defaultdict(list)
    for fake_index, row in enumerate(mapping):
        by_raw[Path(row["raw_file"])].append(fake_index)
    if len(by_raw) != 202:
        raise RuntimeError(f"Expected 202 resolved raw sources, found {len(by_raw)}")

    match_counts = np.zeros(len(mapping), dtype=np.int64)
    schema_failures: dict[str, dict[str, list[str]]] = {}
    unreadable: list[str] = []
    entries_scanned = 0
    for file_number, (raw_path, indices) in enumerate(sorted(by_raw.items()), start=1):
        try:
            with uproot.open(raw_path) as source_file:
                if "Events" not in source_file:
                    raise RuntimeError("no Events tree")
                tree = source_file["Events"]
                missing = schema_missing(tree)
                if any(missing.values()):
                    schema_failures[str(raw_path)] = missing
                expected_keys = np.empty(len(indices), dtype=KEY_DTYPE)
                expected_keys["run"] = [mapping[index]["run"] for index in indices]
                expected_keys["lumi"] = [mapping[index]["lumi"] for index in indices]
                expected_keys["event"] = [mapping[index]["event"] for index in indices]
                sorted_order = np.argsort(expected_keys)
                sorted_expected = expected_keys[sorted_order]
                entry_start = 0
                for arrays in tree.iterate(IDENTITY_FIELDS, library="np", step_size="50 MB"):
                    raw_keys = event_keys(
                        arrays["run"], arrays["luminosityBlock"], arrays["event"]
                    )
                    rows, positions = first_matching_positions(sorted_expected, raw_keys)
                    for row, position in zip(rows, positions, strict=True):
                        mapping_index = indices[int(sorted_order[position])]
                        match_counts[mapping_index] += 1
                        if match_counts[mapping_index] == 1:
                            mapping[mapping_index]["raw_entry"] = entry_start + int(row)
                    entry_start += len(raw_keys)
                    entries_scanned += len(raw_keys)
        except Exception as error:  # noqa: BLE001 - report every selected raw source.
            unreadable.append(f"{raw_path}: {error}")
        if file_number % 25 == 0 or file_number == len(by_raw):
            log.info("Completed raw identity/schema scan for %d/%d sources.", file_number, len(by_raw))

    missing_matches = int(np.count_nonzero(match_counts == 0))
    duplicate_matches = int(np.count_nonzero(match_counts > 1))
    if unreadable or schema_failures or missing_matches or duplicate_matches:
        raise RuntimeError(
            "Raw pseudo-data provenance scan failed: "
            f"unreadable={len(unreadable)}, schema_failures={len(schema_failures)}, "
            f"missing={missing_matches}, duplicates={duplicate_matches}"
        )
    return {
        "resolved_raw_files": len(by_raw),
        "raw_entries_scanned": int(entries_scanned),
        "exact_matches": int(np.count_nonzero(match_counts == 1)),
        "missing_matches": missing_matches,
        "duplicate_matches": duplicate_matches,
        "unreadable_raw_files": unreadable,
        "schema_failures": schema_failures,
    }


def main() -> None:
    """Run every ticketed acceptance check and persist results without replacement."""
    if OUTPUT.exists() or MAPPING_OUTPUT.exists():
        raise FileExistsError("Regeneration evidence already exists; create a new version for a rerun")

    raw_checks: dict[str, dict[str, Any]] = {}
    response_checks: dict[str, dict[str, Any]] = {}
    h4l_checks: dict[str, dict[str, Any]] = {}
    for sample, spec in REPAIRS.items():
        response_checks[sample] = check_response_sidecars(sample, spec)
        raw_checks[sample] = compare_raw_pair(sample, spec)
        h4l_checks[sample] = compare_h4l_pair(sample, spec)
        log.info("Passed response, raw-content, and H4l checks for %s.", sample)

    h4l_sources, original_to_repaired = h4l_sources_with_repairs()
    fake_arrays, mapping = resolve_fake_sources(h4l_sources, original_to_repaired)
    affected = {
        identity
        for spec in REPAIRS.values()
        for identity in spec["affected_fake_identities"]
    }
    resolved_affected = {
        (row["run"], row["lumi"], row["event"])
        for row in mapping
        if row["source_h4l_file"].startswith(str(REPAIR_ROOT))
    }
    if resolved_affected != affected:
        raise RuntimeError("The five affected fake rows did not resolve to the fresh H4l sources")
    raw_scan = repeat_complete_raw_scan(mapping)

    mapping_payload = {
        "schema_version": 1,
        "purpose": "Persistent source-file/source-entry and raw-entry provenance for repaired modified-MC pseudo-data; not collision-data provenance.",
        "fake_data": str(FAKE_DATA),
        "fake_rows": len(mapping),
        "h4l_payload_fields_compared": len(fake_arrays),
        "mapping": mapping,
    }
    result = {
        "schema_version": 1,
        "generated_utc": datetime.now(UTC).isoformat(),
        "status": "pass",
        "scope": "Repair validation for synthetic modified-MC pseudo-data only. This does not supply collision-data NanoAOD provenance or authorize a collision-data VBF/FSR fit.",
        "versioned_repair_root": str(REPAIR_ROOT),
        "response_validation": response_checks,
        "regenerated_raw_validation": raw_checks,
        "regenerated_h4l_validation": h4l_checks,
        "pseudo_data_h4l_payload_join": {
            "fake_rows": len(mapping),
            "exact_full_payload_matches": len(mapping),
            "payload_fields_compared": len(fake_arrays),
            "affected_rows_resolved_to_regenerated_h4l": len(resolved_affected),
            "source_token_manifest": str(MAPPING_OUTPUT),
        },
        "complete_raw_identity_schema_scan": raw_scan,
        "collision_data_status": "unresolved; repaired sources are modified-MC pseudo-data only",
    }
    MAPPING_OUTPUT.write_text(json.dumps(mapping_payload, indent=2) + "\n")
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    log.info("Wrote %s and %s.", OUTPUT, MAPPING_OUTPUT)


if __name__ == "__main__":
    main()
