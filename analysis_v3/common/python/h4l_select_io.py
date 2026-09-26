"""Readers of the final-selection outputs (program h4l_select) for the downstream stages.

The only entry point to the outputs is the full scan
production_v3/h4l_select/<version>/scan.json (never a directory listing).  MC
rows are normalized with the genWeight convention of AGENTS.md:
    w = weight (genWeight) x sigma_eff x L / sum(genEventSumw)
with sigma_eff of the manifest, sum(genEventSumw) over exactly the files the
scan covers and L the luminosity of the data sample (given explicitly for a
scan without data, the N-1 selections); data rows have w = 1.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import uproot

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"


def load_scan(version: str) -> dict:
    scan = json.loads((PRODUCTION / "h4l_select" / version / "scan.json").read_text(encoding="utf-8"))
    if scan.get("schema") != "h4l_v3_select_scan/1":
        raise RuntimeError(f"unexpected scan schema {scan.get('schema')}")
    return scan


def lumi_fb(scan: dict) -> float:
    lumis = {s["lumi_fb"] for s in scan["samples"].values() if s["kind"] == "data"}
    if len(lumis) != 1:
        raise RuntimeError(f"expected one data luminosity, found {lumis}")
    return lumis.pop()


def weight_scale(scan: dict, sample: str, lumi: float | None = None) -> float:
    """sigma_eff x L / sum(genEventSumw) (MC), 1 (data); L from the scan's data sample unless given (MC-only scans)."""
    info = scan["samples"][sample]
    if info["kind"] == "data":
        return 1.0
    if lumi is None:
        lumi = lumi_fb(scan)
    elif any(s["kind"] == "data" for s in scan["samples"].values()) and abs(lumi - lumi_fb(scan)) > 1e-9 * lumi:
        raise RuntimeError(f"the given luminosity {lumi} differs from the scan's data luminosity {lumi_fb(scan)}")
    return info["sigma_eff_pb"] * 1000.0 * lumi / info["genEventSumw"]


def read(scan: dict, samples: list[str], tree: str, branches: list[str], cut=None, lumi: float | None = None) -> dict:
    """Concatenate the branches of one tree over the outputs of the samples.

    Adds "w" (the normalized weight), "sample" (index into samples) and keeps the
    row order of the scan.  cut(arrays) -> boolean mask is applied per output.  lumi (fb^-1) is required for a
    scan without a data sample (the N-1 selections)."""
    wanted = sorted(set(branches) | {"weight"})
    parts: dict[str, list] = {name: [] for name in wanted + ["w", "sample"]}
    for index, sample in enumerate(samples):
        if sample not in scan["samples"]:
            raise KeyError(f"sample {sample} not in the scan")
        scale = weight_scale(scan, sample, lumi)
        for output in scan["outputs"]:
            if output["sample"] != sample or output["trees"].get(tree, 0) == 0:
                continue
            arrays = uproot.open(output["path"])[tree].arrays(wanted, library="np")
            mask = np.ones(len(arrays["weight"]), dtype=bool) if cut is None else cut(arrays)
            for name in wanted:
                parts[name].append(arrays[name][mask])
            weight = arrays["weight"][mask].astype(np.float64)
            parts["w"].append(weight * scale if scan["samples"][sample]["kind"] == "mc" else np.ones_like(weight))
            parts["sample"].append(np.full(int(mask.sum()), index, dtype=np.int16))
    out = {}
    for name, chunks in parts.items():
        if chunks:
            out[name] = np.concatenate(chunks)
        else:
            out[name] = np.zeros(0)
    return out
