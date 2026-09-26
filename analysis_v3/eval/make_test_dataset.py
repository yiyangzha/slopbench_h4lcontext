"""A test dataset in the DATASET.md layout of the evaluation task, built from the analysis inputs by symbolic links
(nothing is copied or modified): data/part_*.root = a subset of the UL16 pseudo-data shards, mc/<process>/*.root = the
h4l_seeds_v2 MC files of the validated manifests (block subdirectories flattened into the file names),
mc/<process>/meta.json, mc/cross_sections.json (the effective cross sections of prompt_renew.md 4.2), lumi.json and
meta.json.

    pixi run py -- analysis_v3/eval/make_test_dataset.py --name test_10fb --data-shards 41 [--mc-fraction 1.0]

The luminosity of a shard subset is the full luminosity times its fraction of the pseudo-data entries (the shards are
equal slices of one dataset).  Writes production_v3/eval_datasets/<name>/ (never overwritten).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
XSEC = {"ZZTo4L": 1.325, "GGZZ4Mu": 0.001575, "GGZZ4E": 0.001619, "GGZZ2E2Mu": 0.003185, "DYJetsToLL": 5396.0, "TTBar": 87.58,
        "GluGluToHToZZ_M125": 0.006024, "VBF_HToZZ_M125": 0.00048794, "VHToZZ_M125": 0.0002726892}
HLT = json.loads((REPO / "analysis_v3/common/config/analysis_ul16_v3.json").read_text(encoding="utf-8"))["triggers"]["analysis_or"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", required=True)
    parser.add_argument("--data-shards", type=int, required=True)
    parser.add_argument("--mc-fraction", type=float, default=1.0, help="fraction of the files of each MC process (every k-th)")
    parser.add_argument("--processes", nargs="*", default=list(XSEC))
    args = parser.parse_args()
    out = PRODUCTION / "eval_datasets" / args.name
    if out.exists():
        raise SystemExit(f"{out} exists")
    data = json.loads((PRODUCTION / "manifests/v2/data.json").read_text(encoding="utf-8"))
    files = sorted(data["files"], key=lambda f: f["relative"])
    total_entries = sum(f["events_entries"] for f in files)
    chosen = files[: args.data_shards]
    lumi = data["lumi_fb"] * sum(f["events_entries"] for f in chosen) / total_entries
    (out / "data").mkdir(parents=True)
    for f in chosen:
        os.symlink(f["path"], out / "data" / f["relative"])
    xsec = {}
    for name in args.processes:
        manifest = json.loads((PRODUCTION / f"manifests/v2/mc_{name}.json").read_text(encoding="utf-8"))
        mfiles = sorted(manifest["files"], key=lambda f: f["relative"])
        step = max(1, round(1.0 / args.mc_fraction))
        mfiles = mfiles[::step]
        d = out / "mc" / name
        d.mkdir(parents=True)
        for f in mfiles:
            os.symlink(f["path"], d / f["relative"].replace("/", "_"))
        meta = {"process": name, "n_files": len(mfiles), "entries": sum(f["events_entries"] for f in mfiles),
                "entries_before_selection": sum(f["entries_before_selection"] for f in mfiles),
                "n_preselection": sum(f["entries_before_selection"] for f in mfiles)}
        (d / "meta.json").write_text(json.dumps(meta, indent=1) + "\n", encoding="utf-8")
        xsec[name] = XSEC[name]
    (out / "mc" / "cross_sections.json").write_text(json.dumps(xsec, indent=1) + "\n", encoding="utf-8")
    (out / "lumi.json").write_text(json.dumps({"lumi_fb": lumi}) + "\n", encoding="utf-8")
    (out / "meta.json").write_text(json.dumps({"n_files": len(chosen), "n_events": sum(f["events_entries"] for f in chosen),
                                               "hlt_paths": HLT}, indent=1) + "\n", encoding="utf-8")
    print(f"[dataset] {out}: {len(chosen)} data shards, lumi {lumi:.4f} fb^-1, processes {args.processes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
