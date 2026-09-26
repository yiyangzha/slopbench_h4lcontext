"""Validation plots of the v3 input manifests (MC bookkeeping only).

    pixi run py -- analysis_v3/framework/scripts/plot_manifests.py --version v1

Writes production_v3/manifests/<version>/plots/: valid files per sample,
selected events per sample, the genWeight- to count-based normalization
ratio per sample, and the events-per-file distribution of every MC sample.
The pseudo-data enter only as their number of valid shards (no event counts
are plotted, as the blinding rules require).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--plot-dir", default="plots", help="output directory name below the manifest version (never reused)")
    args = parser.parse_args()
    base = REPO / "production_v3" / "manifests" / args.version
    summary = json.loads((base / "summary.json").read_text(encoding="utf-8"))
    plots = base / args.plot_dir
    if plots.exists():
        raise SystemExit(f"{plots} exists; choose a new --plot-dir")
    mc = {name: entry for name, entry in summary["samples"].items() if entry["kind"] == "mc"}
    names = list(mc)
    positions = np.arange(len(names))
    written = []

    fig, axis = h4l_style.figure()
    axis.bar(positions, [mc[name]["valid"] for name in names], color="#5790fc", label="valid files")
    axis.bar(positions, [mc[name]["files"] - mc[name]["valid"] for name in names],
             bottom=[mc[name]["valid"] for name in names], color="#e42536", label="invalid or unreadable")
    axis.set_yscale("log")
    top = max(mc[name]["files"] for name in names)
    axis.set_ylim(0.5, top * 30)
    for position, name in zip(positions, names):
        axis.text(position, mc[name]["files"] * 1.15, f"{mc[name]['valid']}/{mc[name]['files']}", ha="center", va="bottom")
    axis.set_ylabel("files (valid/listed)")
    axis.set_xticks(positions, names, rotation=60, ha="right")
    axis.legend(loc="upper right")
    h4l_style.label(axis, data=False)
    written += h4l_style.save(fig, plots / "files_per_sample")

    fig, axis = h4l_style.figure()
    axis.bar(positions, [mc[name]["events_selected"] for name in names], color="#f89c20")
    axis.set_yscale("log")
    axis.set_ylim(0.5, max(mc[name]["events_selected"] for name in names) * 30)
    for position, name in zip(positions, names):
        axis.text(position, mc[name]["events_selected"] * 1.15, f"{mc[name]['events_selected']:.3g}", ha="center",
                  va="bottom", fontsize="small")
    axis.set_ylabel("preselected MC events")
    axis.set_xticks(positions, names, rotation=60, ha="right")
    h4l_style.label(axis, data=False)
    written += h4l_style.save(fig, plots / "events_per_sample")

    fig, axis = h4l_style.figure()
    axis.plot(positions, [mc[name]["ratio_weight_to_count"] for name in names], "o", color="black")
    axis.axhline(1.0, color="grey", linestyle="--")
    axis.set_ylabel("preselected fraction: genWeight-based / count-based")
    axis.set_xticks(positions, names, rotation=60, ha="right")
    h4l_style.label(axis, data=False)
    written += h4l_style.save(fig, plots / "normalization_ratio")

    for name in names:
        manifest = json.loads((base / f"mc_{name}.json").read_text(encoding="utf-8"))
        entries = np.array([item["events_entries"] for item in manifest["files"]])
        fig, axis = h4l_style.figure()
        axis.hist(entries, bins=min(50, max(5, len(entries) // 5)), histtype="step", color="black", linewidth=2)
        axis.set_xlabel(f"{name}: preselected events per file")
        axis.set_ylabel("files")
        axis.text(0.95, 0.9, f"{len(entries)} files, {entries.sum()} events", transform=axis.transAxes, ha="right")
        h4l_style.label(axis, data=False)
        written += h4l_style.save(fig, plots / f"events_per_file_{name}")

    for path in written:
        print(f"[plot] {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
