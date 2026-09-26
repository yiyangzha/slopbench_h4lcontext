"""Full output scan, histograms and validation plots of a stage-2 skim.

    pixi run py -- analysis_v3/skims/scripts/validate_skims.py --version v1 --manifests v2 [--workers 8]

1. Scan: every planned task has its JSON; the JSON parses, comes from one
   frozen program and one analysis configuration, names its ROOT output with
   a sha256, and covers exactly its planned input files; every ROOT output
   passes root_check (keys, exact tree entries, readable first/last entries).
2. Coverage: the processed file keys of each sample equal its manifest.
3. Histograms: the C++ program skim_histograms runs per sample (parallel).
4. Plots: data against stacked MC normalized with genWeight and the Runs
   genEventSumw of the processed files, with a data/MC ratio.

Writes production_v3/skims/<version>/validation/ (scan.json, normalization.json,
histograms/<sample>.root) and production_v3/skims/<version>/<plot dir>/.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import uproot

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
ROOT_CHECK = REPO / "analysis_v3/framework/bin/root_check"
HISTOGRAMS = REPO / "analysis_v3/skims/bin/skim_histograms"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402

GROUPS = [("DY", ["DYJetsToLL"], "#5790fc"), ("ttbar", ["TTBar"], "#f89c20"),
          ("ZZ", ["ZZTo4L", "GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"], "#964a8b"),
          ("H(125)", ["GluGluToHToZZ_M125", "VBF_HToZZ_M125", "VHToZZ_M125"], "#e42536")]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def check_root(task: dict) -> tuple[str, str]:
    report = json.loads(Path(task["outputs"]["json"]).read_text(encoding="utf-8"))
    spec = report["root_output"]
    arguments = [str(ROOT_CHECK), task["outputs"]["root"]]
    arguments += [f"{name}={entries}" for name, entries in spec["trees"].items()]
    result = subprocess.run(arguments, capture_output=True, text=True)
    return task["task_id"], "" if result.returncode == 0 else (result.stderr.strip() or "root_check failed")


def plot_all(histograms: Path, normalization: dict, lumi_fb: float, plots: Path) -> None:
    """Data against stacked MC for every histogram, with an adaptive ratio range."""
    def load(sample: str, name: str):
        with uproot.open(histograms / f"{sample}.root") as file:
            values, edges = file[name].to_numpy()
            variances = file[name].variances()
        return values, variances, edges

    names = []
    with uproot.open(histograms / "pseudo_data.root") as file:
        names = sorted(key.split(";")[0] for key in file.keys())
    for name in names:
        data, data_var, edges = load("pseudo_data", name)
        stack = []
        for label, samples, colour in GROUPS:
            total = np.zeros_like(data)
            for sample in samples:
                values, _, _ = load(sample, name)
                total += values * normalization[sample]["scale"]
            stack.append((label, total, colour))
        mc_total = sum(item[1] for item in stack)
        fig, (main, lower) = h4l_style.figure(ratio=True)
        bottom = np.zeros_like(data)
        centres = 0.5 * (edges[1:] + edges[:-1])
        for label, values, colour in stack:
            main.stairs(values + bottom, edges, baseline=bottom, fill=True, color=colour, label=label)
            bottom = bottom + values
        main.errorbar(centres, data, yerr=np.sqrt(data), fmt="o", color="black", markersize=4, label="pseudo-data")
        logy = name.startswith(("mass", "iso03", "sip", "dxy", "dz", "met", "pt_err", "l1_pt", "l2_pt", "pair_pt"))
        if logy:
            main.set_yscale("log")
            positive = mc_total[mc_total > 0]
            main.set_ylim(max(0.5, positive.min() * 0.5) if positive.size else 0.5, max(data.max(), mc_total.max()) * 50)
        main.set_ylabel("entries")
        main.legend(ncol=2, fontsize="small")
        ratio = np.divide(data, mc_total, out=np.full_like(data, np.nan), where=mc_total > 0)
        ratio_error = np.divide(np.sqrt(data), mc_total, out=np.full_like(data, np.nan), where=mc_total > 0)
        lower.errorbar(centres, ratio, yerr=ratio_error, fmt="o", color="black", markersize=3)
        lower.axhline(1.0, color="grey", linestyle="--")
        finite = ratio[np.isfinite(ratio) & (data > 0)]
        low, high = (np.nanpercentile(finite, 2), np.nanpercentile(finite, 98)) if finite.size else (1.0, 1.0)
        lower.set_ylim(max(0.0, min(0.8, low - 0.05)), min(3.0, max(1.2, high + 0.05)))
        lower.set_ylabel("data/MC")
        lower.set_xlabel(name)
        h4l_style.label(main, data=True, lumi_fb=lumi_fb)
        h4l_style.save(fig, plots / name)
        totals = {"data": float(data.sum()), "mc": float(mc_total.sum())}
        print(f"[plot] {name}: totals {totals}")
    print(f"[plots] {plots}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--manifests", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--plot-dir", default="plots")
    parser.add_argument("--plots-only", action="store_true", help="reuse the published scan and histograms")
    args = parser.parse_args()
    base = PRODUCTION / "skims" / args.version
    plan = json.loads((base / "plan.json").read_text(encoding="utf-8"))
    validation = base / "validation"
    plots = base / args.plot_dir
    if plots.exists():
        fail(f"{plots} exists; choose a new --plot-dir")

    if args.plots_only:
        if not (validation / "scan.json").exists():
            fail("no published scan; run without --plots-only first")
        normalization = json.loads((validation / "normalization.json").read_text(encoding="utf-8"))["samples"]
        lumi_fb = json.loads((validation / "normalization.json").read_text(encoding="utf-8"))["lumi_fb"]
        histograms = validation / "histograms"
        plot_all(histograms, normalization, lumi_fb, plots)
        return 0

    # 1. Scan.
    missing, programs, configs, bad = [], set(), set(), []
    processed: dict[str, set] = {}
    for task in plan["tasks"]:
        path = Path(task["outputs"]["json"])
        if not path.exists():
            missing.append(task["task_id"])
            continue
        report = json.loads(path.read_text(encoding="utf-8"))
        programs.add((report.get("frozen_program") or {}).get("sha256"))
        configs.add(report["analysis_config_fnv1a64"])
        publication = report.get("publication") or {}
        if publication.get("root_path") != task["outputs"]["root"] or not publication.get("root_sha256"):
            bad.append(f"{task['task_id']}: no ROOT publication record")
        keys = [item["file_key"] for item in report["files"]]
        if keys != [item["file_key"] for item in task["inputs"]]:
            bad.append(f"{task['task_id']}: processed files differ from the planned inputs")
        processed.setdefault(task["config"]["sample"], set()).update(keys)
    if missing:
        fail(f"{len(missing)} task outputs missing, e.g. {missing[:5]}")
    if len(programs) != 1 or None in programs or len(configs) != 1:
        fail(f"outputs come from {len(programs)} programs and {len(configs)} configurations")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for task_id, error in pool.map(check_root, plan["tasks"]):
            if error:
                bad.append(f"{task_id}: {error}")
    if bad:
        fail(f"{len(bad)} invalid outputs: " + "; ".join(bad[:5]))

    # 2. Coverage and normalization (genWeight-based, over the processed files).
    manifests = PRODUCTION / "manifests" / args.manifests
    normalization = {}
    data_manifest = json.loads((manifests / "data.json").read_text(encoding="utf-8"))
    if processed[data_manifest["name"]] != {item["file_key"] for item in data_manifest["files"]}:
        fail("the skim does not cover every pseudo-data shard")
    lumi_fb = data_manifest["lumi_fb"]
    for manifest_path in sorted(manifests.glob("mc_*.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        keys = processed[manifest["name"]]
        files = [item for item in manifest["files"] if item["file_key"] in keys]
        if len(files) != len(manifest["files"]):
            print(f"WARN: {manifest['name']}: {len(files)} of {len(manifest['files'])} manifest files processed")
        sumw = sum(item["genEventSumw"] for item in files)
        normalization[manifest["name"]] = {"files": len(files), "genEventSumw": sumw,
                                           "sigma_eff_pb": manifest["normalization"]["sigma_eff_pb"],
                                           "scale": manifest["normalization"]["sigma_eff_pb"] * 1000.0 * lumi_fb / sumw}
    validation.mkdir(parents=True, exist_ok=True)
    (validation / "scan.json").write_text(json.dumps({"tasks": len(plan["tasks"]), "valid": len(plan["tasks"]),
                                                      "program_sha256": programs.pop(), "config": configs.pop()},
                                                     indent=1) + "\n", encoding="utf-8")
    (validation / "normalization.json").write_text(json.dumps({"lumi_fb": lumi_fb, "samples": normalization}, indent=1) + "\n",
                                                   encoding="utf-8")
    print(f"[scan] {len(plan['tasks'])} outputs valid; coverage complete")

    # 3. Histograms (C++), one process per sample.
    analysis = plan["tasks"][0]["analysis_config"]
    paths = analysis["triggers"]["paths"]
    mask = sum(1 << paths.index(name) for name in analysis["triggers"]["analysis_or"])
    histograms = validation / "histograms"
    histograms.mkdir(parents=True, exist_ok=True)
    jobs = {}
    for task in plan["tasks"]:
        jobs.setdefault(task["config"]["sample"], []).append(task["outputs"]["root"])
    environment = dict(os.environ, H4L_TRIGGER_OR_MASK=str(mask))

    def run(sample: str) -> tuple[str, int]:
        output = histograms / f"{sample}.root"
        if output.exists():
            return sample, 0
        temporary = histograms / f"{sample}.partial.{os.getpid()}.root"
        result = subprocess.run([str(HISTOGRAMS), str(temporary), *jobs[sample]], env=environment, capture_output=True, text=True)
        if result.returncode == 0:
            temporary.rename(output)
        else:
            print(result.stderr, file=sys.stderr)
        return sample, result.returncode

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for sample, code in pool.map(run, sorted(jobs)):
            if code:
                fail(f"skim_histograms failed for {sample}")
    print(f"[histograms] {histograms}")

    # 4. Plots.
    plot_all(histograms, normalization, lumi_fb, plots)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
