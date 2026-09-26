"""Task plan of the tag-and-probe pair extraction (program tnp_extract).

    pixi run py -- analysis_v3/tnp/scripts/make_tnp_plan.py --skims v2 --version v2

Inputs are the published stage-2 skim outputs of production_v3/skims/<skims>/,
whose full scan (validation/scan.json) must have passed.  The skims of one
sample are grouped into tasks of about --pairs-per-task Pairs rows.  Every
task carries the "extract" section of the tag-and-probe configuration, the
tag-path masks per flavour resolved against the skim trigger bits, the
planned Pairs count and sha256 of every input and the original NanoAOD file
keys its inputs cover (for the normalization).  Outputs go to
production_v3/tnp/<version>/extract/<sample>/<task_id>.{root,json}.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import sys
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
TNP_CONFIG = REPO / "analysis_v3/tnp/config/tnp_ul16_v1.json"


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skims", required=True, help="skim version, e.g. v2")
    parser.add_argument("--version", required=True, help="tag-and-probe version, e.g. v1")
    parser.add_argument("--pairs-per-task", type=int, default=2_500_000)
    args = parser.parse_args()
    if not re.fullmatch(r"v\d+", args.version) or not re.fullmatch(r"v\d+", args.skims):
        fail("versions must look like v1, v2, ...")
    skims = PRODUCTION / "skims" / args.skims
    skim_plan = json.loads((skims / "plan.json").read_text(encoding="utf-8"))
    scan = json.loads((skims / "validation" / "scan.json").read_text(encoding="utf-8"))
    if scan["valid"] != scan["tasks"] or scan["tasks"] != len(skim_plan["tasks"]):
        fail("the skim scan did not validate every planned task")
    out_dir = PRODUCTION / "tnp" / args.version / "extract"
    plan_path = out_dir / "plan.json"
    if plan_path.exists():
        fail(f"{plan_path} exists; choose a new --version")
    tnp = json.loads(TNP_CONFIG.read_text(encoding="utf-8"))

    analysis = skim_plan["tasks"][0]["analysis_config"]
    if any(task["analysis_config"] != analysis for task in skim_plan["tasks"]):
        fail("the skim tasks do not share one analysis configuration")
    paths = analysis["triggers"]["paths"]
    if len(set(paths)) != len(paths):
        fail("the skim trigger paths are not unique")
    masks = {}
    for flavour in ("muon", "electron"):
        names = tnp["extract"][flavour]["paths"]
        if len(set(names)) != len(names) or not set(names) <= set(paths):
            fail(f"the {flavour} tag paths are not unique paths of the skim trigger mask")
        masks[flavour] = 0
        for name in names:
            masks[flavour] |= 1 << paths.index(name)

    by_sample: dict[str, list[dict]] = {}
    for task in skim_plan["tasks"]:
        report = json.loads(Path(task["outputs"]["json"]).read_text(encoding="utf-8"))
        if report["publication"]["root_path"] != task["outputs"]["root"]:
            fail(f"{task['task_id']}: the publication record names another ROOT file")
        if report["task_id"] != task["task_id"] or report["trigger_bits"] != paths or \
                report["frozen_program"]["sha256"] != scan["program_sha256"]:
            fail(f"{task['task_id']}: the skim report disagrees with the plan or the scan")
        by_sample.setdefault(task["config"]["sample"], []).append({
            "skim_task_id": task["task_id"], "root": task["outputs"]["root"], "json": task["outputs"]["json"],
            "root_sha256": report["publication"]["root_sha256"], "pairs": report["root_output"]["trees"]["Pairs"],
            "file_keys": [item["file_key"] for item in report["files"]], "config": task["config"]})

    tasks = []
    for sample, inputs in sorted(by_sample.items()):
        groups, current, rows = [], [], 0
        for item in inputs:
            if current and rows + item["pairs"] > args.pairs_per_task:
                groups.append(current)
                current, rows = [], 0
            current.append(item)
            rows += item["pairs"]
        if current:
            groups.append(current)
        for number, group in enumerate(groups):
            task_id = f"tnpx_{sample}_{number:04d}"
            tasks.append({
                "task_id": task_id,
                "config": group[0]["config"],
                "extract_config": tnp["extract"],
                "tag_path_masks": masks,
                "inputs": [{key: item[key] for key in ("skim_task_id", "root", "json", "root_sha256", "pairs")} for item in group],
                "source_file_keys": [key for item in group for key in item["file_keys"]],
                "outputs": {"json": str(out_dir / sample / f"{task_id}.json"),
                            "root": str(out_dir / sample / f"{task_id}.root")}})

    payload = {"schema": "h4l_v3_task_plan/1", "stage": f"tnp_extract_{args.version}", "program": "tnp_extract",
               "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "skims": str(skims), "skim_plan_sha256": hashlib.sha256((skims / "plan.json").read_bytes()).hexdigest(),
               "tnp_config": {"path": str(TNP_CONFIG),
                                      "sha256": hashlib.sha256(TNP_CONFIG.read_bytes()).hexdigest(),
                                      "version": tnp["version"]},
               "extract_config": tnp["extract"],
               "tag_path_masks": masks, "tasks": tasks}
    out_dir.mkdir(parents=True, exist_ok=True)
    temporary = plan_path.with_name(plan_path.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    temporary.rename(plan_path)
    counts = {}
    for task in tasks:
        counts[task["config"]["sample"]] = counts.get(task["config"]["sample"], 0) + 1
    print(f"[plan] {plan_path}: {len(tasks)} tasks {counts}, tag path masks {masks}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
