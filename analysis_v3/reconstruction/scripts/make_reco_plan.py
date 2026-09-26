"""Task plan of the 4-lepton event records (program h4l_reco, stage 4a).

    pixi run py -- analysis_v3/reconstruction/scripts/make_reco_plan.py --skims v2 \
        --calibration v3/nominal_r3 --version v1

Inputs are the published stage-2 skim outputs of production_v3/skims/<skims>/,
whose full scan (validation/scan.json) must have passed, and the final lepton
calibration payload production_v3/calibration/<calibration>/payload.json.
The skims of one sample are grouped into tasks of at most --events-per-task
slim events.  Every task embeds the reconstruction configuration
(analysis_v3/reconstruction/config/reco_ul16_v2.json), the skim trigger
configuration (bit order, analysis OR, trigger-object matching) and the
calibration payload, and lists the planned slim-event (and, for signal,
generator-table) counts and sha256 of every input and the original NanoAOD
file keys its inputs cover (for the normalization).  Outputs go to
production_v3/h4l_reco/<version>/<sample>/<task_id>.{root,json}.
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
RECO_CONFIG = REPO / "analysis_v3/reconstruction/config/reco_ul16_v2.json"
REQUIREMENTS = REPO / "analysis_v3/framework/config/input_requirements_v3.json"
THEORY_BRANCHES = ("LHEScaleWeight", "LHEPdfWeight", "PSWeight")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skims", required=True, help="skim version, e.g. v2")
    parser.add_argument("--calibration", required=True, help="<version>/<run> of the calibration payload")
    parser.add_argument("--version", required=True, help="reconstruction version, e.g. v1")
    parser.add_argument("--events-per-task", type=int, default=400_000)
    args = parser.parse_args()
    if not re.fullmatch(r"v\d+", args.version) or not re.fullmatch(r"v\d+", args.skims):
        fail("versions must look like v1, v2, ...")
    skims = PRODUCTION / "skims" / args.skims
    skim_plan = json.loads((skims / "plan.json").read_text(encoding="utf-8"))
    scan = json.loads((skims / "validation" / "scan.json").read_text(encoding="utf-8"))
    if scan["valid"] != scan["tasks"] or scan["tasks"] != len(skim_plan["tasks"]):
        fail("the skim scan did not validate every planned task")
    out_dir = PRODUCTION / "h4l_reco" / args.version
    plan_path = out_dir / "plan.json"
    if plan_path.exists():
        fail(f"{plan_path} exists; choose a new --version")
    reco = json.loads(RECO_CONFIG.read_text(encoding="utf-8"))
    payload_path = PRODUCTION / "calibration" / args.calibration / "payload.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    if set(payload["flavours"]) != {"muon", "electron"}:
        fail(f"{payload_path} does not hold the muon and electron payloads")
    if not payload.get("finalized_by_average") and not payload.get("converged"):
        fail(f"{payload_path} is neither converged nor finalized")
    requirements = json.loads(REQUIREMENTS.read_text(encoding="utf-8"))["events"]["required_role"]

    analysis = skim_plan["tasks"][0]["analysis_config"]
    if any(task["analysis_config"] != analysis for task in skim_plan["tasks"]):
        fail("the skim tasks do not share one analysis configuration")
    triggers = analysis["triggers"]
    if len(set(triggers["paths"])) != len(triggers["paths"]) or not set(triggers["analysis_or"]) <= set(triggers["paths"]):
        fail("the skim trigger paths are not unique or do not contain the analysis OR")

    by_sample: dict[str, list[dict]] = {}
    for task in skim_plan["tasks"]:
        report = json.loads(Path(task["outputs"]["json"]).read_text(encoding="utf-8"))
        if report["publication"]["root_path"] != task["outputs"]["root"]:
            fail(f"{task['task_id']}: the publication record names another ROOT file")
        if report["task_id"] != task["task_id"] or report["trigger_bits"] != triggers["paths"] or \
                report["frozen_program"]["sha256"] != scan["program_sha256"] or \
                report["analysis_config_fnv1a64"] != scan["config"]:
            fail(f"{task['task_id']}: the skim report disagrees with the plan or the scan")
        trees = report["root_output"]["trees"]
        item = {"skim_task_id": task["task_id"], "root": task["outputs"]["root"], "json": task["outputs"]["json"],
                "root_sha256": report["publication"]["root_sha256"], "events": trees["Events"],
                "file_keys": [f["file_key"] for f in report["files"]], "config": task["config"]}
        if task["config"]["role"] == "signal":
            item["gen_rows"] = trees["GenTable"]
        by_sample.setdefault(task["config"]["sample"], []).append(item)

    tasks = []
    for sample, inputs in sorted(by_sample.items()):
        groups, current, events = [], [], 0
        for item in inputs:
            if current and events + item["events"] > args.events_per_task:
                groups.append(current)
                current, events = [], 0
            current.append(item)
            events += item["events"]
        if current:
            groups.append(current)
        for number, group in enumerate(groups):
            task_id = f"reco_{sample}_{number:04d}"
            keys = ("skim_task_id", "root", "json", "root_sha256", "events", "gen_rows", "file_keys")
            role = group[0]["config"]["role"]
            tasks.append({
                "task_id": task_id,
                "config": group[0]["config"],
                "reco_config": reco,
                "triggers": triggers,
                "calibration_payload": payload["flavours"],
                "required_theory_branches": [b for b in requirements.get(role, []) if b in THEORY_BRANCHES],
                "inputs": [{key: item[key] for key in keys if key in item} for item in group],
                "source_file_keys": [key for item in group for key in item["file_keys"]],
                "outputs": {"json": str(out_dir / sample / f"{task_id}.json"),
                            "root": str(out_dir / sample / f"{task_id}.root")}})

    plan = {"schema": "h4l_v3_task_plan/1", "stage": f"h4l_reco_{args.version}", "program": "h4l_reco",
            "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "skims": str(skims), "skim_plan_sha256": sha256(skims / "plan.json"),
            "skim_scan_sha256": sha256(skims / "validation" / "scan.json"),
            "reco_config": {"path": str(RECO_CONFIG), "sha256": sha256(RECO_CONFIG), "version": reco["version"]},
            "calibration_payload": {"path": str(payload_path), "sha256": sha256(payload_path)},
            "tasks": tasks}
    out_dir.mkdir(parents=True, exist_ok=True)
    temporary = plan_path.with_name(plan_path.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(plan, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    try:
        os.link(temporary, plan_path)  # fails instead of replacing a plan created meanwhile
    except FileExistsError:
        fail(f"{plan_path} appeared meanwhile; {temporary} kept")
    counts: dict[str, int] = {}
    for task in tasks:
        counts[task["config"]["sample"]] = counts.get(task["config"]["sample"], 0) + 1
    print(f"[plan] {plan_path}: {len(tasks)} tasks {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
