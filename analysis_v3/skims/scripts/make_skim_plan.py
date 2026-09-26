"""Task plan of the stage-2 skim (program skim_v3), built from the manifests.

    pixi run py -- analysis_v3/skims/scripts/make_skim_plan.py --manifests v2 --version v1

One task per pseudo-data shard; MC files are grouped per sample up to
--mc-bytes of input.  The full analysis configuration is frozen into every
task, so the worker needs no other input.  Outputs go to
production_v3/skims/<version>/<sample>/<task_id>.{root,json}.
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
ANALYSIS_CONFIG = REPO / "analysis_v3/common/config/analysis_ul16_v3.json"


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifests", required=True, help="manifest version, e.g. v2")
    parser.add_argument("--version", required=True, help="skim version, e.g. v1")
    parser.add_argument("--mc-bytes", type=float, default=2.0e9)
    args = parser.parse_args()
    if not re.fullmatch(r"v\d+", args.version) or not re.fullmatch(r"v\d+", args.manifests):
        fail("versions must look like v1, v2, ...")
    manifests = PRODUCTION / "manifests" / args.manifests
    out_dir = PRODUCTION / "skims" / args.version
    plan_path = out_dir / "plan.json"
    if plan_path.exists():
        fail(f"{plan_path} exists; choose a new --version")
    analysis = json.loads(ANALYSIS_CONFIG.read_text(encoding="utf-8"))
    tasks = []

    data = json.loads((manifests / "data.json").read_text(encoding="utf-8"))
    for item in data["files"]:
        task_id = f"skim_{data['name']}_{Path(item['relative']).stem}"
        tasks.append({"task_id": task_id, "config": {"kind": "data", "sample": data["name"], "role": "data"},
                      "analysis_config": analysis,
                      "inputs": [{key: item[key] for key in ("file_key", "path", "relative", "events_entries")}],
                      "outputs": {"json": str(out_dir / data["name"] / f"{task_id}.json"),
                                  "root": str(out_dir / data["name"] / f"{task_id}.root")}})

    for manifest_path in sorted(manifests.glob("mc_*.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        name = manifest["name"]
        groups, current, size = [], [], 0
        for item in manifest["files"]:
            if current and size + item["size_bytes"] > args.mc_bytes:
                groups.append(current)
                current, size = [], 0
            current.append(item)
            size += item["size_bytes"]
        if current:
            groups.append(current)
        for number, group in enumerate(groups):
            task_id = f"skim_{name}_{number:04d}"
            tasks.append({"task_id": task_id, "config": {"kind": "mc", "sample": name, "role": manifest["role"],
                                                         "mode": manifest.get("mode")},
                          "analysis_config": analysis,
                          "inputs": [{key: item[key] for key in ("file_key", "path", "relative", "events_entries")} for item in group],
                          "outputs": {"json": str(out_dir / name / f"{task_id}.json"),
                                      "root": str(out_dir / name / f"{task_id}.root")}})

    payload = {"schema": "h4l_v3_task_plan/1", "stage": f"skims_{args.version}", "program": "skim_v3",
               "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "manifests": str(manifests),
               "manifest_summary_sha256": hashlib.sha256((manifests / "summary.json").read_bytes()).hexdigest(),
               "analysis_config": {"path": str(ANALYSIS_CONFIG), "sha256": hashlib.sha256(ANALYSIS_CONFIG.read_bytes()).hexdigest(),
                                   "version": analysis["version"]},
               "tasks": tasks}
    out_dir.mkdir(parents=True, exist_ok=True)
    temporary = plan_path.with_name(plan_path.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    temporary.rename(plan_path)
    counts = {}
    for task in tasks:
        counts[task["config"]["sample"]] = counts.get(task["config"]["sample"], 0) + 1
    print(f"[plan] {plan_path}: {len(tasks)} tasks {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
