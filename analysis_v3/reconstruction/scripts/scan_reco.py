"""Full output scan of the 4-lepton event records (program h4l_reco).

    pixi run py -- analysis_v3/reconstruction/scripts/scan_reco.py --version v2 --manifests v2

For every planned task of production_v3/h4l_reco/<version>/plan.json:
* the JSON exists, parses and is bound to its plan task: task_id, sample,
  kind, role, the inputs (skim task, ROOT path, planned sha256, slim-event
  count, generator-table rows), the source file keys and the hashes of the
  plan-embedded reconstruction configuration, trigger configuration and
  calibration payload;
* every output comes from one frozen program;
* the publication record names the planned ROOT, whose sha256 is recomputed
  and must equal the record, and root_check passes with the JSON tree counts;
* the event totals equal the planned slim events, Events4l equals kept and
  Inputs the number of inputs; the candidate, lepton and jet overflow
  counters are reported (a sample with overflows is listed, not hidden);
* the per-sample union of the original NanoAOD file keys equals the manifest
  file keys (complete coverage).
Writes production_v3/h4l_reco/<version>/scan.json (never overwritten) with,
per output, the path, sha256 and Events4l count that the selection consumes;
consumers take the paths from this file and never glob (".partial.*" and
".orphan.*" files also end in ".root").
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
ROOT_CHECK = REPO / "analysis_v3/framework/bin/root_check"


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def fnv1a64(text: str) -> str:
    value = 14695981039346656037
    for byte in text.encode("utf-8"):
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def compact(value) -> str:
    """nlohmann::json::dump() of the task JSON: sorted keys, no spaces."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def check_output(task: dict) -> dict:
    record = {"task_id": task["task_id"]}
    path = Path(task["outputs"]["json"])
    if not path.exists():
        return {**record, "error": "JSON output missing"}
    report = json.loads(path.read_text(encoding="utf-8"))
    problems = []
    config = task["config"]
    for key, expected in (("task_id", task["task_id"]), ("sample", config["sample"]), ("kind", config["kind"]),
                          ("role", config["role"]), ("reco_config_version", task["reco_config"]["version"])):
        if report.get(key) != expected:
            problems.append(f"{key} {report.get(key)!r} != {expected!r}")
    for key, value in (("reco_config_fnv1a64", task["reco_config"]), ("triggers_fnv1a64", task["triggers"]),
                       ("calibration_payload_fnv1a64", task["calibration_payload"])):
        if report.get(key) != fnv1a64(compact(value)):
            problems.append(f"{key} differs from the plan-embedded value")
    fields = ("skim_task_id", "root", "root_sha256", "events")
    planned = [tuple(item[k] for k in fields) + (item.get("gen_rows"),) for item in task["inputs"]]
    produced = [tuple(item[k] for k in fields) + (item.get("gen_rows"),) for item in report["inputs"]]
    if planned != produced:
        problems.append("processed inputs differ from the plan")
    if report.get("source_file_keys") != task["source_file_keys"]:
        problems.append("source file keys differ from the plan")
    totals = report["totals"]
    if totals["events"] != sum(item["events"] for item in task["inputs"]):
        problems.append("processed events differ from the planned slim events")
    trees = report["root_output"]["trees"]
    if trees.get("Events4l") != totals["kept"]:
        problems.append("Events4l differs from kept")
    if trees.get("Inputs") != len(task["inputs"]):
        problems.append("Inputs tree count differs from the plan")
    publication = report.get("publication") or {}
    if publication.get("root_path") != task["outputs"]["root"] or not publication.get("root_sha256"):
        problems.append("no publication record of the planned ROOT")
    elif sha256(task["outputs"]["root"]) != publication["root_sha256"]:
        problems.append("ROOT sha256 differs from the publication record")
    arguments = [str(ROOT_CHECK), task["outputs"]["root"]] + [f"{name}={count}" for name, count in trees.items()]
    result = subprocess.run(arguments, capture_output=True, text=True)
    if result.returncode != 0:
        problems.append(f"root_check failed: {result.stderr.strip()}")
    if problems:
        return {**record, "error": "; ".join(problems)}
    return {**record, "sample": config["sample"], "kind": config["kind"], "role": config["role"],
            "path": task["outputs"]["root"], "root_sha256": publication["root_sha256"], "events4l": trees["Events4l"],
            "with_candidate": totals["with_candidate"], "candidates": totals["candidates"],
            "overflow": {k: totals[k] for k in ("candidate_overflow", "lepton_overflow", "jet_overflow")},
            "program_sha256": (report.get("frozen_program") or {}).get("sha256"),
            "source_file_keys": task["source_file_keys"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--manifests", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    base = PRODUCTION / "h4l_reco" / args.version
    plan = json.loads((base / "plan.json").read_text(encoding="utf-8"))
    target = base / "scan.json"
    if target.exists():
        fail(f"{target} exists; it is never overwritten")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(check_output, plan["tasks"]))
    bad = [r for r in records if "error" in r]
    if bad:
        fail(f"{len(bad)} of {len(records)} outputs invalid, e.g. " + "; ".join(f"{r['task_id']}: {r['error']}" for r in bad[:5]))
    programs = {r["program_sha256"] for r in records}
    if len(programs) != 1 or None in programs:
        fail(f"outputs come from {len(programs)} programs")

    manifests = PRODUCTION / "manifests" / args.manifests
    by_name = {}
    for manifest_path in [manifests / "data.json", *sorted(manifests.glob("mc_*.json"))]:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        by_name[manifest["name"]] = manifest
    samples = {}
    for record in records:
        sample = samples.setdefault(record["sample"], {"kind": record["kind"], "role": record["role"], "outputs": [],
                                                       "file_keys": [], "overflow": {}})
        sample["outputs"].append(record["task_id"])
        sample["file_keys"].extend(record["source_file_keys"])
        for key, value in record["overflow"].items():
            sample["overflow"][key] = sample["overflow"].get(key, 0) + value
    for name, sample in samples.items():
        manifest = by_name.get(name)
        if manifest is None:
            fail(f"no manifest for sample {name}")
        keys = sample["file_keys"]
        if len(set(keys)) != len(keys):
            fail(f"{name}: an original file is covered twice")
        if set(keys) != {item["file_key"] for item in manifest["files"]}:
            fail(f"{name}: the reconstruction does not cover exactly the manifest files")
        if sample["kind"] == "mc":
            sample["genEventSumw"] = sum(item["genEventSumw"] for item in manifest["files"])
            sample["sigma_eff_pb"] = manifest["normalization"]["sigma_eff_pb"]
        else:
            sample["lumi_fb"] = manifest["lumi_fb"]
        for key in ("events4l", "with_candidate", "candidates"):
            sample[key] = sum(r[key] for r in records if r["sample"] == name)
        del sample["file_keys"]
        sample["files"] = len(keys)
    payload = {"schema": "h4l_v3_reco_scan/1", "version": args.version, "manifests": str(manifests),
               "plan": str(base / "plan.json"), "plan_sha256": hashlib.sha256((base / "plan.json").read_bytes()).hexdigest(),
               "program_sha256": programs.pop(),
               "outputs": [{key: r[key] for key in ("task_id", "sample", "kind", "role", "path", "root_sha256", "events4l",
                                                    "with_candidate", "candidates", "overflow")} for r in records],
               "samples": samples}
    temporary = target.with_name(target.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    try:
        os.link(temporary, target)
    except FileExistsError:
        fail(f"{target} appeared during the scan; kept {temporary}")
    overflows = {name: s["overflow"] for name, s in samples.items() if any(s["overflow"].values())}
    print(f"[scan] {len(records)} outputs valid, coverage complete; events4l "
          f"{ {name: s['events4l'] for name, s in samples.items()} }\n[scan] overflows: {overflows or 'none'}\n[scan] {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
