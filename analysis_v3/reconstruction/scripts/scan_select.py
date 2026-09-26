"""Full output scan of the final selection (program h4l_select, stage 4a).

    pixi run py -- analysis_v3/reconstruction/scripts/scan_select.py --version v1 --manifests v2

For every planned task of production_v3/h4l_select/<version>/plan.json:
* the JSON exists, parses and is bound to its plan task: schema, task_id,
  sample, kind, role, the selection version, the hashes of the plan-embedded
  selection configuration, calibration payload, lambda payload, refit inputs
  and MELA block, the processed inputs (reco task, ROOT path, planned sha256,
  rows), the source file keys, the MELA files (paths and sha256 of the plan)
  and working copies, and a passed MELA self-test;
* every output comes from one frozen program;
* the publication record names the planned ROOT, whose sha256 is recomputed
  and must equal the record, and root_check passes with the JSON tree counts;
* the processed rows equal the planned rows; every signal-region tree equals
  its selected count, CR the sum of the control-region counts, ZL the Z + l
  count and Inputs the number of inputs; the refit counters add up to the
  filled rows and SR does not exceed SRZ4l (SR events are SRZ4l events);
* the per-sample union of the original NanoAOD file keys equals the manifest
  file keys (complete coverage).
Writes production_v3/h4l_select/<version>/scan.json (never overwritten) with,
per output, the path, sha256 and tree counts, and per sample the totals, the
genEventSumw of the manifest files and sigma_eff (MC) or the luminosity
(data), plus the outputs whose MELA initialization or input opening logged
ROOT errors (recorded, not fatal: the self-test decides); consumers take the
paths from this file and never glob.
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


def canonical(value) -> str:
    """The canonical text of h4l_select's configuration hashes: sorted keys, no spaces, floats as %.17g."""
    if isinstance(value, dict):
        return "{" + ",".join(json.dumps(k, ensure_ascii=False) + ":" + canonical(value[k]) for k in sorted(value)) + "}"
    if isinstance(value, list):
        return "[" + ",".join(canonical(item) for item in value) + "]"
    if isinstance(value, bool) or value is None:
        return json.dumps(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return "%.17g" % value
    return json.dumps(value, ensure_ascii=False)


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def check_output(task: dict, regions: list[str], mela_files: list[dict]) -> dict:
    record = {"task_id": task["task_id"]}
    path = Path(task["outputs"]["json"])
    if not path.exists():
        return {**record, "error": "JSON output missing"}
    report = json.loads(path.read_text(encoding="utf-8"))
    problems = []
    config = task["config"]
    for key, expected in (("schema", "h4l_v3_select_report/2"), ("task_id", task["task_id"]), ("sample", config["sample"]),
                          ("kind", config["kind"]), ("role", config["role"]),
                          ("selection_version", task["selection_config"]["version"]), ("mela_selftest", "passed")):
        if report.get(key) != expected:
            problems.append(f"{key} {report.get(key)!r} != {expected!r}")
    for block in ("selection_config", "calibration_payload", "lambda_payload", "refit_inputs", "mela"):
        if report.get(f"{block}_fnv1a64") != fnv1a64(canonical(task[block])):
            problems.append(f"{block}_fnv1a64 differs from the plan-embedded value")
    if report.get("mela_working_copies") != [{"name": w["name"], "source": w["source"]} for w in task["mela"]["working_copies"]]:
        problems.append("MELA working copies differ from the plan")
    if "calibration" in task and report.get("calibration") != task["calibration"]:
        problems.append("calibration block differs from the plan")
    planned = [(i["reco_task_id"], i["root"], i["root_sha256"], i["rows"]) for i in task["inputs"]]
    produced = [(i["reco_task_id"], i["root"], i["root_sha256"], i["rows"]) for i in report.get("inputs", [])]
    if planned != produced:
        problems.append("processed inputs differ from the plan")
    if report.get("source_file_keys") != task["source_file_keys"]:
        problems.append("source file keys differ from the plan")
    if report.get("mela_files") != mela_files:
        problems.append("MELA libraries differ from the plan")
    counts = report.get("counts", {})
    trees = report.get("root_output", {}).get("trees", {})
    if counts.get("rows") != sum(i["rows"] for i in task["inputs"]):
        problems.append("processed rows differ from the plan")
    for region in regions:
        if trees.get(region) != counts.get(f"selected_{region}"):
            problems.append(f"tree {region} differs from its selected count")
    if trees.get("CR") != sum(v for k, v in counts.items() if k.startswith("cr_")):
        problems.append("tree CR differs from the control-region counts")
    if trees.get("ZL") != counts.get("zl"):
        problems.append("tree ZL differs from the Z + l count")
    if trees.get("Inputs") != len(task["inputs"]):
        problems.append("tree Inputs differs from the number of inputs")
    filled = sum(trees.get(region, 0) for region in regions) + trees.get("CR", 0)
    if counts.get("refit_ok", 0) + counts.get("refit_failed", 0) != filled:
        problems.append("refit counters differ from the filled rows")
    if "SRZ4l" in trees and trees.get("SR", 0) > trees["SRZ4l"]:
        problems.append("more SR than SRZ4l events")
    publication = report.get("publication") or {}
    if publication.get("root_path") != task["outputs"]["root"] or not publication.get("root_sha256"):
        problems.append("no publication record of the planned ROOT")
    elif sha256(task["outputs"]["root"]) != publication["root_sha256"]:
        problems.append("ROOT sha256 differs from the publication record")
    else:
        arguments = [str(ROOT_CHECK), task["outputs"]["root"]] + [f"{name}={count}" for name, count in trees.items()]
        result = subprocess.run(arguments, capture_output=True, text=True)
        if result.returncode != 0:
            problems.append(f"root_check failed: {result.stderr.strip()}")
    if problems:
        return {**record, "error": "; ".join(problems)}
    return {**record, "sample": config["sample"], "kind": config["kind"], "role": config["role"],
            "path": task["outputs"]["root"], "root_sha256": publication["root_sha256"], "trees": trees, "counts": counts,
            "program_sha256": (report.get("frozen_program") or {}).get("sha256"),
            "source_file_keys": task["source_file_keys"],
            "root_errors": {"mela_init": report.get("mela_init_root_error"), "input_open": report.get("input_open_root_errors")}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--manifests", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    base = PRODUCTION / "h4l_select" / args.version
    plan_path = base / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    target = base / "scan.json"
    if target.exists():
        fail(f"{target} exists; it is never overwritten")
    regions = [r["name"] for r in plan["tasks"][0]["selection_config"]["candidate"]["regions"]]
    mela_files = plan["mela"]["files"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(lambda task: check_output(task, regions, mela_files), plan["tasks"]))
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
                                                       "file_keys": [], "trees": {}, "counts": {}})
        sample["outputs"].append(record["task_id"])
        sample["file_keys"].extend(record["source_file_keys"])
        for block in ("trees", "counts"):
            for key, value in record[block].items():
                sample[block][key] = sample[block].get(key, 0) + value
    for name, sample in samples.items():
        manifest = by_name.get(name)
        if manifest is None:
            fail(f"no manifest for sample {name}")
        keys = sample["file_keys"]
        if len(set(keys)) != len(keys):
            fail(f"{name}: an original file is covered twice")
        if set(keys) != {item["file_key"] for item in manifest["files"]}:
            fail(f"{name}: the selection does not cover exactly the manifest files")
        if sample["kind"] == "mc":
            sample["genEventSumw"] = sum(item["genEventSumw"] for item in manifest["files"])
            sample["sigma_eff_pb"] = manifest["normalization"]["sigma_eff_pb"]
        else:
            sample["lumi_fb"] = manifest["lumi_fb"]
        del sample["file_keys"]
        sample["files"] = len(keys)
    logged = {r["task_id"]: r["root_errors"] for r in records if r["root_errors"]["mela_init"] or r["root_errors"]["input_open"]}
    payload = {"schema": "h4l_v3_select_scan/1", "version": args.version, "manifests": str(manifests),
               "plan": str(plan_path), "plan_sha256": sha256(plan_path), "program_sha256": programs.pop(),
               "selection_version": plan["selection_config"]["version"], "root_errors_logged": logged,
               "outputs": [{key: r[key] for key in ("task_id", "sample", "kind", "role", "path", "root_sha256", "trees", "counts")}
                           for r in records],
               "samples": samples}
    temporary = target.with_name(target.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    try:
        os.link(temporary, target)
    except FileExistsError:
        fail(f"{target} appeared during the scan; kept {temporary}")
    summary = {name: {r: s["trees"].get(r, 0) for r in regions} for name, s in samples.items()}
    print(f"[scan] {len(records)} outputs valid, coverage complete; selected {summary}\n"
          f"[scan] outputs with logged ROOT errors: {len(logged)}\n[scan] {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
