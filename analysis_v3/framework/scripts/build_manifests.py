"""Build the validated v3 input manifests.

Two steps, run from the repository root with pixi:

    pixi run py -- analysis_v3/framework/scripts/build_manifests.py plan --version v2
    (stage and run the plan with stage_condor.py and program validate_inputs)
    pixi run py -- analysis_v3/framework/scripts/build_manifests.py assemble --version v2

`plan` is the only place in the workflow that lists the input directories.
It assigns every file a stable key (FNV-1a of "<sample>/<relative path>" for
MC and of "<data directory>/<shard>" for the pseudo-data, so keys differ
between datasets) and groups the files into validation tasks.  `assemble`
reads every task output (the full output scan of this stage), requires every
data shard to be valid and all reports to come from one frozen program,
records unreadable or invalid MC files with the coverage they cost, and
writes one manifest per MC sample, the data manifest and a summary.  Nothing
is published unless every check passes, and re-running `assemble` accepts
files that are already published with identical content.  Later stages take
their inputs only from these manifests.
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
CATALOG = REPO / "analysis_v3/framework/config/samples_ul16_v3.json"
REQUIREMENTS = REPO / "analysis_v3/framework/config/input_requirements_v3.json"
NORMALIZATION_RULE = ("per-event weight = genWeight * sigma_eff_pb * 1000 * lumi_fb / S, where S is the sum of the per-file "
                      "genEventSumw over exactly the files a consumer processed (re-sum it over the file_keys you used; "
                      "genEventSumw_total holds only when every listed file is processed)")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def fnv1a64(text: str) -> str:
    value = 14695981039346656037
    for byte in text.encode():
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(payload) -> str:
    return json.dumps(payload, indent=1, sort_keys=True) + "\n"


def publish_text(path: Path, text: str) -> None:
    """Write beside the target, re-read, then rename.  An existing target is
    accepted only if its content is identical; it is never overwritten."""
    if path.exists():
        if path.read_text(encoding="utf-8") == text:
            return
        fail(f"refusing to overwrite {path} with different content")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    if temporary.read_text(encoding="utf-8") != text:
        fail(f"re-read mismatch for {temporary}")
    temporary.rename(path)


def list_root_files(directory: Path) -> list[str]:
    found = []
    for root, _, names in os.walk(directory):
        for name in names:
            if name.endswith(".root") and ".partial." not in name:
                found.append(str(Path(root, name).relative_to(directory)))
    return sorted(found)


def plan(args: argparse.Namespace) -> None:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    requirements = json.loads(REQUIREMENTS.read_text(encoding="utf-8"))
    events = requirements["events"]
    out_dir = PRODUCTION / "manifests" / args.version
    plan_path = out_dir / "plan.json"
    if plan_path.exists():
        fail(f"{plan_path} exists; choose a new --version")
    tasks = []
    listing = {}

    data = catalog["data"]
    data_dir = Path(data["directory"])
    pattern = re.escape(data["pattern"]).replace(r"\*", r"(\d+)")
    shards = sorted(name for name in os.listdir(data_dir) if re.fullmatch(pattern, name))
    numbers = [int(re.fullmatch(pattern, name).group(1)) for name in shards]
    if not shards or numbers != list(range(len(shards))):
        fail("pseudo-data shards are not a contiguous part_00000... sequence")
    listing["pseudo_data"] = len(shards)
    data_config = {"kind": "data", "sample": data["name"], "provenance_object": data["provenance_object"],
                   "expected_schema": data["expected_schema"], "expected_lumi_fb": data["lumi_fb"],
                   "required_events": events["required_all"], "recorded_events": events["recorded"],
                   "full_read": True, "open_attempts": 5, "check_attempts": 3}
    for start in range(0, len(shards), args.data_per_task):
        chunk = shards[start:start + args.data_per_task]
        task_id = f"validate_pseudo_data_{start // args.data_per_task:04d}"
        tasks.append({"task_id": task_id, "config": data_config,
                      "inputs": [{"file_key": fnv1a64(f"{data_dir}/{name}"), "path": str(data_dir / name), "relative": name}
                                 for name in chunk],
                      "outputs": {"json": str(out_dir / "validation" / f"{task_id}.json")}})

    mc_root = Path(catalog["mc"]["root"])
    for sample in catalog["samples"]:
        name = sample["name"]
        files = list_root_files(mc_root / name)
        if not files:
            fail(f"no ROOT files for {name}")
        listing[name] = len(files)
        required = events["required_all"] + events["required_mc"] + events["required_role"].get(sample["role"], [])
        config = {"kind": "mc", "sample": name, "role": sample["role"], "provenance_object": catalog["mc"]["provenance_object"],
                  "expected_schema": catalog["mc"]["expected_schema"], "required_events": required,
                  "recorded_events": events["recorded"], "full_read": False, "open_attempts": 5, "check_attempts": 3}
        for start in range(0, len(files), args.mc_per_task):
            chunk = files[start:start + args.mc_per_task]
            task_id = f"validate_{name}_{start // args.mc_per_task:04d}"
            tasks.append({"task_id": task_id, "config": config,
                          "inputs": [{"file_key": fnv1a64(f"{name}/{relative}"), "path": str(mc_root / name / relative),
                                      "relative": relative} for relative in chunk],
                          "outputs": {"json": str(out_dir / "validation" / f"{task_id}.json")}})

    keys = [item["file_key"] for task in tasks for item in task["inputs"]]
    if len(set(keys)) != len(keys):
        fail("file-key collision")
    payload = {"schema": "h4l_v3_task_plan/1", "stage": f"manifests_{args.version}", "program": "validate_inputs",
               "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "catalog": {"path": str(CATALOG), "sha256": sha256(CATALOG)},
               "requirements": {"path": str(REQUIREMENTS), "sha256": sha256(REQUIREMENTS)},
               "listing": listing, "tasks": tasks}
    publish_text(plan_path, dump(payload))
    print(f"[plan] {plan_path}: {len(tasks)} tasks; files per input: {listing}")


def add_vectors(total: list[float] | None, values: list[float] | None, what: str) -> list[float] | None:
    if values is None:
        return total
    if total is None:
        return list(values)
    if len(total) != len(values):
        fail(f"{what} changes length between files")
    return [a + b for a, b in zip(total, values)]


def assemble(args: argparse.Namespace) -> None:
    out_dir = PRODUCTION / "manifests" / args.version
    plan_path = out_dir / "plan.json"
    plan_data = json.loads(plan_path.read_text(encoding="utf-8"))
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    if sha256(CATALOG) != plan_data["catalog"]["sha256"]:
        fail("the catalogue changed since the plan was made")
    samples = {sample["name"]: sample for sample in catalog["samples"]}

    # Full output scan: every task output exists, parses, comes from the one
    # frozen program, and covers exactly its planned inputs.
    records: dict[str, list] = {}
    programs, configs, missing_outputs = set(), {}, []
    for task in plan_data["tasks"]:
        output = Path(task["outputs"]["json"])
        if not output.exists():
            missing_outputs.append(task["task_id"])
            continue
        report = json.loads(output.read_text(encoding="utf-8"))
        if report.get("schema") != "h4l_v3_validation_report/2":
            fail(f"{output} has schema {report.get('schema')}")
        programs.add((report.get("frozen_program") or {}).get("sha256"))
        configs.setdefault(task["config"]["sample"], set()).add(report["config_fnv1a64"])
        planned = [item["file_key"] for item in task["inputs"]]
        found = [item["file_key"] for item in report["files"]]
        if planned != found:
            fail(f"{output} does not cover exactly the planned inputs")
        records.setdefault(task["config"]["sample"], []).extend(report["files"])
    if missing_outputs:
        fail(f"{len(missing_outputs)} task outputs missing, e.g. {missing_outputs[:5]}")
    if len(programs) != 1 or None in programs:
        fail(f"reports come from {len(programs)} different frozen programs: {sorted(map(str, programs))}")
    for name, values in configs.items():
        if len(values) != 1:
            fail(f"{name}: reports come from different task configurations")

    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    provenance = {"plan": str(plan_path), "plan_sha256": sha256(plan_path), "catalog": plan_data["catalog"],
                  "requirements": plan_data["requirements"], "validator_sha256": programs.pop()}
    outputs: dict[Path, dict] = {}
    summary = {"schema": "h4l_v3_manifest_summary/2", "provenance": provenance, "samples": {}}

    data = catalog["data"]
    shards = records[data["name"]]
    bad = [item for item in shards if item["status"] != "ok" or not item.get("full_read_ok")]
    if bad:
        fail(f"{len(bad)} pseudo-data shards are not valid: " + "; ".join(f"{item['relative']}: {item.get('reason')}" for item in bad[:5]))
    signatures = {item["branch_signature_fnv1a64"] for item in shards}
    if len(signatures) != 1:
        fail("pseudo-data shards differ in their branch list")
    # Per-shard entry counts serve only the completeness checks of later
    # stages; they are never summed or displayed (blinding rules).
    outputs[out_dir / "data.json"] = {
        "schema": "h4l_v3_manifest/1", "kind": "data", "name": data["name"], "directory": data["directory"],
        "lumi_fb": data["lumi_fb"], "provenance": provenance, "branch_signature_fnv1a64": signatures.pop(),
        "files": [{key: item[key] for key in ("file_key", "path", "relative", "events_entries")} for item in shards]}
    summary["samples"][data["name"]] = {"kind": "data", "files": len(shards), "valid": len(shards), "full_read": True}

    for name, sample in samples.items():
        items = records[name]
        valid = [item for item in items if item["status"] == "ok"]
        if not valid:
            fail(f"{name}: no valid file")
        skipped = [{"relative": item["relative"], "status": item["status"], "reason": item.get("reason")} for item in items
                   if item["status"] != "ok"]
        sums = {key: sum(item["mc"][key] for item in valid) for key in
                ("genEventCount", "genEventSumw", "genEventSumw2", "entries_before_selection", "entries_after_selection",
                 "sum_genWeight", "sum_genWeight2", "n_negative_genWeight")}
        lhe = {}
        for key in ("LHEScaleSumw_absolute", "LHEPdfSumw_absolute"):
            total = None
            for item in valid:
                total = add_vectors(total, item["mc"].get(key), f"{name} {key}")
            if total is not None:
                lhe[key] = total
        weight_fraction = sums["sum_genWeight"] / sums["genEventSumw"]
        count_fraction = sums["entries_after_selection"] / sums["entries_before_selection"]
        record_total = sample.get("record_total_events")
        normalization = {
            "method": "genWeight (user decision 2026-09-23)",
            "rule": NORMALIZATION_RULE,
            "sigma_eff_pb": sample["sigma_eff_pb"],
            "genEventSumw_total": sums["genEventSumw"],
            "count_cross_check": "1 * sigma_eff_pb * 1000 * lumi_fb / (sum of entries_before_selection over the processed files)",
            "preselected_fraction_weight": weight_fraction,
            "preselected_fraction_count": count_fraction,
            "ratio_weight_to_count": weight_fraction / count_fraction,
        }
        coverage = {"files_listed": len(items), "files_valid": len(valid), "fraction_files_valid": len(valid) / len(items),
                    "skipped": skipped, "record_total_events": record_total,
                    "fraction_of_record": sums["genEventCount"] / record_total if record_total else None}
        signature_count = len({item["branch_signature_fnv1a64"] for item in valid})
        vector_lengths = sorted({json.dumps({key: item["mc"].get(key) for key in
                                            ("first_entry_nLHEScaleWeight", "first_entry_nLHEPdfWeight", "first_entry_nPSWeight")},
                                           sort_keys=True) for item in valid})
        outputs[out_dir / f"mc_{name}.json"] = {
            "schema": "h4l_v3_manifest/1", "kind": "mc", "name": name, "role": sample["role"],
            "mode": sample.get("mode"), "mass_GeV": sample.get("mass_GeV"), "note": sample.get("note"),
            "provenance": provenance, "normalization": normalization, "totals": sums, "lhe_sums_before_selection": lhe,
            "coverage": coverage, "branch_signatures": signature_count, "weight_vector_lengths": vector_lengths,
            "files": [{"file_key": item["file_key"], "path": item["path"], "relative": item["relative"],
                       "size_bytes": item["size_bytes"], "events_entries": item["events_entries"],
                       "genEventSumw": item["mc"]["genEventSumw"], "genEventCount": item["mc"]["genEventCount"],
                       "entries_before_selection": item["mc"]["entries_before_selection"],
                       "sum_genWeight": item["mc"]["sum_genWeight"],
                       **{key: item["mc"][key] for key in ("LHEScaleSumw_absolute", "LHEPdfSumw_absolute") if key in item["mc"]},
                       "recorded_present": item["recorded_present"]}
                      for item in valid]}
        summary["samples"][name] = {"kind": "mc", "files": len(items), "valid": len(valid), "skipped": skipped,
                                    "events_selected": sums["entries_after_selection"],
                                    "entries_before_selection": sums["entries_before_selection"],
                                    "fraction_of_record": coverage["fraction_of_record"],
                                    "genEventSumw": sums["genEventSumw"], "sum_genWeight_selected": sums["sum_genWeight"],
                                    "negative_fraction_selected": sums["n_negative_genWeight"] / sums["entries_after_selection"],
                                    "ratio_weight_to_count": normalization["ratio_weight_to_count"],
                                    "size_GB": sum(item["size_bytes"] for item in valid) / 1e9,
                                    "branch_signatures": signature_count, "weight_vector_lengths": vector_lengths,
                                    "lhe_sums": sorted(lhe)}
    outputs[out_dir / "summary.json"] = summary

    lines = [f"# v3 input manifests {args.version}", "", f"Assembled from `{plan_path}`.", "",
             "| sample | valid / listed files | preselected events | fraction of record | weight/count ratio | negative w | size [GB] |",
             "|---|---|---|---|---|---|---|"]
    for name, entry in summary["samples"].items():
        if entry["kind"] == "data":
            lines.append(f"| {name} | {entry['valid']} / {entry['files']} (full read) | - | - | - | - | - |")
        else:
            record = f"{entry['fraction_of_record']:.4f}" if entry["fraction_of_record"] else "-"
            lines.append(f"| {name} | {entry['valid']} / {entry['files']} | {entry['events_selected']} | {record} | "
                         f"{entry['ratio_weight_to_count']:.5f} | {entry['negative_fraction_selected']:.4f} | {entry['size_GB']:.2f} |")
    for path, payload in outputs.items():
        publish_text(path, dump(payload))
    publish_text(out_dir / "summary.md", "\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[assembled] {out_dir} ({stamp})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("plan")
    one.add_argument("--version", required=True)
    one.add_argument("--mc-per-task", type=int, default=100)
    one.add_argument("--data-per-task", type=int, default=5)
    two = sub.add_parser("assemble")
    two.add_argument("--version", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"v\d+", args.version):
        fail("--version must look like v1, v2, ...")
    (plan if args.command == "plan" else assemble)(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
