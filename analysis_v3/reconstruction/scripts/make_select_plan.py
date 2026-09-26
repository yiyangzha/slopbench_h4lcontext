"""Task plan of the final selection (program h4l_select, stage 4a).

    pixi run py -- analysis_v3/reconstruction/scripts/make_select_plan.py --reco v2 --version v5 \
        --calibration v4/nominal --lambda-run v4/nominal/lambda_v2_r2 --refit v3 \
        [--selection-config analysis_v3/reconstruction/config/nm1/selection_ul16_v2_nm1_iso.json \
         --samples GluGluToHToZZ_M125 ... --version v5_nm1_iso]

Inputs: the full scan of the 4-lepton event records
production_v3/h4l_reco/<reco>/scan.json (never a directory listing), the
selection configuration analysis_v3/reconstruction/config/selection_ul16_v2.json,
(or --selection-config: the N-1 variants of analysis_v3/reconstruction/config/nm1/),
the lepton calibration payload (applied by the program to the stored raw pT of
every lepton; the event records may have been made with another payload), the lambda
payload, the refit inputs (Z1 line shape, FSR photon resolution) and the
MELA block: the sha256 of every MELA library, dictionary and data file that
MELA reads (the list of the production review's strace and XRootD traces),
verified by every task before MELA is initialized; the data files MELA opens
for update, which every task copies into its empty working directory; the
reference values of the fixed MELA self-test, computed here by the program
(--mela-selftest, in a fresh directory under production_v3/tmp/ with the same
working copies) and reproduced by every task.  One task per event-record
output (of the --samples only, when given: every output of each such sample).
Outputs go to
production_v3/h4l_select/<version>/<sample>/<task_id>.{root,json}; a version is
v<n>, optionally with a suffix (v5_nm1_iso).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
SELECTION_CONFIG = REPO / "analysis_v3/reconstruction/config/selection_ul16_v2.json"
EXTERNAL = PRODUCTION / "external"
PROGRAM = REPO / "analysis_v3/reconstruction/bin/h4l_select"
MELA_DATA = EXTERNAL / "JHUGenMELA/MELA/data"
# The data files MELA opens for update through the links it makes in the working directory (Mela::build):
# every task copies them into its empty working directory first.
MELA_WORKING = ["input.DAT", "process.DAT", "br.sm1", "br.sm2", "ffwarn.dat", "Pdfdata/cteq6l1.tbl", "Pdfdata/cteq6l.tbl",
                "Pdfdata/NNPDF30_lo_as_0130.LHgrid"]
# Everything MELA loads or reads (the shared objects the loader resolves, the dictionary, the data files).
MELA_FILES = ([EXTERNAL / "mela_wrapper/lib/libMelaWrapper.so"] +
              [MELA_DATA / "el9_amd64_gcc11" / name for name in ("libJHUGenMELAMELA.so", "libjhugenmela.so", "libmcfm_711.so",
                                                                 "libcollier.so", "libMG_SMEFTsim_v2.so", "LinkDef_out_rdict.pcm")] +
              [MELA_DATA / name for name in MELA_WORKING] +
              [MELA_DATA / name for name in ("HiggsTotalWidth_YR3.txt", "CombinationInputs/SM_inputs_8TeV/inputs_4mu.txt")] +
              [MELA_DATA / f"pAvgSmooth_{name}.root" for name in (
                  "JHUGen_ZZGG_HSMHiggs", "MCFM_ZZQQB_bkgZZ", "MCFM_ZZGG_HSMHiggs", "MCFM_ZZGG_bkgZZ",
                  "JHUGen_JJVBF_HSMHiggs_13TeV", "JHUGen_JJQCD_HSMHiggs_13TeV", "JHUGen_JQCD_HSMHiggs_13TeV",
                  "JHUGen_Had_WH_HSMHiggs_13TeV", "JHUGen_Had_ZH_HSMHiggs_13TeV", "MCFM_JJVBF_S_HSMHiggs_13TeV",
                  "MCFM_JJVBF_bkgZZ_13TeV", "MCFM_JJQCD_bkgZZ_13TeV", "MCFM_JJQCD_bkgZJets_13TeV_2l2q",
                  "MCFM_Had_WH_S_HSMHiggs_13TeV", "MCFM_Had_WH_bkgZZ_13TeV", "MCFM_Had_ZH_S_HSMHiggs_13TeV",
                  "MCFM_Had_ZH_bkgZZ_13TeV")] +
              [MELA_DATA / f"resolution_mJJ_recoVStrue_{v}_13TeV.root" for v in ("WH", "ZH")])


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reco", required=True, help="event-record version, e.g. v2")
    parser.add_argument("--version", required=True, help="selection version, e.g. v1")
    parser.add_argument("--calibration", required=True, help="<version>/<run> of the calibration payload")
    parser.add_argument("--lambda-run", required=True, help="<version>/<run>/<label> of the lambda payload")
    parser.add_argument("--refit", required=True, help="refit-input version, e.g. v3")
    parser.add_argument("--selection-config", type=Path, default=SELECTION_CONFIG)
    parser.add_argument("--samples", nargs="+", default=None, help="plan only these samples (all their outputs)")
    args = parser.parse_args()
    if not re.fullmatch(r"v\d+(_[a-z0-9]+)*", args.version) or not re.fullmatch(r"v\d+", args.reco):
        fail("versions must look like v1, v2, ... (the selection version may carry a suffix, v5_nm1_iso)")
    selection_config = args.selection_config if args.selection_config.is_absolute() else REPO / args.selection_config
    scan_path = PRODUCTION / "h4l_reco" / args.reco / "scan.json"
    scan = json.loads(scan_path.read_text(encoding="utf-8"))
    if scan.get("schema") != "h4l_v3_reco_scan/1":
        fail(f"unexpected scan schema {scan.get('schema')}")
    out_dir = PRODUCTION / "h4l_select" / args.version
    plan_path = out_dir / "plan.json"
    if plan_path.exists():
        fail(f"{plan_path} exists; choose a new --version")
    selection = json.loads(selection_config.read_text(encoding="utf-8"))
    payload_path = PRODUCTION / "calibration" / args.calibration / "payload.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    lambda_path = PRODUCTION / "calibration" / args.lambda_run / "lambda.json"
    lam = json.loads(lambda_path.read_text(encoding="utf-8"))
    if lam.get("run") != str(PRODUCTION / "calibration" / args.calibration):
        fail(f"the lambda payload belongs to {lam.get('run')}, not to the calibration {args.calibration}")
    refit_path = PRODUCTION / "h4l_reco" / "refit" / args.refit / "refit_inputs.json"
    refit = json.loads(refit_path.read_text(encoding="utf-8"))
    reco_plan = json.loads(Path(scan["plan"]).read_text(encoding="utf-8"))
    if hashlib.sha256(Path(scan["plan"]).read_bytes()).hexdigest() != scan["plan_sha256"]:
        fail("the event-record plan changed after its scan")
    # h4l_select recalibrates every lepton from its stored raw pT with this payload; with the payload of the event
    # records it must reproduce their stored pT (checked by the program).
    same_payload = reco_plan["calibration_payload"]["sha256"] == sha256(payload_path)
    reco_tasks = {t["task_id"]: t for t in reco_plan["tasks"]}
    consumed = [o["task_id"] for o in scan["outputs"]]
    if sorted(consumed) != sorted(reco_tasks) or len(set(consumed)) != len(consumed):
        fail("the scan outputs are not the reco plan tasks, each once")
    known = {o["sample"] for o in scan["outputs"]}
    if args.samples and not set(args.samples) <= known:
        fail(f"unknown samples {sorted(set(args.samples) - known)}")
    mela = {"sqrts_tev": 13.0, "mh": 125.0, "files": [{"path": str(p), "sha256": sha256(p)} for p in MELA_FILES],
            "working_copies": [{"source": str(MELA_DATA / name), "name": name} for name in MELA_WORKING]}
    # Reference values of the MELA self-test, in a fresh directory holding the working copies.
    selftest_dir = PRODUCTION / "tmp" / f"mela_selftest_select_{args.version}_{os.getpid()}"
    selftest_dir.mkdir(parents=True, exist_ok=False)
    for name in MELA_WORKING:
        (selftest_dir / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(MELA_DATA / name, selftest_dir / name)
    selftest_path = selftest_dir / "selftest.json"
    env = {**os.environ, "TMPDIR": str(PRODUCTION / "tmp")}
    result = subprocess.run([str(PROGRAM), "--mela-selftest", str(selftest_path), "--sqrts", str(mela["sqrts_tev"]),
                             "--mh", str(mela["mh"])], cwd=selftest_dir, env=env, capture_output=True, text=True)
    if result.returncode != 0 or not selftest_path.exists():
        fail(f"the MELA self-test failed: {result.stderr.strip()[-2000:]}")
    selftest = json.loads(selftest_path.read_text(encoding="utf-8"))
    if not all(case["ok"] for case in selftest["cases"].values()):
        fail("a MELA self-test case failed")
    mela["selftest"] = {"program": str(PROGRAM), "program_sha256": sha256(PROGRAM), "cases": selftest["cases"]}
    tasks = []
    for output in scan["outputs"]:
        if ".partial." in output["path"] or ".orphan." in output["path"]:
            fail(f"the scan lists a temporary file {output['path']}")
        if args.samples and output["sample"] not in args.samples:
            continue
        reco_task = reco_tasks[output["task_id"]]
        if output["sample"] != reco_task["config"]["sample"]:
            fail(f"{output['task_id']}: the scan sample differs from the reco task")
        task_id = output["task_id"].replace("reco_", "sel_", 1)
        sample = output["sample"]
        tasks.append({
            "task_id": task_id,
            "config": reco_task["config"],
            "selection_config": selection,
            "calibration_payload": payload["flavours"],
            "calibration": {"payload": str(payload_path), "payload_sha256": sha256(payload_path),
                            "event_records_payload_sha256": reco_plan["calibration_payload"]["sha256"],
                            "same_as_event_records": same_payload},
            "lambda_payload": lam,
            "refit_inputs": refit,
            "mela": mela,
            "inputs": [{"reco_task_id": output["task_id"], "root": output["path"], "root_sha256": output["root_sha256"],
                        "rows": output["events4l"]}],
            "source_file_keys": reco_task["source_file_keys"],
            "outputs": {"json": str(out_dir / sample / f"{task_id}.json"), "root": str(out_dir / sample / f"{task_id}.root")}})
    plan = {"schema": "h4l_v3_task_plan/1", "stage": f"h4l_select_{args.version}", "program": "h4l_select",
            "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "reco_scan": str(scan_path), "reco_scan_sha256": sha256(scan_path),
            "selection_config": {"path": str(selection_config), "sha256": sha256(selection_config),
                                 "version": selection["version"]},
            "samples_filter": sorted(args.samples) if args.samples else None,
            "calibration_payload": {"path": str(payload_path), "sha256": sha256(payload_path)},
            "lambda_payload": {"path": str(lambda_path), "sha256": sha256(lambda_path)},
            "refit_inputs": {"path": str(refit_path), "sha256": sha256(refit_path)},
            "mela": mela, "tasks": tasks}
    out_dir.mkdir(parents=True, exist_ok=True)
    temporary = plan_path.with_name(plan_path.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(plan, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    try:
        os.link(temporary, plan_path)
    except FileExistsError:
        fail(f"{plan_path} appeared meanwhile; {temporary} kept")
    counts: dict[str, int] = {}
    for task in tasks:
        counts[task["config"]["sample"]] = counts.get(task["config"]["sample"], 0) + 1
    print(f"[plan] {plan_path}: {len(tasks)} tasks {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
