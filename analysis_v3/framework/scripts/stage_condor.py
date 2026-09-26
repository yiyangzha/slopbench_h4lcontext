"""Stage one v3 task plan on HTCondor through eossubmit.

A plan (JSON) lists independent tasks; every task names its program, its
inputs and its final outputs below production_v3/.  This script

* freezes the program, root_check and the worker wrapper by content under
  production_v3/program_inputs/<stage>/<sha12>/;
* writes the run directory production_v3/condor/<stage>/<UTC stamp>/ with one
  task file per pending task, the queue, the JDL, logs/ and records/;
* skips every task whose final JSON output exists, and (for --dry-run and
  --submit) every task that is still in the queue of an earlier submission
  of the same stage;
* runs one task locally through the frozen wrapper with PATH=/usr/bin:/bin
  and records its peak memory (--local-smoke); a passed smoke of the same
  frozen programs is required before --submit;
* checks the JDL with condor_submit -dry-run (--dry-run) or submits it
  (--submit), verifying Iwd, Arguments, Environment and the log paths of
  every job, and records the cluster ID.

Run from the repository root with pixi, after
`source analysis_v3/framework/scripts/eossubmit_env.sh` for --dry-run/--submit:

    pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan PLAN.json --stage NAME \
        --program analysis_v3/framework/bin/validate_inputs --memory-mb 1000 --flavour espresso --submit
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
WRAPPER = REPO / "analysis_v3/framework/condor/run_task.sh"
ROOT_CHECK = REPO / "analysis_v3/framework/bin/root_check"
QUEUE_MACROS = ("task_id", "task_path", "task_file", "program_name")
RESERVED = {"output", "error", "log", "path", "input", "arguments", "executable", "universe", "environment"}
PROXY = "/afs/cern.ch/user/y/yiyangz/private/x509up_u165165.backup_20260923"
# Program time limit per job flavour (seconds), leaving room for the publication steps.
PROGRAM_TIMEOUT = {"espresso": 900, "microcentury": 2400, "longlunch": 5400, "workday": 24000, "tomorrow": 72000}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def utc(fmt: str = "%Y-%m-%dT%H:%M:%SZ") -> str:
    return dt.datetime.now(dt.timezone.utc).strftime(fmt)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_path(path: Path) -> Path:
    """Absolute /eos/user form; never resolve(), which yields /eos/home-y."""
    return path if path.is_absolute() else REPO / path


def inside_production(path: str) -> bool:
    return path.startswith(str(PRODUCTION) + "/")


def freeze(stage: str, program: Path) -> dict:
    """Content-addressed copies of the program, root_check and the wrapper."""
    files = {"program": program, "root_check": ROOT_CHECK, "wrapper": WRAPPER}
    digests = {name: sha256(path) for name, path in files.items()}
    key = hashlib.sha256("".join(digests[name] for name in sorted(digests)).encode()).hexdigest()[:12]
    directory = PRODUCTION / "program_inputs" / stage / key
    directory.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for name, source in files.items():
        target = directory / source.name
        if target.exists():
            if sha256(target) != digests[name]:
                fail(f"frozen copy differs from its content address: {target}")
        else:
            temporary = target.with_name(target.name + f".partial.{os.getpid()}")
            shutil.copy2(source, temporary)
            if sha256(temporary) != digests[name]:
                fail(f"copy of {source} does not reproduce its sha256")
            temporary.rename(target)
        os.chmod(target, 0o755)
        frozen[name] = {"source": str(source), "frozen": str(target), "sha256": digests[name]}
    return {"key": key, "directory": str(directory), "files": frozen}


def write_jdl(path: Path, run_dir: Path, frozen: dict, args: argparse.Namespace, queue: Path) -> None:
    program = frozen["files"]["program"]["frozen"]
    root_check = frozen["files"]["root_check"]["frozen"]
    wrapper = frozen["files"]["wrapper"]["frozen"]
    text = f"""universe = vanilla
initialdir = {run_dir}
executable = {wrapper}
arguments = "$(task_file) $(program_name)"
environment = "H4L_CLUSTER=$(ClusterId) H4L_PROC=$(ProcId)"
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_executable = True
transfer_input_files = {program},{root_check},$(task_path)
transfer_output_files = ""
output = {run_dir}/logs/$(task_id).$(ClusterId).$(ProcId).out
error = {run_dir}/logs/$(task_id).$(ClusterId).$(ProcId).err
log = {run_dir}/logs/cluster.$(ClusterId).log
use_x509userproxy = true
+JobFlavour = "{args.flavour}"
MY.WantOS = "el9"
request_cpus = 1
request_memory = {args.memory_mb}
request_disk = {args.disk_mb * 1024}
max_retries = 2
queue {",".join(QUEUE_MACROS)} from {queue}
"""
    path.write_text(text, encoding="utf-8")


def in_flight_tasks(stage: str) -> set[str]:
    """Task identifiers still queued by earlier submissions of this stage."""
    clusters = []
    for submission in sorted((PRODUCTION / "condor" / stage).glob("*/submission.json")):
        clusters.append(str(json.loads(submission.read_text(encoding="utf-8"))["cluster"]))
    if not clusters:
        return set()
    result = subprocess.run(["condor_q", *clusters, "-af:t", "Arguments", "Args", "JobStatus"], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        fail("condor_q failed; cannot check for in-flight tasks")
    tasks = set()
    for line in result.stdout.splitlines():
        for field in line.split("\t")[:2]:
            first = field.split()[0] if field.split() else ""
            if first.endswith(".json"):
                tasks.add(first[: -len(".json")])
    return tasks


def local_smoke(run_dir: Path, frozen: dict, task_path: Path) -> int:
    """Run one task through the frozen wrapper in a sandbox with a minimal PATH."""
    sandbox = PRODUCTION / "tmp" / "smoke" / run_dir.relative_to(PRODUCTION / "condor") / task_path.stem
    sandbox.mkdir(parents=True, exist_ok=False)
    for name in ("program", "root_check", "wrapper"):
        shutil.copy2(frozen["files"][name]["frozen"], sandbox / Path(frozen["files"][name]["frozen"]).name)
    shutil.copy2(task_path, sandbox / task_path.name)
    program_name = Path(frozen["files"]["program"]["frozen"]).name
    environment = {"PATH": "/usr/bin:/bin", "HOME": str(sandbox), "TMPDIR": str(sandbox), "LANG": "C"}
    command = ["/usr/bin/time", "-v", "/bin/bash", Path(frozen["files"]["wrapper"]["frozen"]).name, task_path.name,
               program_name]
    print(f"[smoke] sandbox {sandbox}\n[smoke] {' '.join(command)} (PATH=/usr/bin:/bin)")
    result = subprocess.run(command, cwd=sandbox, env=environment, stderr=subprocess.PIPE, text=True)
    print(result.stderr, file=sys.stderr)
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", result.stderr)
    peak_mb = int(match.group(1)) / 1024 if match else None
    print(f"[smoke] exit code {result.returncode}, peak memory {peak_mb:.0f} MB" if peak_mb else
          f"[smoke] exit code {result.returncode}")
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", result.stderr)
    task = json.loads(task_path.read_text(encoding="utf-8"))
    if result.returncode == 0 and Path(task["outputs"]["json"]).exists():
        root = task["outputs"].get("root")
        marker = Path(frozen["directory"]) / f"smoke_passed.{task['task_id']}.json"
        marker.write_text(json.dumps({"task_id": task["task_id"], "utc": utc(), "peak_memory_mb": peak_mb,
                                      "wall_clock": wall.group(1) if wall else None,
                                      "output_bytes": Path(root).stat().st_size if root and Path(root).exists() else None,
                                      "sandbox": str(sandbox)}, indent=1) + "\n", encoding="utf-8")
        print(f"[smoke] passed; marker {marker}")
    return result.returncode


def check_dry_run(text: str, run_dir: Path, pending: list[dict], program_name: str) -> None:
    procs = set(re.findall(r"^ProcId=(\d+)", text, flags=re.M))
    if len(procs) != len(pending):
        fail(f"dry run has {len(procs)} distinct ProcIds for {len(pending)} tasks")
    arguments = re.findall(r'^(?:Arguments|Args)="([^"]*)"', text, flags=re.M)
    expected = {f"{task['task_id']}.json {program_name}" for task in pending}
    if set(arguments) != expected:
        missing = sorted(expected - set(arguments))[:3]
        fail(f"dry-run Arguments differ from the queue (e.g. missing {missing})")
    for iwd in re.findall(r'^Iwd="([^"]*)"', text, flags=re.M):
        if iwd != str(run_dir):
            fail(f"dry-run Iwd is {iwd}, expected {run_dir}")
    if not re.search(r"^Iwd=", text, flags=re.M):
        fail("dry-run classads lack Iwd")
    if "H4L_CLUSTER=" not in text:
        fail("dry-run classads lack the H4L_CLUSTER environment")
    for path in re.findall(r'^(?:Out|Err|UserLog)="([^"]+)"', text, flags=re.M):
        if not path.startswith(str(run_dir) + "/logs/"):
            fail(f"log path outside the run directory: {path}")
    for needed in ("Cmd=", "Out=", "Err=", "UserLog=", "TransferInput=", "RequestMemory=", "RequestDisk=", "JobFlavour="):
        if needed not in text:
            fail(f"dry-run classads lack {needed}")
    proxies = set(re.findall(r'^x509userproxy="([^"]*)"', text, flags=re.M))
    if proxies != {PROXY}:
        fail(f"dry-run x509userproxy is {sorted(proxies)}, expected {PROXY}")
    for transfer in re.findall(r'^TransferInput="([^"]*)"', text, flags=re.M):
        if program_name not in transfer or "root_check" not in transfer or ".json" not in transfer:
            fail(f"dry-run TransferInput lacks the program, root_check or the task: {transfer}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--stage", required=True, help="stage name, used for the program_inputs and condor directories")
    parser.add_argument("--program", type=Path, required=True, help="built C++ program implementing --task/--out-json")
    parser.add_argument("--memory-mb", type=int, default=2000, help="request_memory; keep <= 2000 for one slot")
    parser.add_argument("--disk-mb", type=int, default=2000)
    parser.add_argument("--flavour", default="espresso")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--local-smoke", metavar="TASK_ID", help="run this task locally through the wrapper")
    mode.add_argument("--dry-run", action="store_true", help="condor_submit -dry-run of the staged JDL")
    mode.add_argument("--submit", action="store_true", help="submit the staged JDL")
    args = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9_]+", args.stage):
        fail("--stage must be alphanumeric/underscore")
    for macro in QUEUE_MACROS:
        if macro.lower() in RESERVED:
            fail(f"reserved queue macro {macro}")
    plan_path = repo_path(args.plan)
    program = repo_path(args.program)
    for path in (program, ROOT_CHECK, WRAPPER):
        if not path.is_file() or not os.access(path, os.X_OK):
            fail(f"missing or non-executable: {path}")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("stage") != args.stage:
        fail(f"--stage {args.stage} differs from the plan stage {plan.get('stage')}")
    if args.flavour not in PROGRAM_TIMEOUT:
        fail(f"unknown job flavour {args.flavour}")
    tasks = plan["tasks"]
    identifiers = [task["task_id"] for task in tasks]
    if len(set(identifiers)) != len(identifiers):
        fail("duplicate task identifiers in the plan")

    cluster_mode = args.dry_run or args.submit
    if cluster_mode and os.environ.get("_myschedd_POOL") != "eossubmit":
        fail("source analysis_v3/framework/scripts/eossubmit_env.sh first (eossubmit pool not loaded)")
    if cluster_mode:
        if os.environ.get("X509_USER_PROXY") != PROXY:
            fail(f"X509_USER_PROXY must be {PROXY}")
        result = subprocess.run(["voms-proxy-info", "-file", PROXY, "-timeleft"], capture_output=True, text=True)
        left = int(result.stdout.strip() or 0) if result.returncode == 0 and result.stdout.strip().isdigit() else 0
        if left < 86400:
            fail(f"the proxy has {left} s left (< 86400); ask the user to renew it")
        print(f"[proxy] {PROXY}: {left} s left")
    flying = in_flight_tasks(args.stage) if cluster_mode else set()
    program_sha = sha256(program)
    pending = []
    complete = 0
    foreign = []
    for task in tasks:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", task["task_id"]):
            fail(f"task identifier with unsafe characters: {task['task_id']}")
        for key in ("json", "root"):
            value = task["outputs"].get(key)
            if value and not inside_production(value):
                fail(f"output outside production_v3: {value}")
        output = Path(task["outputs"]["json"])
        if output.exists():
            # A completed output counts only if the same program produced it.
            produced_by = (json.loads(output.read_text(encoding="utf-8")).get("frozen_program") or {}).get("sha256")
            if produced_by != program_sha:
                foreign.append(task["task_id"])
            complete += 1
        elif task["task_id"] not in flying:
            pending.append(task)
    if foreign:
        fail(f"{len(foreign)} existing outputs come from another program version (e.g. {foreign[:3]}); "
             "use a new plan version instead of mixing programs")
    print(f"[plan] {len(tasks)} tasks: {complete} complete, {len(flying)} in flight, {len(pending)} to stage")
    if not pending:
        return 0

    # Create every output directory here, once: many workers creating the same
    # EOS directory at the same time fail to see it through their FUSE mounts.
    for directory in sorted({str(Path(task["outputs"][key]).parent) for task in pending for key in ("json", "root")
                             if task["outputs"].get(key)}):
        Path(directory).mkdir(parents=True, exist_ok=True)
        if not Path(directory).is_dir():
            fail(f"cannot create the output directory {directory}")
    frozen = freeze(args.stage, program)
    if args.submit and not list(Path(frozen["directory"]).glob("smoke_passed.*.json")):
        fail(f"no passed local smoke for the frozen programs {frozen['key']}; run --local-smoke first")
    stamp = utc("%Y%m%dT%H%M%SZ")
    run_dir = PRODUCTION / "condor" / args.stage / stamp
    for sub in ("tasks", "logs", "records"):
        (run_dir / sub).mkdir(parents=True, exist_ok=False)
    program_name = Path(frozen["files"]["program"]["frozen"]).name
    rows = []
    for task in pending:
        payload = dict(task)
        payload["record_dir"] = str(run_dir / "records")
        payload["frozen_program"] = frozen["files"]["program"]
        payload["program_timeout_s"] = PROGRAM_TIMEOUT[args.flavour]
        task_path = run_dir / "tasks" / f"{task['task_id']}.json"
        task_path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        rows.append(f"{task['task_id']} {task_path} {task_path.name} {program_name}\n")
    queue = run_dir / "queue.txt"
    queue.write_text("".join(rows), encoding="utf-8")
    jdl = run_dir / "job.jdl"
    write_jdl(jdl, run_dir, frozen, args, queue)
    staging = {
        "schema": "h4l_v3_staging/2",
        "stage": args.stage,
        "created_utc": stamp,
        "plan": str(plan_path),
        "plan_sha256": sha256(plan_path),
        "frozen": frozen,
        "resources": {"request_memory_mb": args.memory_mb, "request_disk_mb": args.disk_mb, "job_flavour": args.flavour},
        "tasks_total": len(tasks),
        "tasks_complete": complete,
        "tasks_in_flight": sorted(flying),
        "tasks_staged": len(pending),
        "run_dir": str(run_dir),
        "jdl": str(jdl),
        "queue": str(queue),
    }
    (run_dir / "staging.json").write_text(json.dumps(staging, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[staged] {run_dir}")

    if args.local_smoke:
        if args.local_smoke not in {task["task_id"] for task in pending}:
            fail(f"task {args.local_smoke} is not pending in this plan")
        return local_smoke(run_dir, frozen, run_dir / "tasks" / f"{args.local_smoke}.json")
    if cluster_mode:
        dry = run_dir / "dry_run.classads"
        result = subprocess.run(["condor_submit", "-dry-run", str(dry), str(jdl)], cwd=run_dir, capture_output=True, text=True)
        print(result.stdout.strip())
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            fail("condor_submit -dry-run failed")
        check_dry_run(dry.read_text(encoding="utf-8"), run_dir, pending, program_name)
        print(f"[dry-run] ok: {len(pending)} jobs, Iwd and logs in {run_dir}")
        if args.submit:
            result = subprocess.run(["condor_submit", str(jdl)], cwd=run_dir, capture_output=True, text=True)
            print(result.stdout.strip())
            if result.returncode != 0:
                print(result.stderr, file=sys.stderr)
                fail("condor_submit failed")
            match = re.search(r"(\d+) job\(s\) submitted to cluster (\d+)", result.stdout)
            if match:
                # Recorded whenever a cluster exists, so that the in-flight check sees its jobs.
                submission = {"cluster": int(match.group(2)), "jobs": int(match.group(1)), "submitted_utc": utc()}
                (run_dir / "submission.json").write_text(json.dumps(submission, indent=1) + "\n", encoding="utf-8")
            if not match or int(match.group(1)) != len(pending):
                fail("could not confirm the number of submitted jobs")
            print(f"[submitted] cluster {submission['cluster']} with {submission['jobs']} jobs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
