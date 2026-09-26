#!/usr/bin/env bash
# Generic v3 worker: run one frozen C++ program on one task and publish its
# outputs atomically on EOS.
#
# Usage: run_task.sh TASK_JSON PROGRAM
#
# The program runs in the local sandbox as
#   PROGRAM --task TASK_JSON --out-json LOCAL.json [--out-root LOCAL.root]
# and must validate its own outputs before returning 0.  A ROOT output must
# be described in the JSON output as {"root_output": {"trees": {NAME: N},
# "objects": [...]}} with at least one tree.  The wrapper then
#   * copies each output to an attempt-qualified temporary name beside its
#     final EOS path and verifies the copy (sha256, and root_check for ROOT);
#   * adds a "publication" block (attempt, ROOT path and sha256) to the JSON;
#   * renames the ROOT, then the JSON, onto their final paths without ever
#     overwriting.  The JSON marks a complete task: a final ROOT without its
#     JSON is an orphan of an interrupted attempt, which is moved aside under
#     an attempt-qualified name (never deleted) before publication.
# Every attempt leaves a record with its result, exit code, message and the
# tails of the program output in the task's record directory.
# Every EOS copy, checksum and ROOT check, and the program itself (the task's
# program_timeout_s, set by the stager), runs under a time limit: a stalled
# XRootD or FUSE read (stage-2 job 1152510.193) fails the attempt instead of
# holding the slot until the flavour limit.  A job terminated by the batch
# system (SIGTERM) still leaves its attempt record.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 TASK_JSON PROGRAM" >&2
  exit 64
fi
TASK=$1
PROGRAM=$2
[[ "${PROGRAM}" == */* ]] || PROGRAM="./${PROGRAM}"
ROOT_CHECK=./root_check
STARTS=""
if [[ -n "${_CONDOR_JOB_AD:-}" && -r "${_CONDOR_JOB_AD}" ]]; then
  STARTS=$(sed -n 's/^NumJobStarts = //p' "${_CONDOR_JOB_AD}")
fi
ATTEMPT="${H4L_CLUSTER:-local}.${H4L_PROC:-0}.s${STARTS:-0}.$(date -u +%Y%m%dT%H%M%SZ).$(hostname -s 2>/dev/null || echo host).$$.${RANDOM}${RANDOM}"
STARTED="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
MESSAGE=""

[[ -r "${TASK}" ]] || { echo "ERROR: task JSON unreadable: ${TASK}" >&2; exit 66; }
[[ -x "${PROGRAM}" ]] || { echo "ERROR: program not executable: ${PROGRAM}" >&2; exit 66; }
command -v python3 >/dev/null || { echo "ERROR: python3 is required" >&2; exit 66; }

FIELDS_TEXT=$(python3 - "${TASK}" <<'PY'
import json, sys
task = json.load(open(sys.argv[1]))
outputs = task["outputs"]
print(task["task_id"])
print(outputs["json"])
print(outputs.get("root") or "")
print(task["record_dir"])
print(int(task.get("program_timeout_s", 0)))
PY
) || { echo "ERROR: cannot read the task fields of ${TASK}" >&2; exit 66; }
readarray -t FIELDS <<< "${FIELDS_TEXT}"
TASK_ID=${FIELDS[0]}
OUT_JSON=${FIELDS[1]}
OUT_ROOT=${FIELDS[2]}
RECORD_DIR=${FIELDS[3]}
PROGRAM_TIMEOUT=${FIELDS[4]}
[[ "${PROGRAM_TIMEOUT}" =~ ^[0-9]+$ && "${PROGRAM_TIMEOUT}" -gt 0 ]] ||
  { echo "ERROR: the task has no positive program_timeout_s" >&2; exit 64; }
for path in "${OUT_JSON}" "${RECORD_DIR}" ${OUT_ROOT:+"${OUT_ROOT}"}; do
  [[ "${path}" == /eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/* ]] ||
    { echo "ERROR: path outside production_v3: ${path}" >&2; exit 64; }
done
LOG="./${TASK_ID}.program"

make_directory() {
  # The stager creates the output directories before submission; this only
  # covers the record directory and late cases.  EOS FUSE can hide a
  # directory created concurrently by another worker for a while.
  local directory=$1 attempt
  for attempt in 1 2 3 4 5 6 7 8; do
    mkdir -p "${directory}" 2>/dev/null || true
    [[ -d "${directory}" ]] && return 0
    sleep $((5 * attempt))
  done
  echo "ERROR: cannot create ${directory}" >&2
  return 1
}

record() {
  local result=$1 code=$2
  make_directory "${RECORD_DIR}" || return 0
  (
    set -o noclobber
    {
      printf 'task_id=%s\nattempt=%s\nresult=%s\nexit_code=%s\nmessage=%s\n' "${TASK_ID}" "${ATTEMPT}" "${result}" "${code}" "${MESSAGE}"
      printf 'started_utc=%s\nfinished_utc=%s\nhost=%s\n' "${STARTED}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname -f 2>/dev/null || hostname)"
      printf 'out_json=%s\nout_root=%s\n' "${OUT_JSON}" "${OUT_ROOT}"
      if [[ -f "${LOG}.err" ]]; then printf -- '--- program stderr (tail) ---\n'; tail -n 40 "${LOG}.err"; fi
      if [[ -f "${LOG}.out" ]]; then printf -- '--- program stdout (tail) ---\n'; tail -n 20 "${LOG}.out"; fi
    } > "${RECORD_DIR}/${TASK_ID}.${ATTEMPT}.txt"
  ) || echo "WARN: attempt record not written" >&2
}

fail() {
  local code=$1
  shift
  MESSAGE="$*"
  echo "ERROR: ${MESSAGE}" >&2
  record failed "${code}"
  exit "${code}"
}

trap 'MESSAGE="terminated by SIGTERM"; record terminated 143; exit 143' TERM

# Every shared library of the program must resolve on this worker (a loader failure would otherwise
# surface only as exit 127 of the program).
if command -v ldd >/dev/null; then
  UNRESOLVED=$(ldd "${PROGRAM}" 2>&1 | grep 'not found' || true)
  [[ -z "${UNRESOLVED}" ]] || fail 66 "unresolved shared libraries of ${PROGRAM}: ${UNRESOLVED//$'\n'/; }"
fi

if [[ -e "${OUT_JSON}" ]]; then
  echo "[task] ${TASK_ID}: final output exists, left untouched: ${OUT_JSON}"
  record skipped_existing 0
  exit 0
fi

LOCAL_JSON="./${TASK_ID}.out.json"
LOCAL_ROOT="./${TASK_ID}.out.root"
ARGS=(--task "${TASK}" --out-json "${LOCAL_JSON}")
[[ -n "${OUT_ROOT}" ]] && ARGS+=(--out-root "${LOCAL_ROOT}")
echo "[task] ${TASK_ID} attempt ${ATTEMPT}: ${PROGRAM} ${ARGS[*]}"
set +e
timeout --kill-after=60 "${PROGRAM_TIMEOUT}" "${PROGRAM}" "${ARGS[@]}" > "${LOG}.out" 2> "${LOG}.err"
code=$?
set -e
cat "${LOG}.out"
cat "${LOG}.err" >&2
[[ "${code}" -ne 124 ]] || fail 124 "program exceeded its ${PROGRAM_TIMEOUT} s limit"
[[ "${code}" -eq 0 ]] || fail "${code}" "program exited with ${code}"
[[ -s "${LOCAL_JSON}" ]] || fail 65 "program returned 0 without its JSON output"
python3 -c 'import json,sys; json.load(open(sys.argv[1]))' "${LOCAL_JSON}" || fail 65 "JSON output does not parse"

# Copy one local output beside its final path, verify the copy and rename it
# into place.
publish() {
  local local_file=$1 final=$2 kind=$3
  local temporary="${final}.partial.${ATTEMPT}.${kind}"
  make_directory "$(dirname "${final}")" || fail 73 "cannot create the output directory of ${final}"
  timeout 1800 cp -n "${local_file}" "${temporary}" || fail 74 "copy to ${temporary} failed or timed out"
  local local_sum copy_sum
  local_sum=$(sha256sum "${local_file}" | cut -d' ' -f1) || fail 74 "sha256 of ${local_file} failed"
  copy_sum=$(timeout 900 sha256sum "${temporary}" | cut -d' ' -f1) || fail 74 "sha256 of ${temporary} failed or timed out"
  [[ "${local_sum}" == "${copy_sum}" ]] || fail 74 "sha256 mismatch after copy: ${temporary}"
  if [[ "${kind}" == root ]]; then
    local spec
    spec=$(python3 - "${LOCAL_JSON}" <<'PY'
import json, sys
spec = json.load(open(sys.argv[1]))["root_output"]
trees = spec.get("trees", {})
if not trees:
    raise SystemExit("root_output lists no tree")
for name, entries in trees.items():
    print(f"{name}={int(entries)}")
for name in spec.get("objects", []):
    print(name)
PY
) || fail 65 "cannot read root_output from the JSON output"
    local checks attempt checked=0
    readarray -t checks <<< "${spec}"
    # ROOT reads /eos paths through XRootD, whose view of a file just written
    # through the FUSE mount can lag; retry before declaring the copy bad.
    for attempt in 1 2 3 4 5 6; do
      if timeout 600 "${ROOT_CHECK}" "${temporary}" "${checks[@]}"; then
        checked=1
        break
      fi
      echo "WARN: root_check attempt ${attempt} failed on ${temporary}; retrying" >&2
      sleep $((10 * attempt))
    done
    [[ "${checked}" -eq 1 ]] || fail 65 "root_check failed on the EOS copy ${temporary}"
    ROOT_SHA256=${local_sum}
  fi
  if [[ -e "${final}" ]]; then
    if [[ "${kind}" == json || -e "${OUT_JSON}" ]]; then
      MESSAGE="a complete competing output exists; kept ${temporary}"
      echo "WARN: ${MESSAGE}" >&2
      record skipped_competing 0
      exit 0
    fi
    # A ROOT without its JSON: an interrupted earlier attempt.  Keep it aside.
    mv -n "${final}" "${final}.orphan.${ATTEMPT}.root" || fail 74 "cannot move the orphan ${final} aside"
    [[ ! -e "${final}" ]] || fail 74 "the orphan ${final} is still in place"
    echo "WARN: moved an orphan ROOT aside: ${final}.orphan.${ATTEMPT}.root" >&2
  fi
  mv -n "${temporary}" "${final}" || fail 74 "rename onto ${final} failed"
  [[ -e "${final}" && ! -e "${temporary}" ]] || fail 74 "rename onto ${final} did not complete"
}

ROOT_SHA256=""
if [[ -n "${OUT_ROOT}" ]]; then
  [[ -x "${ROOT_CHECK}" ]] || fail 66 "root_check is required for a ROOT output"
  [[ -s "${LOCAL_ROOT}" ]] || fail 65 "program returned 0 without its ROOT output"
  publish "${LOCAL_ROOT}" "${OUT_ROOT}" root
fi
python3 - "${LOCAL_JSON}" "${ATTEMPT}" "${OUT_ROOT}" "${ROOT_SHA256}" <<'PY' || fail 65 "cannot add the publication block"
import datetime, json, sys
path, attempt, root_path, root_sha = sys.argv[1:5]
payload = json.load(open(path))
payload["publication"] = {"attempt": attempt, "root_path": root_path or None, "root_sha256": root_sha or None,
                          "published_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
with open(path + ".tmp", "x") as stream:
    json.dump(payload, stream, indent=1, sort_keys=True)
    stream.write("\n")
PY
mv "${LOCAL_JSON}.tmp" "${LOCAL_JSON}"
publish "${LOCAL_JSON}" "${OUT_JSON}" json
record completed 0
echo "[task] ${TASK_ID}: published ${OUT_JSON}"
