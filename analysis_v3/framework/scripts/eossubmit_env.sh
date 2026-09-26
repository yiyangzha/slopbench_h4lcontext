#!/usr/bin/env bash
# Source before every eossubmit/Condor command, in the same shell:
#   source analysis_v3/framework/scripts/eossubmit_env.sh
# It keeps submit-side temporary state inside the repository, loads the
# eossubmit pool, sets the grid proxy and refuses a proxy with less than one
# day left.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Source this file: source ${BASH_SOURCE[0]}" >&2
  exit 64
fi
H4L_V3_REPO=/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l
: "${H4L_CONDOR_TMP:=${H4L_V3_REPO}/production_v3/tmp/condor_submit}"
mkdir -p "${H4L_CONDOR_TMP}"
[[ -d "${H4L_CONDOR_TMP}" && -w "${H4L_CONDOR_TMP}" ]] || {
  echo "ERROR: H4L_CONDOR_TMP is not writable: ${H4L_CONDOR_TMP}" >&2
  return 1
}
export H4L_CONDOR_TMP TMPDIR="${H4L_CONDOR_TMP}" TMP="${H4L_CONDOR_TMP}" TEMP="${H4L_CONDOR_TMP}"
export RANDFILE="${H4L_CONDOR_TMP}/.openssl-rand"
export EOS_MGM_URL=root://eosuser.cern.ch
export X509_USER_PROXY=/afs/cern.ch/user/y/yiyangz/private/x509up_u165165.backup_20260923
source /usr/share/Modules/init/bash
module load lxbatch/eossubmit || { echo "ERROR: module load lxbatch/eossubmit failed" >&2; return 1; }
H4L_PROXY_LEFT=$(voms-proxy-info -file "${X509_USER_PROXY}" -timeleft 2>/dev/null || echo 0)
if [[ "${H4L_PROXY_LEFT}" -lt 86400 ]]; then
  echo "ERROR: proxy ${X509_USER_PROXY} has ${H4L_PROXY_LEFT} s left (< 86400); ask the user to renew it" >&2
  return 1
fi
echo "[eossubmit_env] TMPDIR=${TMPDIR} proxy_left=${H4L_PROXY_LEFT}s pool=${_myschedd_POOL:-?}"
