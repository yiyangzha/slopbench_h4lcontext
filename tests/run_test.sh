#!/usr/bin/env bash
# Run the evaluation submission on a dataset in the DATASET.md layout and check its outputs.
#   tests/run_test.sh <dataset_dir> <output_dir>
# A test dataset in that layout can be built from the analysis inputs on lxplus with
#   pixi run py -- analysis_v3/eval/make_test_dataset.py --name test_10fb --data-shards 41
# (symbolic links under production_v3/eval_datasets/<name>; about 10 fb^-1 with all the MC).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
start=$(date +%s)
"$HERE/run.sh" "$1" "$2"
echo "run.sh finished in $(( $(date +%s) - start )) s"
python "$HERE/tests/check_outputs.py" "$2"
