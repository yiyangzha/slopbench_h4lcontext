#!/usr/bin/env bash
# Evaluation entry point from the repository root: forwards to the submission directory my_analysis/
# (EVAL_CONTRACT.md: ./run.sh <dataset_dir> <output_dir>).
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/my_analysis/run.sh" "$@"
