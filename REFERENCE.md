# Reference analysis: H -> ZZ* -> 4l on UL16 PFNano pseudo-data (branch `reference`)

## Evaluation submission (EVAL_CONTRACT v0.4, task h4l_ntuple)

    ./run.sh <dataset_dir> <output_dir>                  # from the repository root (forwards to my_analysis/run.sh)
    my_analysis/run.sh <dataset_dir> <output_dir>        # the submission directory itself
    tests/run_test.sh <dataset_dir> <output_dir>         # run + check_outputs.py (contract checks, pyhf build of MODEL.json)

`RESULT.json` and `MODEL.json` at the repository root are the scoring outputs of the main analysis on the 20 fb^-1
pseudo-data (identical to results/main_final/); results/eval_test_10fb/ holds the outputs of ./run.sh on the ~10 fb^-1
test dataset.

`my_analysis/` is the self-contained submission directory (run.sh, pixi.toml + pixi.lock, EVAL_CONTRACT.md with the
per-flavour calibration block, the code in h4l_eval/ and src/).  `<dataset_dir>` has the DATASET.md layout
(data/part_*.root, mc/<process>/*.root + meta.json, mc/cross_sections.json, lumi.json, meta.json) with ROOT files in
the UL16 PFNano format.  Outputs: RESULT.json, MODEL.json (pyhf workspace), MODEL_meta.json, summary.txt,
results_full.json, logs/run.log.  See my_analysis/README.md.

## Contents

* `my_analysis/` — the evaluation submission.
* `analysis_v3/` — the main analysis code (C++ programs with Makefiles, Python stages, configurations); the binaries are
  built with `make`.  `analysis_v3/eval/` builds the submission's method constants and test datasets.
* `results/main_final/` — the main analysis results on the 20 fb^-1 pseudo-data (RESULT.json in the contract format
  with the per-flavour calibration, MODEL.json, results_full.json, summary.md, tables, plots, the freeze record).
* `results/eval_test_10fb/` — the submission run on a ~10 fb^-1 test dataset in the DATASET.md layout (41 of the 83
  pseudo-data shards plus the MC; 18 min on 16 cores).
* `results/truth/` — the truth of the pseudo-data used for the comparison after unblinding.
* `deliverables/AN_h4l_ul16_pfnano_v2/` — the analysis note (main.pdf, main.tex; regenerated with
  analysis_v3/documentation/make_an.py).
* `AGENTS.md`, `PLAN.md`, `README.md`, `experiment_log.md` — rules, plan, workflow and the full decision log.

## Key numbers

Main analysis (20 fb^-1, set b, 3D fit with m_H floating): mu = 0.903 +0.189/-0.169, m_H = 124.96 +0.40/-0.39 GeV,
observed significance 8.5 (expected 10.1).  Truth: mu = 1, m_H = 125 GeV.  Calibration (data-weighted, muon /
electron): scale_shift -0.01993 / -0.01987 (truth -0.01994 / -0.01991), smear 0.0116 / 0.0113 (truth 0.010 in pT plus
0.01 angular smears), sel_eff 0.994 / 1.000 (truth 1 for this selection).
