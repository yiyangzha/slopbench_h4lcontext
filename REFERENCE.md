# Reference analysis: H -> ZZ* -> 4l on UL16 PFNano pseudo-data (branch `reference`)

## Evaluation submission (EVAL_CONTRACT v0.4, task h4l_ntuple)

    ./run.sh <dataset_dir> <output_dir>                  # from the repository root (forwards to my_analysis/run.sh)
    my_analysis/run.sh <dataset_dir> <output_dir>        # the submission directory itself
    tests/run_test.sh <dataset_dir> <output_dir>         # run + check_outputs.py (contract checks, pyhf build of MODEL.json)

`my_analysis/` is the self-contained submission (run.sh, pixi.toml + pixi.lock with root_base and the C++ compiler,
EVAL_CONTRACT.md with the per-flavour calibration block, the code in h4l_eval/ and src/).  It carries out the main
analysis of this repository with the same methods and parameters in every step; the only simplification is the
tag-and-probe nominal fit model (see my_analysis/README.md and the AN section "Evaluation submission").
`<dataset_dir>` has the DATASET.md layout (data/part_*.root, mc/<process>/*.root + meta.json,
mc/cross_sections.json, lumi.json, meta.json) with ROOT files in the UL16 PFNano format.  Outputs: RESULT.json,
MODEL.json (pyhf workspace), MODEL_meta.json, summary.txt, results_full.json, logs/run.log (and the intermediate
reader files in work/).

`RESULT.json` and `MODEL.json` at the repository root are the results of the main analysis on the 20 fb^-1
pseudo-data (identical to results/main_final/).

## Contents

* `my_analysis/` — the evaluation submission.
* `analysis_v3/` — the main analysis code (C++ programs with Makefiles, Python stages, configurations); the binaries
  are built with `make`.  `analysis_v3/eval/` builds the submission's method constants and the test datasets.
* `results/main_final/` — the main analysis on the 20 fb^-1 pseudo-data (RESULT.json with the per-flavour
  calibration, MODEL.json, results_full.json, summary.md, tables, plots, the freeze record).
* `results/eval_test_10fb_A/` — the submission run on half A of the pseudo-data (9.99 fb^-1 in the DATASET.md layout; 13 min on 16 cores with the input on EOS).
* `results/eval_test_10fb_B/` — the submission run on half B of the pseudo-data (10.01 fb^-1 in the DATASET.md layout; 17 min on 16 cores with the input on EOS).
* `results/eval_full_20fb/` — the submission run on the full 20 fb^-1 pseudo-data (20.00 fb^-1 in the DATASET.md layout; 25 min on 16 cores with the input on EOS).
* `results/truth/` — the truth of the pseudo-data (injection record and generator profiles, make_truth.py).
* `deliverables/AN_h4l_ul16_pfnano_v4/` — the analysis note (main.pdf, main.tex; regenerated with
  analysis_v3/documentation/make_an.py), including the truth comparison and the evaluation submission.
* `AGENTS.md`, `PLAN.md`, `README.md`, `experiment_log.md` — rules, plan, workflow and the full decision log.

## Key numbers

Truth: mu = 1, m_H = 125 GeV; scale_shift muon -0.01994, electron -0.01991; smear 0.01 (pT part; the injected angular smears are absorbed by
the measured smear); sel_eff muon 0.99352, electron 0.99998 (data-weighted).

| analysis | L [fb^-1] | mu | m_H [GeV] | Z obs (exp) | GoF p | scale_shift mu / e | smear mu / e | sel_eff mu / e |
|---|---|---|---|---|---|---|---|---|
| main analysis (3D, categories) | 20.00 | 0.903 +- 0.179 | 124.96 +- 0.39 | 8.46 (10.08) | 0.42 | -0.01993 / -0.01987 | 0.0116 / 0.0113 | 0.9941 / 1.0004 |
| submission, half A of the pseudo-data | 9.99 | 0.890 +- 0.246 | 125.41 +- 0.47 | 5.77 (6.48) | 0.04 | -0.01993 / -0.01992 | 0.0115 / 0.0113 | 0.9945 / 1.0023 |
| submission, half B of the pseudo-data | 10.01 | 0.813 +- 0.251 | 124.58 +- 0.55 | 5.08 (6.03) | 0.74 | -0.01993 / -0.01989 | 0.0116 / 0.0112 | 0.9945 / 1.0021 |
| submission, the full 20 fb^-1 pseudo-data | 20.00 | 0.833 +- 0.175 | 125.09 +- 0.38 | 7.47 (8.80) | 0.42 | -0.01994 / -0.01990 | 0.0114 / 0.0113 | 0.9945 / 1.0022 |

Main analysis: mu = 0.903 +- 0.179, m_H = 124.96 +- 0.39 GeV (set b, 3D fit with m_H floating); the submission fits m4l x D_mass binned without MELA categories.

