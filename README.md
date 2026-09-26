# H->ZZ*->4l reference analysis on the UL16 PFNano v3 samples

This directory holds the reference analysis of the SlopBench H->4l benchmark.
It applies the method of CMS AN-16-442 v8 and JHEP 11 (2017) 047 to the
20 fb^-1 UL16 PFNano v3 pseudo-data and to the matching MC.

- The durable rules are in `AGENTS.md`.
- The living plan and the decisions are in `PLAN.md`.
- The dated record is in `experiment_log.md`.
- The 2017 Open Data workflow and its README are historical; they are
  preserved in `deliverables/pre_v3_rules_20260923/README.md`.

## Inputs

| input | location |
|---|---|
| pseudo-data | `/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/pseudo_data/917f9326efa44ea9/part_000{00..82}.root` (read-only) |
| MC | `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l_seeds_v2/mc/<sample>/` (read-only) |
| YR4 cross sections | `/eos/home-y/yiyangz/codex/slopbench_code_fork/results/h4l_seeds/v3/inputs/Higgs_XSBR_YR4_update.xlsx` |

The data path and the luminosity are parameters of the workflow, so another
dataset of the same family runs without code changes.

## Layout

| path | content |
|---|---|
| `analysis_v3/` | all v3 code (C++ scans and fits, Python plots and orchestration), one folder per stage |
| `production_v3/` | every generated product and every Condor artifact (git-ignored) |
| `deliverables/` | documents: the analysis note and the preserved 2017 rules |
| `production_framework/`, `calibration_skim/`, `h4l_reconstruction/`, `reducible_background/`, `inference/`, `fake_data_nanoaod/`, `ntuplizer/` | 2017 code, reference only |

## Environment

Run everything from the repository root:

    cd /eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l
    export TMPDIR=$PWD/production_v3/tmp PIXI_CACHE_DIR=$PWD/production_v3/tmp/pixi_cache
    pixi install

- C++ programs build against the system ROOT (`/usr/bin/root-config`).
- Python scripts run with `pixi run py -- <script> [args]`.
- For Condor, see `AGENTS.md` ("Before every eossubmit submission").

## Workflow

The stages follow `PLAN.md` section 3.  Each stage below gets its runnable
commands, working directory, completion gate and plot directory when it is
implemented.

| stage | status |
|---|---|
| 0. preservation and rules | done |
| 1. infrastructure: catalogue, manifests, Condor framework | done (manifests v2) |
| 2. object definitions and skims | done (skims v2) |
| 3. lepton calibration and tag-and-probe | done (calibration v4/nominal, lambda v4, T&P v7 fits_i1c; closures) |
| 4. 4l reconstruction and optimization | done (h4l_select v5; N-1 selections v5_nm1d_*) |
| 5. backgrounds | done (zx_v2) |
| 6. signal model | done (sm_v2, gen_v1) |
| 7. statistical model | done (models model_*_r1: efficiency-uncertainty split, renormalized Z+X fractions) |
| 8. results | done (production_v3/results/v5/r2; r1 superseded by the review fixes) |
| 9. validation | done (toys, paired injection, closures) |
| 10. freeze, unblind, AN | final reviews; then the truth question and the AN |

### Stage 1. Input manifests

Working directory: the repository root.

1. Build the framework programs:

       make -C analysis_v3/framework

2. Make the validation plan.  This is the only directory listing in the
   workflow; choose a new version for a new plan:

       pixi run py -- analysis_v3/framework/scripts/build_manifests.py plan --version v2

3. Run one task locally through the frozen worker with a minimal `PATH`
   (smoke test):

       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/manifests/v2/plan.json \
           --stage manifests_v2 --program analysis_v3/framework/bin/validate_inputs \
           --memory-mb 1000 --disk-mb 1000 --flavour espresso --local-smoke validate_ZZTo4L_0000

4. Submit the remaining tasks:

       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/manifests/v2/plan.json \
           --stage manifests_v2 --program analysis_v3/framework/bin/validate_inputs \
           --memory-mb 1000 --disk-mb 1000 --flavour espresso --submit

   Tasks whose output already exists are skipped, so re-running the command
   resubmits only missing tasks.  The staging directory is printed as
   `[staged]`, and the cluster ID is in its `submission.json`.

5. Completion gate: every task output exists in
   `production_v3/manifests/v2/validation/`.  Then assemble the manifests
   (the full output scan of this stage) and plot:

       pixi run py -- analysis_v3/framework/scripts/build_manifests.py assemble --version v2
       pixi run py -- analysis_v3/framework/scripts/plot_manifests.py --version v2 --plot-dir plots_r2

The manifests are in `production_v3/manifests/v2/` (`data.json`,
`mc_<sample>.json`, `summary.json`, `summary.md`), and the plots are in
`production_v3/manifests/v2/plots_r2/`.

### Stage 2. Skims

Working directory: the repository root.  Inputs: the stage-1 manifests.

1. Build and plan (the analysis configuration
   `analysis_v3/common/config/analysis_ul16_v3.json` is frozen into every task):

       make -C analysis_v3/skims
       pixi run py -- analysis_v3/skims/scripts/make_skim_plan.py --manifests v2 --version v2

2. Smoke one signal task and one pseudo-data shard, then submit the rest:

       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/skims/v2/plan.json \
           --stage skims_v2 --program analysis_v3/skims/bin/skim_v3 --memory-mb 2000 --disk-mb 3000 \
           --flavour microcentury --local-smoke skim_VHToZZ_M125_0000
       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/skims/v2/plan.json \
           --stage skims_v2 --program analysis_v3/skims/bin/skim_v3 --memory-mb 2000 --disk-mb 3000 \
           --flavour microcentury --submit

3. Completion gate: all 283 task outputs exist.  Then run the full output
   scan, the coverage check, the histograms and the plots:

       pixi run py -- analysis_v3/skims/scripts/validate_skims.py --version v2 --manifests v2 --plot-dir plots_r2

The skims are in `production_v3/skims/v2/<sample>/`, with trees `Pairs`,
`PhotonPairs`, `TrigObjPairs`, `Events` (slim NanoAOD with the original
`h4l_file_key` and `h4l_entry`), `GenTable` (signal) and `Files`.  The validation output is in
`production_v3/skims/v2/validation/`, and the plots are in
`production_v3/skims/v2/plots_r2/`.


### Stage 3a. Lepton momentum scale and resolution

Working directory: the repository root.  Inputs: the stage-2 skims (v2) and
the stage-1 manifests (v2).  Method: `PLAN.md` stage 3a (factorized per-lepton
model, event-level MC template fits, iterated); configurations
`analysis_v3/calibration/config/calibration_ul16_v4.json` (calib_v4: the
template-fit fix of 2026-09-25, see experiment_log.md) and
`closure_ul16_v2.json`.  Current nominal run: `v4/nominal`.

1. Build, plan and smoke the Condor extraction of the calibration pairs
   (43 tasks), then submit the rest:

       make -C analysis_v3/calibration
       pixi run py -- analysis_v3/calibration/scripts/make_extract_plan.py --skims v2 --version v2
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/calibration/v2/extract/plan.json \
           --stage calib_extract_v2 --program analysis_v3/calibration/bin/calib_extract --memory-mb 2000 --disk-mb 1000 \
           --flavour microcentury --local-smoke calx_TTBar_0000
       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/calibration/v2/extract/plan.json \
           --stage calib_extract_v2 --program analysis_v3/calibration/bin/calib_extract --memory-mb 2000 --disk-mb 1000 \
           --flavour microcentury --submit

2. Completion gate: all 43 outputs in `production_v3/calibration/v2/extract/`.
   Full scan (plan binding, sha256, root_check, row counts, coverage):

       pixi run py -- analysis_v3/calibration/scripts/scan_extract.py --version v2 --manifests v2

3. The iterative calibration (local, about 1.5 h; run it in tmux).  It is
   restartable: finished jobs are reused only if their job is identical.

       pixi run py -- analysis_v3/calibration/scripts/run_calibration.py --version v4 --extract v2 --run nominal

   Completion gate: `production_v3/calibration/v4/<name>/summary.json` with
   `"converged": true` and `payload.json`; a run that reached the iteration
   limit at the statistical plateau is finalized with the mean of the payloads
   of its last iterations (their rms added to the uncertainty):

       pixi run py -- analysis_v3/calibration/scripts/finalize_calibration.py --version v4 --run nominal --last 4

   Then the plots and the report values (s and r at 45 GeV, |eta| 1.2):

       pixi run py -- analysis_v3/calibration/scripts/plot_calibration.py --version v4 --run nominal

   Plots: `production_v3/calibration/v3/<name>/plots/`.  The charge,
   signed-eta and phi diagnostics of a payload:

       pixi run py -- analysis_v3/calibration/scripts/diagnose_calibration.py --version v3 --extract v2 --run <name> \
           --payload-iteration <n> --label diag_iter<n+1>

4. Closure runs (`--closure <scenario>` of `closure_ul16_v1.json`) and their
   report (`<run>/closure/`):

       pixi run py -- analysis_v3/calibration/scripts/run_calibration.py --version v4 --extract v2 \
           --run closure_<scenario> --closure <scenario> \
           --closure-config analysis_v3/calibration/config/closure_ul16_v2.json --workers 5
       pixi run py -- analysis_v3/calibration/scripts/closure_report.py --version v4 --run closure_<scenario> \
           --closure-config analysis_v3/calibration/config/closure_ul16_v2.json

5. The final-payload diagnostics (charge, signed eta, phi; FullPairs for
   stage 3b):

       pixi run py -- analysis_v3/calibration/scripts/diagnose_calibration.py --version v4 --extract v2 --run nominal \
           --final --label diag_final

### Stage 3b. Per-event mass uncertainty (lambda)

Method: `PLAN.md` stage 3b; configuration
`analysis_v3/calibration/config/lambda_ul16_v2.json` (nominal; the core-window
variant `lambda_ul16_v2b.json` gives the method uncertainty).  Local, a few
minutes:

    make -C analysis_v3/calibration bin/lambda_histograms bin/fit_lambda
    pixi run py -- analysis_v3/calibration/scripts/run_lambda.py --version v4 --run nominal --label lambda_v2_r2
    pixi run py -- analysis_v3/calibration/scripts/run_lambda.py --version v4 --run nominal --label lambda_v2b_r3 \
        --config analysis_v3/calibration/config/lambda_ul16_v2b.json

Outputs: `production_v3/calibration/v4/nominal/<label>/lambda.json` and
`plots/` (closure per region, fit galleries).  (The v4 label lambda_v2 is the
superseded first attempt, stopped by two non-converged electron fits before
the fit_lambda freezing and multi-start fix.)

### Stage 3c. Tag-and-probe

Configuration `analysis_v3/tnp/config/tnp_ul16_v7.json` (the extraction
section is shared with v2-v6; model: DY-MC template (x) Gaussian with adaptive
template smoothing and a nearest-neighbour floor, the fail template mixed with
a free pass-like fraction, CMSShape background with the EGM ranges,
effective-weight errors for weighted MC); the efficiency chain is id | loose,
sip | id, iso | id and sip (conditional), the full step is measured directly as
a cross-check.  Current run: `v7` (histograms with the calibration v4 payload).

1. Condor extraction (43 tasks) and its full scan:

       make -C analysis_v3/tnp
       pixi run py -- analysis_v3/tnp/scripts/make_tnp_plan.py --skims v2 --version v2
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/tnp/v2/extract/plan.json \
           --stage tnp_extract_v2 --program analysis_v3/tnp/bin/tnp_extract --memory-mb 1000 --disk-mb 2000 \
           --flavour espresso --local-smoke tnpx_TTBar_0000
       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py ... --submit
       pixi run py -- analysis_v3/tnp/scripts/scan_tnp_extract.py --version v2 --manifests v2

2. Histograms (about 2 min) and one fit round per inspection label (about
   1 h, in tmux; the histograms are reused by later labels of the same run),
   then the plots, the gallery pages and the report (SF at 45 GeV, |eta| 1.2):

       pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal --run v6 --fit-label fits_i1
       pixi run py -- analysis_v3/tnp/scripts/plot_tnp.py --extract v2 --run v6 --fit-label fits_i1

   The fits run as Condor chunks of 4 bins (each refits the inclusive bin
   first, so the chunks reproduce the serial fits exactly), then merge:

       pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal --run v7 \
           --fit-label fits_i1c --condor plan
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/tnp/v2/v7/fits_i1c/condor_plan.json \
           --stage tnp_fits_v7_fits_i1c --program analysis_v3/tnp/bin/fit_tnp --memory-mb 2000 --disk-mb 1000 \
           --flavour microcentury --local-smoke fit_mm_iso_c046          # then --dry-run, --submit
       pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal --run v7 \
           --fit-label fits_i1c --condor merge
       pixi run py -- analysis_v3/tnp/scripts/plot_tnp.py --extract v2 --run v7 --fit-label fits_i1c --workers 16

   Every gallery page (`production_v3/tnp/v2/v7/<label>/plots/gallery/` and
   `gallery_alternatives/`) is inspected; later rounds refit and re-inspect
   only the bins that are not fine.  FINAL: run v7, label fits_i1c (user
   decision 2026-09-25).  Superseded:
   runs `nominal_r2` (v2 binning), `v3` (shared BW (x) CB model, electron
   WP90; fits_i1 stopped), `v4` (analytic model), `v5` (templates with a fixed
   0.5 GeV smoothing, calibration v3 payload), `v6` (v6 model: no fail
   admixture, sandwich errors).  An existing run keeps the frozen
   tnp_histograms recorded in its run.json.

3. Efficiency closure (thinned MC; MC halves by original-file-key parity, the
   even half with a known identification loss keep_id per flavour plays the
   data; coarse bins; three points; about 1.5 h per point):

       make -C analysis_v3/tnp bin/tnp_histograms
       for k in k1 k2 k3; do pixi run py -- analysis_v3/tnp/scripts/run_tnp.py --extract v2 --calibration v4/nominal \
           --run closure_v1_$k --closure analysis_v3/tnp/config/tnp_closure_ul16_v1.json --closure-point $k \
           --fit-label fits_i1 --workers 4; done

   Outputs: `production_v3/tnp/v2/closure_v1_<k>/fits_i1/` (sf.json, closure.json).

### Stage 4a. 4-lepton reconstruction and selection

1. Refit inputs (true Z1 line shape, FSR photon resolution):

       pixi run py -- analysis_v3/reconstruction/scripts/make_refit_inputs.py --skims v2 --version v3

2. Event records (h4l_reco, Condor): plan, smokes (one data, one signal and
   one MC task), dry run, submission, full scan:

       make -C analysis_v3/reconstruction
       pixi run py -- analysis_v3/reconstruction/scripts/make_reco_plan.py --skims v2 --calibration v3/nominal_r3 --version v2
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/h4l_reco/v2/plan.json \
           --stage h4l_reco_v2 --program analysis_v3/reconstruction/bin/h4l_reco --memory-mb 2000 --disk-mb 2000 \
           --flavour microcentury --local-smoke <task_id>
       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py ... --dry-run   # then --submit
       pixi run py -- analysis_v3/reconstruction/scripts/scan_reco.py --version v2 --manifests v2

3. Final selection (h4l_select, Condor, MELA; the plan records the sha256 of
   every MELA library, dictionary and data file, the working copies and the
   MELA self-test reference, all verified by every task), smokes (a data and
   an MC task), dry run, submission, full scan:

       make -C analysis_v3/reconstruction
       pixi run py -- analysis_v3/reconstruction/scripts/make_select_plan.py --reco v2 --version v5 \
           --calibration v4/nominal --lambda-run v4/nominal/lambda_v2_r2 --refit v3
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/h4l_select/v5/plan.json \
           --stage h4l_select_v5 --program analysis_v3/reconstruction/bin/h4l_select --memory-mb 2000 --disk-mb 2000 \
           --flavour longlunch --local-smoke sel_pseudo_data_0000     # and sel_GluGluToHToZZ_M125_0000
       source analysis_v3/framework/scripts/eossubmit_env.sh
       pixi run py -- analysis_v3/framework/scripts/stage_condor.py ... --dry-run   # then --submit
       pixi run py -- analysis_v3/reconstruction/scripts/scan_select.py --version v5 --manifests v2

   Versions: v1 and v2 superseded (v2 run incomplete: one TTBar muon with a
   NaN momentum error); v3 = SR Z2 > 12 GeV (reference of the stage-4b
   scans); v4 = the stage-4b selection (SR Z2 > 10 GeV) on the event-record
   calibration (v3 payload); v5 = the final selection, every lepton
   recalibrated from its raw pT with the calibration v4 payload and lambda v4.

   The v1 plan of h4l_reco (production_v3/h4l_reco/v1, condor/h4l_reco_v1)
   is SUPERSEDED (use-after-free fixed, 256 candidates); it must not be
   submitted.

### Stage 4b. Optimization on MC

    pixi run py -- analysis_v3/reconstruction/scripts/optimize_selection.py --select v3 --label scan_i1

Outputs: `production_v3/optimization/v3/<label>/` (scan.json, scan_<variable>.png).

N-1 expected significances of the declared cuts (eval_selfreport.json): seven
N-1 selections of the MC (one requirement relaxed to the event-record floor
each; the same frozen h4l_select, event records, payloads and MELA block as
v5; do not run make in analysis_v3/reconstruction before staging, it would
rebuild h4l_select), then the MC counting significance in 118-130 GeV:

    pixi run py -- analysis_v3/reconstruction/scripts/make_nm1_configs.py
    for cut in iso sip z1 z2 pt leadpt osmass; do
      pixi run py -- analysis_v3/reconstruction/scripts/make_select_plan.py --reco v2 --version v5_nm1_$cut \
          --calibration v4/nominal --lambda-run v4/nominal/lambda_v2_r2 --refit v3 \
          --selection-config analysis_v3/reconstruction/config/nm1/selection_ul16_v2_nm1_$cut.json \
          --samples GluGluToHToZZ_M125 VBF_HToZZ_M125 VHToZZ_M125 ZZTo4L GGZZ4Mu GGZZ4E GGZZ2E2Mu DYJetsToLL TTBar
      pixi run py -- analysis_v3/framework/scripts/stage_condor.py --plan production_v3/h4l_select/v5_nm1_$cut/plan.json \
          --stage h4l_select_v5_nm1_$cut --program analysis_v3/reconstruction/bin/h4l_select --memory-mb 2000 \
          --disk-mb 2000 --flavour longlunch --local-smoke sel_GGZZ4E_0000     # then --dry-run, --submit
    done
    pixi run py -- analysis_v3/reconstruction/scripts/scan_select.py --version v5_nm1_<cut> --manifests v2
    pixi run py -- analysis_v3/reconstruction/scripts/nm1_significance.py --nominal v5 \
        --variants iso sip z1 z2 pt leadpt osmass --label nm1_v1

Outputs: `production_v3/optimization/v5/nm1_v1/` (nm1.json, nm1.csv, nm1.md).
Categories: `analysis_v3/reconstruction/config/categories_ul16_v2.json` (paper-convention production
discriminants recomputed from the stored raw MELA probabilities; constants and working points optimized
on MC in `production_v3/tmp/category_opt/`, adopted by the user on 2026-09-25; see experiment_log.md).

### Stage 5. Backgrounds

    pixi run py -- analysis_v3/backgrounds/scripts/zx_estimate.py --select <v> --label <zx label>

Z+X OS and SS methods, the MC-closure systematic (user decision 2026-09-25), combination and Landau
shapes: `production_v3/backgrounds/<v>/<label>/` (zx.json, plots/).  The detailed closure tests
(MC closure per final state and m4l window, 3P1F predicted from 2P2F in data and MC):

    pixi run py -- analysis_v3/backgrounds/scripts/zx_closure.py --select <v> --label <closure label>

qqZZ and ggZZ come from the MC inside the model builder (stage 7).

### Stage 6. Signal model

    pixi run py -- analysis_v3/signal_model/scripts/yr4_xsbr.py --version v1
    pixi run py -- analysis_v3/signal_model/scripts/signal_model.py --select <v> --yr4 v1 --label <signal label>

YR4 sigma x BR table (`production_v3/signal_model/yr4/v1/`), DCB shapes, per-event width model
(conditional on D_mass = dm/m4l: width s x D_mass x m_H), VH non-resonant part, rest-frame morphing and
acceptance ratios (`production_v3/signal_model/<v>/<label>/`; current sm_v2, sm_v1 superseded by the
D_mass-conditional width).  SM fiducial cross sections from the skim generator tables:

    pixi run py -- analysis_v3/signal_model/scripts/fiducial_xsec.py --select <v> --label gen_v1

### Stage 7. Statistical model and fits

    pixi run py -- analysis_v3/inference/scripts/build_model.py --select <v> --signal <signal label> --zx <zx label> \
        --sf production_v3/tnp/v2/<run>/<fit label>/sf.json --label <model label>
    pixi run py -- analysis_v3/inference/scripts/fit_model.py --model <v>/<model label> --yr4 v1 --label <fit label> \
        [--dimension 1D|2D|3D] [--no-refit] [--set a|b] [--dataset data|asimov] \
        [--poi-scheme inclusive|final_state|category|fv|mode|stxs0|fid_fs|fid_int|fid_<obs>] [--fix-mh 125.09] \
        [--final-state 4mu|4e|2e2mu] [--fiducial production_v3/signal_model/<v>/gen_v1/fiducial_xsec.json]
    pixi run py -- analysis_v3/inference/scripts/scan_likelihood.py --model <v>/<model label> --yr4 v1 --label <scan label> \
        --x mu 0 2.5 26 [--y <name> <low> <high> <n>] [--poi-scheme fv] [--stat-only] [--dataset asimov]
    pixi run py -- analysis_v3/inference/scripts/fit_width.py --model <v>/<model label> --yr4 v1 --label <width label>
    pixi run py -- analysis_v3/inference/scripts/z4l_mass.py --select <v> --label <z4l label>
    pixi run py -- analysis_v3/inference/scripts/gof.py --model <v>/<model label> --fit <fit label> --yr4 v1 --label <gof label>
    pixi run py -- analysis_v3/inference/scripts/make_pyhf.py --model <v>/<model label> --fit <fit label> --yr4 v1 --out MODEL.json

Models: `--channels categories` (default; the mass, signal-strength and mode fits), `inclusive` (per final
state, the fiducial fits: 1D m4l without the refit, as in the paper), `pt4l`, `njets`, `ptj1` (differential
fiducial fits).

### Stage 8. Results

    pixi run py -- analysis_v3/inference/scripts/plot_distributions.py --select v5 --zx zx_v2 \
        --sf production_v3/tnp/v2/v7/fits_i1c/sf.json --label dist_r1
    production_v3/tmp/build_models.sh r1 production_v3/tnp/v2/v7/fits_i1c/sf.json     # the five final models
    pixi run py -- analysis_v3/inference/scripts/run_results.py --select v5 --tag r2 --model-cat model_cat_r1 \
        --model-incl model_incl_r1 --model-pt4l model_pt4l_r1 --model-njets model_njets_r1 --model-ptj1 model_ptj1_r1 \
        --fiducial production_v3/signal_model/v5/gen_v1/fiducial_xsec.json --yr4 v1 --workers 12
    pixi run py -- analysis_v3/inference/scripts/validate_dmass_model.py --select v5 --model model_cat_r1 --label dmass_v1
    pixi run py -- analysis_v3/inference/scripts/validate_mc_closure.py --select v5 --model model_cat_r1 --yr4 v1 --label closure_mc_v1
    pixi run py -- analysis_v3/inference/scripts/make_pyhf.py --model v5/model_cat_r1 --fit r2_b3d_data --yr4 v1 \
        --out production_v3/results/v5/<out tag>/MODEL.json
    pixi run py -- analysis_v3/inference/scripts/make_results.py <the run_results arguments> --calibration v4/nominal \
        --sf production_v3/tnp/v2/v7/fits_i1c/sf.json --nm1 production_v3/optimization/v5/nm1_v1/nm1.json \
        --pyhf production_v3/results/v5/r1/MODEL.json

Outputs: `production_v3/results/v5/<out tag>/` (RESULT.json, MODEL.json, eval_selfreport.json, results_full.json,
summary.md, key_numbers.csv, cuts.csv, systematics_coverage.csv, plots/paper_comparison.png).

### Stage 9. Validation with pseudo-experiments

    pixi run py -- analysis_v3/inference/scripts/run_toys.py --model <v>/<model label> --yr4 v1 --label <toy label> \
        --mu 0 0.5 1 2 3 --mh 121 123 125 127 129 --ntoys 500 --workers 8

Outputs: `production_v3/inference/<v>/<model label>/<label>/` (fit.json, toys.json).  Paired injection
(`--paired`: dataset B = dataset A + an independent mu = 1 signal, every background event shared).

### Stage 10. Analysis note

    pixi run py -- analysis_v3/documentation/make_an.py --select v5 --tag r1 --model-cat model_cat_r0 --tnp v7/fits_i1c \
        --calibration v4/nominal --zx zx_v2 --signal sm_v2 --dist dist_r1 --out deliverables/AN_h4l_ul16_pfnano_v1 \
        --toys v5/model_cat_r0/r1_toys_mu v5/model_cat_r0/r1_toys_mh --paired v5/model_cat_r0/r1_paired \
        --tnp-closure production_v3/tnp/v2/closure_v1_k1/fits_i1 ... [--truth <truth json after unblinding>]

### Stage 11. Evaluation submission (EVAL_CONTRACT v0.4, task h4l_ntuple)

    pixi run py -- analysis_v3/eval/make_eval_constants.py          # my_analysis/h4l_eval/data/constants.json
    pixi run py -- analysis_v3/eval/make_test_dataset.py --name test_10fb --data-shards 41            # DATASET.md layout
    pixi run py -- analysis_v3/eval/make_test_dataset.py --name test_10fb_b --data-shards 42 --data-offset 41
    pixi run py -- analysis_v3/eval/make_test_dataset.py --name full_20fb --data-shards 83
    my_analysis/run.sh production_v3/eval_datasets/<name> <output_dir>    # RESULT.json, MODEL.json, ...
    python tests/check_outputs.py <output_dir>                           # contract checks, pyhf build and refit

`my_analysis/` is self-contained (its own pixi.toml / pixi.lock with root_base and the C++ compiler; run.sh uses an
active environment with ROOT and the Python packages, else `pixi run --frozen`).  It follows the main analysis in every
step except the simplified tag-and-probe nominal fit (user rule 2026-09-26); see my_analysis/README.md.  Development
only: H4L_DEBUG_DUMP=<file> saves the compact reader outputs, H4L_DEBUG_LOAD=<file> reuses them.
