# H->ZZ*->4l renewal on the UL16 PFNano v3 samples: working plan

Version 2, 2026-09-23T22:10Z.  All decision points of version 1 are answered
(section 6); implementation starts with stage 0.

This is the living plan.  It is updated whenever the plan changes and re-read
before every stage, submission and report.  `experiment_log.md` keeps the
dated history (plans, decisions, questions and answers, submissions, checks,
results); `AGENTS.md` holds the durable rules and `README.md` the runnable
workflow.

## 1. Principles

- **Method.** AN-16-442 v8 and HIG-16-041 (JHEP 11 (2017) 047).
  - Where the inputs make an identical implementation impossible, use a
    reasonable, accurate and documented substitute.
  - Whenever several methods are possible, ask the user which to use.
- **Generality.** Nothing is specific to one dataset:
  - every detector effect is measured from the data at run time;
  - m_H always floats and mu is unbounded;
  - the data path and the luminosity are parameters.
- **Optimization.** Selection thresholds and lepton IDs start from the AN
  values and are optimized for the joint expected precision of mu and m_H at
  20 fb^-1, on MC only (never on the data).
- **Division of labour.**
  - C++ against system ROOT 6.40 reads every ROOT file and performs every
    fit (RooFit/RooStats/Minuit2).
  - Python in pixi handles orchestration, plots (mplhep) and the pyhf
    export.
- **Code and products.**
  - New code lives in `analysis_v3/`; the 2017 folders stay untouched.
  - Every product and every Condor artifact lives under `production_v3/`.
  - Condor is used only through eossubmit.
- **Integrity.**
  - Stage inputs come only from validated manifests.
  - Every ROOT output is closed, reopened and checked (keys, trees, entries)
    before its atomic publication.
  - Every multi-job stage gets a full output scan before downstream use.
  - Failed artifacts are kept under attempt-qualified names.
  - Staging is idempotent.
  - Every output records its input manifest, configuration sha and program
    sha.
- **Statistics.** Every data and MC entry is an independent event: Poisson
  data, ordinary MC statistics, no repetition correction.
- **Blinding.** The rules of prompt section 4.3 hold until the user unblinds.
- **Stage gates.** Each stage passes these gates in order:
  1. a local smoke test on a representative real input, with the worker
     wrapper and a minimal `PATH`;
  2. inspection of the outputs and plots;
  3. two independent reviews, one on physics and one on production and Condor
     safety;
  4. submission (the user allows submitting directly);
  5. monitoring: every few minutes until the queue is stable, then every
     15 minutes;
  6. a full output scan;
  7. numerical validation of the plots and a visual inspection;
  8. a report to the user: counts, validation, plot directories, issues.

## 2. Layouts

### 2.1 Code tree

```
analysis_v3/
  common/          shared C++ headers: objects, FSR, isolation, calibration
                   and SF lookup, seeded smearing, atomic ROOT output;
                   shared Python: plot style, JSON I/O;
                   config/analysis_ul16_v3.json, the frozen analysis
                   configuration
  framework/       sample catalogue, C++ input validator, manifest builder,
                   Condor stager, worker wrapper, output scanner, monitor,
                   eossubmit_env.sh
  skims/           the one-pass C++ skimmer: dilepton, multilepton and
                   signal-generator tables
  calibration/     Z-peak BW(x)DCB fits, per-lepton scale and smear,
                   D_mass lambda calibration
  tnp/             tag-and-probe: reconstruction, ID, isolation, SIP and
                   trigger legs; the 4l trigger method; SF payloads
  reconstruction/  4l selection, FSR, candidates, Z1 refit, D_mass, jets,
                   categories, discriminants (MELA), generator/fiducial/STXS;
                   selection optimization
  backgrounds/     qqZZ/ggZZ templates and theory variations; Z+X OS, SS and
                   combination
  signal_model/    DCB parameterization in m_H, morphing, YR4 sigma x BR,
                   A x eps(m_H), non-resonant VH
  inference/       likelihood, fits, scans, impacts, GoF, results; writers
                   of RESULT.json, MODEL.json and eval_selfreport.json
  validation/      MC pseudo-experiments
```

- Every stage folder has `src/` (C++), `scripts/` (Python), `config/`,
  `condor/` and a `Makefile`.
- Build products go to `analysis_v3/<stage>/bin/`, git-ignored.

### 2.2 Products

```
production_v3/
  pseudo_data/917f9326efa44ea9/  input, read-only
  work_reading/                  reading-phase inspections
  tmp/                           TMPDIR and tool caches of local runs
  tmp/condor_submit/             eossubmit submission-time TMP
  external/                      downloaded third-party code (JHUGen-MELA)
  program_inputs/<stage>/<sha>/  content-frozen executables, configurations
                                 and payloads for workers
  manifests/<ver>/               catalogue, MC and data manifests, integrity
                                 reports
  skims/<ver>/<sample>/          dilepton, multilepton and generator skims
                                 with sidecars
  calibration/<ver>/             Z-peak fits, scale/smear and D_mass
                                 payloads, plots
  tnp/<ver>/                     efficiency histograms and fits, SF
                                 payloads, trigger study, plots
  optimization/<ver>/            selection and electron-ID optimization
                                 (MC only)
  h4l_reco/<ver>/<sample>/       4l candidate trees (signal and control
                                 regions), plots
  backgrounds/<ver>/             qqZZ/ggZZ templates, fake rates, OS and SS
                                 predictions, combination, plots
  signal_model/<ver>/            DCB fits against m_H, morphed templates,
                                 YR4 inputs, plots
  fits/<ver>/                    workspaces, fits, scans, impacts, GoF toys,
                                 plots
  results/<ver>/                 RESULT.json, MODEL.json,
                                 eval_selfreport.json, results_full.json,
                                 summary.md, key-number text files
  validation/<ver>/              pseudo-experiment outputs and summaries
  condor/<stage>/<UTC stamp>/    JDL, queue, task JSONs, stdout, stderr,
                                 UserLog
```

Everything a Condor job or eossubmit reaches uses
`/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/...` or
`root://eosuser.cern.ch//eos/user/...`, never `/eos/home-y/...`.

## 3. Stages

### Stage 0. Preservation and rules

- Copy `AGENTS.md` and `README.md` to `deliverables/pre_v3_rules_20260923/`.
- Update `AGENTS.md` per prompt section 12, adding the user rules logged on
  2026-09-23.
- Rewrite `README.md` for the v3 workflow; the old notes stay in the
  preserved copy.
- Update `.analysis_config` to the v3 paths.
- Apply the `.gitignore` changes and run the check.

### Stage 1. Infrastructure (`framework/`, `common/`)

**Components**
- **Sample catalogue** `samples_ul16_v3.json`. For each of the nine samples:
  - role, sigma_eff, directory and blocks;
  - weight type, and the coverage of partial records.

  It also holds the data: the path, and the luminosity taken as a parameter
  (20 fb^-1, checked against the provenance `lumi_fb`).
- **C++ input validator**, per file:
  - open through XRootD;
  - check the keys (`Events`, `Runs`, provenance);
  - check that every branch any stage needs is present with its type;
  - read the first and the last entry;
  - sum the `Runs` bookkeeping (count, sum w, sum w^2) and read
    `entries_before_selection`;
  - MC only: sum genWeight and genWeight^2 over `Events`;
  - data only: compare the schema, `lumi_fb` and branch list across shards,
    and record entries only for coverage;
  - write a JSON sidecar.
- **Manifest builder** (Python). It makes the only directory listing and runs
  the validator.
  - An unreadable data shard is a hard failure.
  - An unreadable MC file, after a bounded retry, is recorded together with
    its coverage loss.
- **Condor framework:**
  - stager: splits tasks by bytes, builds the JDL from a template with the
    eossubmit settings, dry-run-checks the arguments and all three log paths,
    records the cluster ID;
  - worker wrapper: minimal `PATH`, `./NAME` normalization, output
    validation, atomic publication, attempt-qualified preservation,
    idempotent `mkdir -p`;
  - output scanner and monitor.
- **Common C++ library:**
  - event reader with optional event identifiers;
  - smearing seeds = hash(manifest file id, entry, object index, variation);
  - atomic ROOT publisher;
  - payload readers: half-open bins, unbounded tail bin.

**Normalization.** genWeight-based: L sigma_eff (sum of genWeight over the
selected events) / (sum of `genEventSumw` over the files present).  The
count-based value is a documented cross-check.

**Condor.** About 65 jobs, 1 CPU, 2 GB, flavour microcentury.

**Validation**
- smoke test on one data shard and one file per sample;
- eossubmit dry run;
- manifest completeness report;
- count- vs genWeight-based comparison on all files.

### Stage 2. Object definitions and skims (`skims/`)

**Analysis configuration** `analysis_ul16_v3.json`:
- triggers: AN Table 2 in NanoAOD names, plus the single-lepton paths for
  tags;
- the object definitions, and the FSR (muons only, with the AN criteria);
- jets: pT > 30, abs(eta) < 4.7, `jetId >= 2`, dR > 0.4 cleaning, DeepJet
  medium 0.2489;
- calibration and tag-and-probe binning.

**Skim floors.** Every skim keeps raw values above the NanoAOD collection
floors (muon pT > 3, electron pT > 5 GeV, loose impact parameters), so that
the optimization of stage 4 can move any threshold without re-skimming.

**One C++ pass** over every data shard and MC file writes three outputs:
1. **Dilepton table.**
   - Every OS same-flavour pair with 40 < m_ll < 140 GeV (mass without FSR;
     the FSR-corrected mass is stored too).
   - For both legs: kinematics, momentum errors, ID flags and raw MVA
     scores, isolation, SIP, dxy, dz, trigger-object matches.
   - Per event: trigger bits, MET, PV, weights.
   - Also tag-electron plus `Photon` pairs, and tag-muon plus HLT
     muon-object pairs, for the reconstruction efficiency.
2. **Multilepton slim NanoAOD.**
   - Events with at least 3 leptons above the floors.
   - The branches the later stages need: leptons, FSR photons, jets, MET, PV,
     analysis trigger bits and `TrigObj`.
   - For MC also the generator, LHE/PS-weight and HTXS branches.
3. **Signal generator table** (signal MC, every event):
   - generator 4l, dressed leptons (photons within dR < 0.3), isolation
     (dressing photons excluded), fiducial flags;
   - HTXS stage 0 and y_H, generator jets, VH decay class.

The volume is measured in the smoke test (target: at most about 150 GB).

**Condor.** About 470 jobs, 1 CPU, 2-3 GB, flavour longlunch.

**Validation**
- data/MC comparisons: lepton kinematics, ID variables, m_ll, N_PV and rho
  (a pileup-sensitive check: any difference is reported), object counts,
  trigger rates;
- full output scan.

### Stage 3. Lepton calibration and tag-and-probe (local, on the dilepton tables)

**3a. Momentum scale and resolution** (`calibration/`; implemented as calibration v3,
user decisions 2026-09-24: MC template fits, factorized model)
- **Control pairs:** analysis trigger OR; both legs AN tight (PF muon with
  global-or-tracker-with-station; electron mvaFall17V2noIso WP90 until
  stage 4b; FSR-subtracted iso < 0.35, SIP < 4), thresholds on the
  transformed pT.  Condor extraction `calib_extract` (v2: 43 tasks, scanned).
- **Model per flavour:** ln(1 + s) = a(|eta|) + b_R(pT), r^2 = c(|eta|) +
  d_R(pT); a, c piecewise constant on fine |eta| bins (muons 11 bins, 1.2 at
  a bin centre; electrons 13 bins in |eta_SC| keeping 1.4442/1.566); b, d per
  coarse region (muons 3, electrons 2), linear in pT between the mean-pT
  nodes, zero at the 40-50 GeV reference bin.
- **Categories:** family A = pairs of fine |eta| bins, both legs pT > 20 GeV
  (no pT-bin sculpting); family B = pairs of (pT, region) bins.  Only
  categories whose data mode lies in 80-100 GeV enter.
- **Template fit per category:** event-level MC templates (each MC pair moved
  to k m (1 + sqrt(delta^2 + D) eps) with frozen pair deviates, fixed 0.3 %
  Gaussian kernel), data smeared once by the common delta = 0.5 %, likelihood
  groups of >= 20 MC entries, free normalization, Barlow-Beeston lite.
- **Joint least squares** of both families with the measured leg
  compositions and the per-category response (measured at iteration 1 by
  injecting -1 % scale and 1.2 % smear into the MC itself; accepted range
  0.5-1.5), outlier rejection (|pull| > 5, <= 10 %), a pT smoothness prior
  (1e-3 on ln(1+s), 3e-5 on r^2), category set frozen from iteration 2;
  errors inflated by chi2/ndf.
- **Iteration and finalization:** data pT exp(-u), MC pT (1 + r N); stop at
  0.2 sigma or, after 10 iterations, the mean of the last 4 applied payloads
  with their rms added to the uncertainty (finalize_calibration.py).
- **Nominal result:** production_v3/calibration/v3/nominal_r3 (payload.json,
  plots_r2/, report.json).
- **Still to do:** closure runs (sloped_same, null_halves, uniform_halves,
  large_halves) -> method non-closure systematic; variants (window, kernel,
  delta, prior x2, family-A pT threshold) -> fit-model systematic;
  charge/eta/phi diagnostics of the final payload; a second pass with the
  tag-and-probe SFs applied to the MC templates.

**3b. Per-event mass uncertainty** (AN 5.3.1; implemented as lambda v2:
`calibration/src/lambda_histograms.cpp`, `fit_lambda.cpp`, `scripts/run_lambda.py`,
config `lambda_ul16_v2.json`; nominal payload
`production_v3/calibration/v3/nominal_r3/lambda_v2_r2/lambda.json`)
- Input: FullPairs of the final calibration payload (data corrected, MC
  smeared), written by diagnose_calibration.py --final.
- Regions as the AN: muons |eta| 0-0.9-1.8-2.4; electrons |eta_SC| with
  d = dpT/pT < 0.03 in 0-0.8, 0.8-1.0 and d < 0.07 in 1.0-1.2-1.44-1.57-2.0-
  2.5, plus the high-error classes |eta| < 1, d > 0.03 and |eta| > 1,
  d > 0.07 (AN reference-electron method: one leg in |eta| < 0.8, d < 0.03
  with its lambda fixed).  No ecalDriven flag in UL NanoAOD (no separate
  tracker-driven corrections).
- The binned form of the AN conditional fit: all bins of the predicted
  relative mass error of one pair class fitted simultaneously with
  BW (x) DCB + exponential, shared tails, sigma_k = (m_Z/2) sqrt(lambda_a^2
  <d_a^2> + lambda_b^2 <d_b^2>); MC first, data with the MC tails; closure
  per bin (AN Figure 18) and fit galleries.  A core-window variant (82-100
  GeV) gives the lambda method uncertainty (1-3 %).
- Usage: every lepton error d is multiplied by lambda of its region (data
  and MC each their own); D_mass and the Z1 refit use it.

**3c. Tag-and-probe** (`tnp/`; config `tnp_ul16_v5.json`, run `production_v3/tnp/v2/v5`)
- **Extraction** `tnp_extract` (v2: 43 tasks, scanned): pairs with raw mass
  50-130 GeV and a tag-candidate leg, no probe requirement.
- **Tag** (calibrated pT, path-consistent): IsoMu24 with the Iso bit or
  IsoTkMu24 with the IsoTkMu bit, trigger pT >= 24, pT > 26, |eta| < 2.4, POG
  tightId, pfRelIso04 < 0.15; Ele27_WPTight (trigger pT >= 27) or
  Ele25_eta2p1_WPTight (>= 25, |eta_SC| < 2.1), pT > 30, |eta_SC| < 2.1
  outside the gap, cut-based tight.  (Substitution: the unprescaled 2016 tag
  paths instead of the AN IsoMu20/22.)
- **Probe:** AN loose lepton (muons pT > 5, |eta| < 2.4, global or tracker
  with a station; electrons pT > 7, |eta| < 2.5), calibrated; dR(tag, probe)
  > 0.2; per tag and step only the probe closest to m_Z.
- **Steps (the selection chain, user rule):** id | loose (AN tight muon;
  electrons mvaFall17V2noIso WPL, the stage-4b starting point); sip (< 4) |
  id; iso (AN FSR-subtracted < 0.35) | id, sip; full | loose measured
  directly as a cross-check of the product.  Fine bins after
  efficiency_tnp_an: muons pT 5-10-15-20-22.5-...-45-50-60-80-inf x |eta|
  0-0.2-...-1.6-1.85-2.1-2.4; electrons pT 7-12-17-20-22.5-...-80-inf x
  |eta_SC| 0-0.3-0.6-0.9-1.2-1.479-1.7-1.9-2.1-2.3-2.5; report point (45, 1.2)
  by bilinear interpolation between bin centres.
- **Fits** (binning and quality control after efficiency_tnp_an; line shape
  by the user decision of 2026-09-24): pass and fail simultaneously with
  the efficiency as a parameter, window 60-120 GeV, the narrowest of
  0.5-3 GeV bins with >= 40 entries per bin; nominal signal = the template
  of the prompt-prompt DY-MC pairs of the bin and category (fallback: pT
  row, then inclusive, below 500 effective entries; smoothed 0.5 GeV)
  convolved with a Gaussian per category (shared below 300 entries), which
  carries the genuine fail signal with degraded mass (FSR photons inside the
  isolation cone: the 60-80 GeV hump of the fits_i1 inspection);
  background CMSShape (erf x exp; a category below 100 entries keeps the
  inclusive shape); alternatives: analytic stand-alone DSCB (the fail one
  with a Gaussian low-mass component) and a Bernstein-3 background (fit
  systematic = largest deviation among the alternatives that the inspection
  accepts); efficiency errors from MINOS (sandwich-scaled for weighted MC)
  with the binomial floor; rejection of narrow background maxima near m_Z;
  categories below 20 entries counting only; counting cross-check; the MC
  generator truth is a diagnostic only (user: efficiencies from
  tag-and-probe only).  Rounds: v4 fits_i1 (analytic model, inspected;
  superseded), v5 fits_i1 (template model).
- **Inspection (user rule):** every fit of every round is looked at on the
  gallery pages (8 fits per page, pass | fail, linear scale, pulls); the
  classes of bad fits are fixed by model changes or recorded overrides in a
  new fit label (fits_i1, fits_i2, ...) until no bad fit remains; never
  judged by chi2 or pulls alone.
- **Reconstruction efficiency:** muons: Z yield ratio, TrigObj probes and
  MC; electrons: Photon probes; best agreeing method, differences as a
  systematic.  **Trigger:** AN 4l method with 4l events plus per-leg T&P.
- **Outputs:** SF payload per step and bin with stat and fit-model errors;
  per-event weight; eff_correction; muEffScale/elEffScale at (45, 1.2).
  The final run uses the stage-4b working points.

### Stage 4. 4l reconstruction and optimization (`reconstruction/`)

**4a. Reconstruction** (implemented in two programs).
- **h4l_reco** (Condor, config `reco_ul16_v1.json`): per event with >= 3 AN
  loose leptons and an OS SF pair of 30-120 GeV, every loose lepton
  (calibrated pT, raw pT and smearing deviate, error, FSR-subtracted
  isolation, IDs, AN FSR photon, HLT match, generator match), the jets
  (tight ID, CSVv2 and DeepJet) and every ZZ candidate under loose floors
  (both orderings of OS-OS pairs, SS Z2, smart-cut masses); theory weights
  for MC; generator-table fields for signal.
- **h4l_select** (Condor, config `selection_ul16_v1.json`, MELA): the AN
  selection: trigger OR, good PV, selected leptons (tight ID, SIP < 4,
  iso < 0.35), e-mu cross cleaning, SR and SRZ4l (Z2 > 12 / 4 GeV), best
  candidate by D_bkg^kin (same four leptons: Z1 closest to m_Z); the
  control regions 2P2F, 3P1F (best by D_bkg^kin), SS, and Z + 1 loose
  lepton; per candidate the lambda-corrected per-event mass error, the Z1
  refit (true Z1 line shape of the ggH generator record, adaptive KDE,
  refit inputs v3), MELA probabilities and discriminants (MELA
  normalization; D_WH/D_ZH constants from MC downstream), jets cleaned from
  the selected leptons and FSR photons, DeepJet medium b tags, additional
  leptons, MET, the MC m4l responses to the lepton scale and smearing.
- **Downstream** (templates): normalization, the SF weights (chained
  efficiency SFs), the categories (paper working points 0.5), the
  variations.
- **Validation:** cut flows; m4l in 70-170 and 70-800 GeV per final state
  and category; m_Z1, m_Z2 and m_Z1 vs m_Z2 in 118-130 GeV; the Z->4l peak;
  paper Tables 1 and 2; the refit and D_mass closure on signal MC.

**4b. Optimization** (MC only; `reconstruction/scripts/optimize_selection.py`,
first scan `production_v3/optimization/v3/scan_i1/`; decisions: electron
noIso WPL, SIP < 4, isolation < 0.35, SR Z2 > 10 GeV, Z1 > 40 GeV, pT 20/10;
the categories and discriminant working points still to do).
- **Scanned:**
  - lepton pT thresholds, isolation, SIP, m_Z1 and m_Z2 minima, the OS-pair
    minimum;
  - the muon ID WP;
  - the electron ID: cuts on the raw Fall17V2 BDT per (pT, eta) region as in
    AN Table 7, or a BDT/NN trained on MC (TMVA, or a pixi package), prompt
    against non-prompt electrons, with separate training and evaluation
    samples.
- **Objective:** minimize (sigma_mu/sigma_mu,AN)^2 + (sigma_mH/sigma_mH,AN)^2
  from Asimov fits of a simplified model (mu = 1, m_H = 125 GeV,
  20 fb^-1).  Z+X comes from the fake-rate method applied to MC.
- **Outcome:** the results are shown to the user; if several choices remain
  close, the user decides.  The chosen values are frozen in the analysis
  configuration, and tag-and-probe is re-run for them.

### Stage 5. Backgrounds (`backgrounds/`)

- **qqZZ, ggZZ.**
  - MC templates with the calibration and the SFs; sigma_eff with no
    k-factors.
  - Set (a) adds the qqZZ QCD scale, PDF and PS, the AN k-factor/EW
    uncertainty, and 10 % on ggZZ.
  - Set (b) keeps the MC statistics and the detector nuisances.
- **Z+X, OS method** (AN 7.2.2): fake rates from Z1 + 1 loose lepton, with
  the lessons of AGENTS.md; 2P2F and 3P1F with the prompt ZZ subtracted; a
  non-negative template that preserves every signed integral.
- **Z+X, SS method** (AN 7.2.1): conversion correction; r_OS/SS from DY + TT
  MC.
- **Combination:** inverse variance per final state; Landau shapes; closure
  on MC.

### Stage 6. Signal model (`signal_model/`)

- **Samples.** Per production mode (ggH, VBF, and VH split by HTXS stage 0
  into VH-had, WH-lep, ZH-lep), final state and category: DCB fits at
  125 GeV.
- **m_H dependence** by rest-frame scaling of the 4l system (user decision):
  - scale by k = m_H/125 in the generator H rest frame, keeping the H pT and
    rapidity;
  - DCB parameters on a grid, parameterized continuously in m_H.
- **Systematics of the morphing:**
  - the alternative m4l-shift morphing;
  - the variation of the A x eps slope;
  - the effect of the unphysical Z1-mass scaling, estimated by reweighting
    the generator (m_Z1, m_Z2) distribution.
- **Normalization:** sigma_eff(mode, 125) x YR4 ratios x A x eps(m_H).
- **Non-resonant VH:** a Landau shape.

### Stage 7. Statistical model (`inference/`)

- **Likelihood.** Extended, unbinned, per final state x category:
  P(m4l) x P(D_bkg^kin | m4l); the mass fits in 1D, 2D and 3D, each with and
  without the Z1 refit.
- **Parameter ranges** (user decision):
  - mu unbounded;
  - m_H free in a wide range (initial 110-140 GeV, with the m4l window
    widened accordingly, initial 100-160 GeV);
  - final values from toys, with no boundary artifacts.
- **Discovery:** a local significance of at least 5 sigma.
- **Nuisances.**
  - Set (b): lepton scale and resolution per flavour, efficiency, trigger,
    Z+X, MC statistics, A x eps(m_H), shapes, D_mass, morphing.
  - Set (a) adds: luminosity 2.5 %, BR 2 %, qqZZ theory, ggZZ 10 %, signal
    QCD/PDF (YR4), PS, JES/b tag.
- **Intervals.** From the profile likelihood; the asymptotics are checked with
  toys; Feldman-Cousins where needed.
- **Stat/syst.** Fixed-nuisance impacts.
- **Goodness of fit.** Saturated test with toys; pulls and impacts.
- **MODEL.json:** binned pyhf templates at the fitted m_H, with the set-(b)
  nuisances and a consistency check.

### Stage 8. Results

- Every item of prompt sections 6.1 and 6.2, in both result sets.
- The fiducial fit with m_H profiled, plus a variant at 125.09 GeV.
- The files `RESULT.json`, `MODEL.json`, `eval_selfreport.json`,
  `results_full.json` and `summary.md`, plus judge-readable text tables.
- The paper-comparison figure.

### Stage 9. Validation (`validation/`)

- **Detector chain:** two independent halves of the MC; the half with
  injected generic effects plays the data.
- **Fit level:** toys at m_H = 121, 123, 125, 127 and 129 GeV and
  mu = 0, 0.5, 1, 2 and 3, with at least 500 toys per point.
- **Tests:** every test of prompt section 8.

### Stage 10. Freeze, unblind, AN, handoff

- **Freeze.** Freeze the code, the configuration and the results with
  checksums.
- **Unblind.** Ask the user where the truth comes from (user decision).
- **AN.** Written in LaTeX and compiled with pdflatex, in
  `deliverables/<versioned name>/`.
- **Handoff.** A final independent review, then the handoff report.

## 4. Resources and volume (estimates, re-measured in every smoke test)

| stage | where | jobs | per job | new volume |
|---|---|---|---|---|
| 1 manifests | Condor | about 65 | 1 CPU, 2 GB, under 1 h | under 50 MB |
| 2 skims | Condor | about 470 | 1 CPU, 2-3 GB, 1-3 h | at most about 160 GB |
| 3 calibration, tag-and-probe | local (Condor if slow) | - | 8 processes | under 5 GB |
| 4 reconstruction, optimization | local or Condor about 50 | - | - | 2-10 GB |
| 5-6 backgrounds, signal model | local | - | - | under 3 GB |
| 7-8 fits, results, toys | local + Condor about 50 | - | 1 CPU, 2 GB | under 1 GB |
| 9 validation | Condor about 200 | - | 1 CPU, 2 GB, 1-4 h | under 5 GB |

## 5. Failure modes of earlier reference pipelines and the design answer

| failure mode | design answer |
|---|---|
| grid-edge attractor | continuous m_H parameterization; MLE by Minuit; wide m_H range; toys near the edges |
| truncated-Gaussian width correction | none anywhere; mu unbounded |
| FSR-driven left shoulder | DCB tails from MC per final state; low-statistics toys |
| loose background normalization plus free slopes | normalizations fixed within their uncertainties; only constrained shape nuisances; paired-injection test |
| bounded mu plus mass profiling | mu unbounded; null datasets studied with toys |
| fixed Z-peak window | data-seeded windows, tested with injected +-5 % shifts |
| bootstrap/linearity without the measured shift | none used; toys with the measured calibration |
| signal left at the nominal efficiency | SFs applied to all MC; thinned-MC validation |
| smear in the dilepton convention | per-lepton convention by construction; injected-smear closure |

## 6. Decisions (2026-09-23; details in `experiment_log.md`)

| item | decision | source |
|---|---|---|
| code tree | `analysis_v3/` | user rule 21:32Z, layout mine |
| products, Condor artifacts, TMP | section 2.2 | user rule 21:19Z, layout mine |
| AGENTS/README copies | `deliverables/pre_v3_rules_20260923/` | mine (21:48Z permission) |
| selection thresholds | AN start, optimized on MC | user, round 1 |
| electron ID | study the thresholds; train a BDT/NN if useful | user, round 1 |
| muon ID | AN PF (`looseId`) baseline, in the optimization | mine |
| muon reconstruction efficiency | yield ratio + TrigObj T&P + MC; best agreeing method, differences as a systematic | user, round 2 |
| morphing | rest-frame scaling, with a systematic | user, round 1 |
| mu, m_H range, discovery | mu unbounded; wider m_H range; local 5 sigma | user, round 2 |
| MC normalization | genWeight-based | user, round 2 |
| trigger efficiency | 4l method + T&P, Z-peak per-leg if data are short | user, round 2 |
| electron FSR | muons only | mine (no alternative) |
| discriminants | JHUGen-MELA | mine (AN method; install allowed) |
| jets and b tag | DeepJet medium, AN categories | mine (only usable tagger) |
| width | BW (x) resolution, no interference | mine (no alternative) |
| fiducial isolation | dressing photons excluded | user, round 3 |
| fiducial m_H | profiled, plus a 125.09 variant | user, round 3 |
| MODEL.json | binned pyhf at the fitted m_H | user, round 3 |
| result sets | prompt proposal | user, round 3 |
| optimization objective | joint mu and m_H precision | user, round 4 |
| unblinding | ask where the truth comes from | user, round 4 |
| calibration block | per-flavour values separately; RESULT.json required single-value calibration fields null, per-flavour values as extra keys; ask whenever a required format conflicts with our results | user, 2026-09-25 (supersedes round 4) |
| efficiency SF reported | full single-lepton selection | user, round 4 |
| fits | compiled C++ RooFit/Minuit2 | mine |
| AN build | LaTeX + pdflatex, "Private work" | mine |
| lepton calibration fit | MC template fits per leg-bin category, iterated per lepton; BW (x) DCB for validation and lambda | user, 2026-09-24 |
| fit quality | judged on the plotted fits, not chi2/pulls alone | user, 2026-09-24 |
| fit functions and limits | after `efficiency_tnp_an` where applicable | user, 2026-09-24 |
| convolution | own numerical BW (x) DCB (no FFTW plugin in the lxplus ROOT) | mine (no alternative) |
| efficiency chains | measured directly or as conditional chains (reco -> id -> iso -> sip), never independent products | user, 2026-09-24 |
| discriminant constants | paper convention D' = 1/(1 + c P_bkg/P_sig), constants and working points optimized on MC (categories_ul16_v2: c_2jet 0.1, c_1jet 0.0014, c_WH 0.16, c_ZH 0.4, WP 0.5, VH-MET MET > 90 GeV) | user, 2026-09-25 (supersedes MELA normalization of 2026-09-24) |
| b tagging | DeepJet medium 0.2489 (UL16 post-VFP) | user, 2026-09-24 |
| lambda method | binned simultaneous AN conditional fit, AN regions and reference method | mine (AN method) |
| Z1 refit line shape | ggH generator Z1 (dressed), adaptive KDE | mine (AN method) |
| selection thresholds (stage 4b) | electrons WPL, SIP < 4, iso < 0.35 kept, SR Z2 > 10 GeV | user, 2026-09-24/25 |
| T&P fit model | DY-MC template (x) Gaussian + CMSShape (EGM ranges); the fail template mixed with a free pass-like fraction phi (the bin's pass template); alternatives DSCB (+ low-mass Gaussians in pass and fail) and Bernstein-3 | user, 2026-09-24 and 2026-09-25 (admixture) |
| T&P template smoothing | adaptive (Abramson) kernel, Silverman pilot of the Z-peak entries, nearest-neighbour floor of 100 effective entries; fallback |eta| column (pT within a factor 2, then 4) -> pT row -> inclusive; min_effective_entries 1000 | mine (v7, after the v5/v6 inspections) |
| weighted-MC fit errors | effective-weight scaling sqrt(sum w^2 / sum w) of HESSE and MINOS (T&P) | mine (the numerical sandwich failed in weakly constrained fits) |
| lepton scale / resolution nuisances | from the calibration, its halves closures and the residual dilepton widths (make_systematics.py): scale mu 3.7e-4, e 5.0e-4; resolution mu 2.6 %, e 3.4 % | mine (derived; set b) |
| Z+X systematic | combined-final-state MC closure of both methods, max(|1 - r|, sigma_r) of the larger, relative, all final states, with the bootstrap stat; AN composition estimate as a cross-check | user, 2026-09-25 |
| 3D per-event resolution | signal m4l density conditional on D_mass = dm/m4l, DCB width s x D_mass x m_H (s fitted on MC with D_mass x 125 GeV); the earlier s x dm (dm growing with m4l at fixed D_mass) left the signal pdf unnormalized by 0.2-0.7 % | mine (root-cause fix 2026-09-25, CMS L(m4l given m_H, D_mass)) |
| eval_selfreport resSmear | value and unc null, per-flavour extra keys resSmear_muon, resSmear_electron | user, 2026-09-25 |
| sub-POI ranges | per-category, per-mode, (mu_F, mu_V), STXS and fiducial POIs in [0, 20] (paper); inclusive mu unbounded | user, 2026-09-25 |
| N-1 significances | cuts relaxed to the event-record floor (relaxed_to stated), data-driven Z+X per selection, ID variants, one program with a checked reference | user, 2026-09-25 |
| final T&P | run v7 fits_i1c; later iterations only for bins not yet fine | user, 2026-09-25 |
| final reviews | two reviewer subagents | user, 2026-09-25 |

## 7. Status (updated 2026-09-25T21:05Z)

| item | status |
|---|---|
| Stages 0-6 | done (T&P final: v7 fits_i1c; Z+X zx_v2; signal sm_v2; N-1 nm1_v1) |
| Stage 7 | done: models model_*_r1 (efficiency-uncertainty split, renormalized Z+X fractions, Z+X -1 sigma fixed, 2Dmass dimension, Asimov D_mass sub-points, systematics v3) |
| Stage 8 | results r2 (production_v3/results/v5/r2); the paper-style fits rerun without the refit in set (a); final assembly as a new output tag |
| Stage 9 | done: toys (mu 0-3, 500 per point; m_H 121-129, 400; paired), MC-event closure, D_mass validation, T&P and calibration closures |
| Stage 10 | final reviews: round 1 (2 reviewers) fixed, round 2 no blocker; freeze, the truth question, the AN (make_an.py) |

## 8. Execution order of the remaining work (2026-09-24, user: complete everything without stopping)

1. MELA (JHUGen-MELA in production_v3/external/, built against the system
   ROOT; a C++ interface for D_bkg^kin, D_2jet, D_1jet, D_WH, D_ZH) - in
   parallel with 2-4.
2. Tag-and-probe: tnp_histograms (data and MC, calibrated), fit_tnp
   (simultaneous pass/fail), SF payload, galleries; reconstruction and
   trigger efficiencies.
3. Lambda calibration on the final FullPairs.
4. Calibration closures and variants (systematics), two at a time.
5. 4l reconstruction (C++, Condor on the skim Events): objects, FSR, iso,
   SIP, cleaning, Z1/Z2, smart cut, m4l > 70, candidate choice by D_bkg^kin
   (MELA; documented substitute if unavailable), jets and categories,
   D_mass, Z1 refit, control regions (Z1+l, 2P2F, 3P1F, SS), MC weights with
   the calibration, SFs and variations, generator information; cut flows and
   validation plots; stage-4b optimization on MC.
6. Backgrounds: qqZZ/ggZZ templates; Z+X OS and SS, combination.
7. Signal model: DCB per final state and category vs m_H (rest-frame
   morphing), YR4 sigma x BR(m_H), A x eps(m_H).
8. Statistical model and every section-6 result (mu inclusive/final
   state/category/mode, 2D contours, STXS stage 0, fiducial and
   differential, m_H 1D/2D/3D with and without the refit, width, Z->4l,
   yield tables, systematics and impacts, GoF).
9. Benchmark files (RESULT.json, MODEL.json, eval_selfreport.json, full
   results JSON, summary) and validation toys.
10. The detailed AN (pdflatex), final review, handoff.
