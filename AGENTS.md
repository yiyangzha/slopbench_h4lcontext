# H->4l UL16 PFNano v3 reference analysis: persistent working rules

Keep this file current.  Record every new durable user requirement here
immediately, and replace superseded rules instead of keeping contradictory
history.  The 2017 rules are preserved in
`deliverables/pre_v3_rules_20260923/AGENTS.md`.

Before every new stage, Condor submission and report, re-read this file,
`PLAN.md` and the recent part of `experiment_log.md`.

## Documents

- `AGENTS.md`: the durable rules (this file).
- `README.md`: the concise, ordered, runnable workflow.
- `PLAN.md`: the living plan and the table of decisions.
- `experiment_log.md`: the append-only, dated record of plans, decisions,
  questions and answers, submissions, checks and results.

Update all of them in the same work item as any change to code,
configuration, workflow or outputs.

## User working rules (highest priority)

- **Asking.** Ask the user directly (interactive questions) whenever
  something is unclear, ambiguous or uncertain.  Never fill a gap with an
  assumption.  Whenever several methods are possible, ask which one to use.
  Questions to the user are written in Chinese (user, 2026-09-24).
- **Correctness over precedent (user, 2026-09-24).** Use what is correct.
  Earlier requirements, rules and code are not necessarily right: fix what
  is certainly wrong (and record it), ask the user whenever unsure, and
  never keep a wrong definition only for consistency with earlier work.
  The goal is a correct analysis that earns full marks on every score of
  the benchmark (`ref/Slopbench H4l v0 (Copy).html`, prompt section 9).
- **Permissions.** These may be done directly, without asking:
  - install what is needed: add missing pixi packages to `pixi.toml` and run
    `pixi install`; put other downloaded tools in `production_v3/external/`;
  - change code and configuration;
  - submit Condor jobs.
- **Limits.**
  - Change only files inside this repository directory, except the
    authorized output locations below.
  - **Read access (user, 2026-09-25).** Outside this repository, open,
    list or search only what the user has authorized: the data and MC
    inputs of "Physics scope", the truth record the user named for the
    comparison (and the program that generated it, only to interpret the
    format and only after asking the user for its location), and the paths
    named explicitly elsewhere in this file.  Never browse or search other
    directories of the user (other repositories, the group area, other seed
    folders); ask the user first whenever something outside is needed.
  - Never delete files or directories.  If a deletion seems really
    necessary, ask.
  - Keep failed or superseded artifacts under attempt-qualified names, and
    write new versioned outputs instead of overwriting.
- **Temporary files.** Temporary files and tool caches go inside the
  repository (`production_v3/tmp/`, e.g. `TMPDIR`, `PIXI_CACHE_DIR`,
  `APPTAINER_CACHEDIR`), never into shared `/tmp` or `$HOME`.
- **Old and new code.** Keep them in separate folders.
  - The new code lives only in `analysis_v3/`.
  - The 2017 folders (`production_framework/`, `calibration_skim/`,
    `h4l_reconstruction/`, `reducible_background/`, `inference/`,
    `fake_data_nanoaod/`, `ntuplizer/`) are reference only and are never
    modified.
- **Correctness over reuse.** The existing code is only a reference.  What
  matters is correctness, completeness and consistency with the paper and
  the AN.
- **Fidelity.** Follow the paper as closely as possible.  Where the inputs
  make an identical implementation impossible, a reasonable and accurate
  substitute is enough; document every substitution.
- **Division of labour.** Read large ROOT files with C++ (faster).  Use
  Python for plotting, orchestration and similar work.
- **Condor.** Submit only through eossubmit, with correctly set parameters.
  Every output, scientific or Condor, goes into a folder inside
  `production_v3/`.
- **Records.** Record each user requirement promptly (here or in the log).
  Keep the log and the plan up to date, and review them regularly.
- **Unattended mode.** It applies only when the user explicitly declares it.
- **Complete the whole analysis without stopping** (user, 2026-09-24): the
  task is the full prompt_renew.md scope, every stage through the detailed
  AN and the benchmark files, not a single stage.  Do not pause after partial
  results; report progress briefly and keep working.  A genuine method
  choice is asked while other work continues.
- **Records and review** (user, 2026-09-24; reaffirmed 2026-09-25): record
  every requirement (here or in the log), the log and the plan promptly, in
  the same work item as the decision; re-read AGENTS.md, PLAN.md and the
  recent log before every stage, at every stage boundary and regularly while
  long work runs, so that nothing is forgotten or done wrongly.
- **Git.**
  - Never commit, push, create or switch branches, stash, reset or clean.
  - Never run `git add -A` or `git add .`.
  - Never restore the deleted tracked files (`CLAUDE.md`, `phase*/CLAUDE.md`,
    `phase5_documentation/outputs/references.bib`) without the user's
    confirmation.
  - Ignore generated products with exact anchored patterns.  After any
    `.gitignore` edit, `git ls-files | git check-ignore --stdin --no-index
    --verbose` must print nothing, and no source directory may become
    ignored.
- **Never write to:**
  - any slopbench repository;
  - `h4l_seeds*` or `data_premix/`;
  - the delivered MC or the pseudo-data directory;
  - the old group area `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/`;
  - the old Condor area `/eos/user/y/yiyangz/h4l_2017_nanoaod_analysis/condor/`.

## Order of authority

1. The user; the newest explicit instruction wins.
2. This file.
3. `slopbench_code_fork/AGENTS.md`.
4. The jfc framework (`agents/`, `conventions/`, `methodology/`).

On the physics method, AN-16-442 v8 and JHEP 11 (2017) 047 are the
authority.  Never edit the framework symlinks and never re-run the
scaffolder.

## Blinding (until the user unblinds)

- **Never list, open or read:**
  - `slopbench_code_fork/injection/records/`;
  - `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/ref_v3/injection_tmp/`;
  - anything under `slopbench_code_fork/results/h4l_seeds/` other than
    `v3/an/` and `v3/inputs/Higgs_XSBR_YR4_update.xlsx`;
  - any `eval/truth/` directory.
- **User-authorized exception** (2026-09-24): the fit functions, parameters
  and limits of `slopbench_code_fork/results/h4l_seeds/efficiency_tnp_an/`
  (its note's method sections, its fit configuration
  `slopbench_code/benchmark/configs/h4l_tnp.json`, the DY-MC fit records)
  may be consulted.  Its seed results, paired closures and predicted
  injection ratios are never read or used.
- **Never infer the dataset's configuration** from its entry count, shard
  sizes, the injection catalogue or other metadata.  Read
  `PFnanoFakeDataProvenance` only for schema and luminosity checks.
- **Injection definitions** (functional forms, profiles) may be read to
  understand the benchmark.  They are never analysis inputs, fit models,
  priors or snapping targets.
- **Independence.** Every entry of the data and of the MC is an independent
  event; the data and the MC are unrelated.
  - Use ordinary Poisson statistics for the data and ordinary MC statistics
    for the MC.
  - Apply no repetition correction (no `k_eff`).
  - Do no event fingerprinting or matching between entries or samples.

## Physics scope

- **Inputs.**
  - Data: `/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/pseudo_data/917f9326efa44ea9/part_00000.root`
    ... `part_00082.root` (83 shards, 20 fb^-1, UL16), read-only.  Treat
    them exactly as collision data.  The data path and the luminosity are
    parameters.
  - MC: `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l_seeds_v2/mc/`
    only, with the effective cross sections of `prompt_renew.md` 4.2
    (TTBar 87.58 pb).
  - Never use modified or intermediate samples.
- **Normalization.** genWeight-based (user decision): L sigma_eff (sum of
  genWeight over the selected events) / (sum of `genEventSumw` over the
  files present).  The count-based value is a documented cross-check.  No
  k-factors or EW corrections in the central prediction.
  - A consumer re-sums the per-file `genEventSumw` over exactly the files it
    processed, and records that file set.
  - Theory variations are normalized with the Runs `LHE*Sumw` sums of the
    same files, because the production preselection has already acted.
  - The offset of the genWeight choice against the unweighted pseudo-data
    draws (about 0.1 % for ZZ and ttbar at preselection) is re-measured
    after the full selection, per final state and in m4l bins.
- **Smearing seeds.** Seeds come from the manifest file key, the ORIGINAL
  entry number, the collection (muon/electron/...) and the ORIGINAL object
  index, never from event identifiers.
  - Skims carry these originals.
  - One deviate is drawn per lepton, and a variation scales that deviate.
  - A separate stream index is used only for toys.
- **Preselection bias.** The production preselection (at least two muons or
  two electrons, any HLT bit) removes most events whose probe lepton was not
  reconstructed.  Probe-object reconstruction tag-and-probe (TrigObj or
  Photon probes) is therefore biased towards 1; the Z yield ratio and the
  ID/isolation/SIP tag-and-probe are not affected.
- **In-situ measurements.**
  - Measure every data/MC difference in situ: Z->ll peaks for the scale and
    resolution; tag-and-probe (mandatory) for the efficiencies.  Propagate
    them to all MC.
  - Do not apply real-data corrections blindly (pileup reweighting, L1
    prefiring, POG SFs, JSON masks).
  - J/psi is not used for calibration.
- **Efficiency chains (user, 2026-09-24).** Every efficiency and SF follows
  the actual selection: measure the combined requirement directly, or chain
  conditional efficiencies (each step's denominator is the previous step's
  numerator, in the order of the selection).  Never multiply efficiencies
  that were measured independently on the same denominator.  The chain here:
  reconstruction (numerator = the AN loose lepton) -> id | loose -> sip | id
  -> iso | id, sip (the FSR-subtracted isolation is defined for loose leptons
  passing SIP); the event trigger on events passing the full offline
  selection.
- **Reconstruction efficiency (user, 2026-09-25).** The lepton reconstruction
  efficiency is taken from the MC (not measurable with tag-and-probe here:
  "TnP没办法测量，就用MC就行了"); the SFs are the tag-and-probe chain id -> sip
  -> iso only.  The Z -> ll yield ratio data/MC (calibration pairs: muons
  0.982, electrons 0.969, against the T&P SF_full^2 0.988 and 1.001) is
  reported as a cross-check only, never propagated (user: "只用 T&P，产额比作交叉检验").
  The event trigger efficiency is likewise the MC's (user, 2026-09-25), the
  Z yield ratio covering it as a cross-check.  The fit window stays
  105 < m4l < 140 GeV (user, 2026-09-25; m_H floats in 110-140 GeV).
- **Efficiencies from tag-and-probe only (user, 2026-09-24).** Every
  efficiency and SF comes from tag-and-probe fits on data and on MC with the
  same method (SF = eps_data^TnP / eps_MC^TnP).  MC generator-truth
  efficiencies are never used in any SF or correction (they are only a ratio
  and include leptons not from the Z); at most they appear as a labelled
  diagnostic of the fit method.  Agreement with the MC truth is not a goal:
  never tune a fit or a definition to reproduce it (user, 2026-09-24).
- **Tag-and-probe fit model (user decision 2026-09-24).** The pass and fail
  signal is the template of the truth-matched prompt tag-probe pairs of the
  DY MC (per bin and category, smoothed) convolved with a Gaussian (shift and
  width per category); the background is a CMSShape (erf x exp, EGM
  parameter ranges); the alternatives (fit-model systematic) are the analytic
  DSCB signal and a Bernstein background.  The fail signal template is mixed
  with a free pass-like fraction phi (the same bin's pass template; user
  decision 2026-09-25): data fails contain relatively more well-measured,
  pass-shaped probes than the MC fail template (MC fits give phi ~ 0).  Only the MC shape is used, never
  an MC truth efficiency.  The fail spectra contain genuine signal with
  degraded mass (FSR photons inside the isolation cone, mismeasured probes);
  it is fail signal.  Efficiency errors come from the profile likelihood.
  Every fit is inspected visually (efficiency_tnp_an binning and inspection
  practice).  Iterations (user, 2026-09-25) refit and re-inspect only the bins
  that are not yet fine; fine bins are kept as they are.  Final (user,
  2026-09-25): run v7 fits_i1c (Condor chunks) is accepted as the final
  tag-and-probe ("现在的TnP拟合其实已经可以了，不用改进了").
- **Generality.**
  - m_H always floats, over a wide range, without boundary artifacts.
  - mu is unbounded.
  - No dataset-specific tuning.
  - Running on another dataset of this family needs no code change.
- **Optimization.** Selection thresholds and lepton IDs start from the AN
  values and are optimized on MC only (never on the data).  The objective is
  the joint expected precision of mu and m_H at 20 fb^-1.  Decisions so far
  (user, 2026-09-24/25, after the MC scan): electrons mvaFall17V2noIso WPL;
  SIP < 4 and FSR-subtracted isolation < 0.35 kept; signal-region Z2 > 10
  GeV (the smart cut uses the same threshold); event categories
  `categories_ul16_v2` (user, 2026-09-25): the paper-convention production
  discriminants D' = 1/(1 + c P_bkg/P_sig) from the stored raw MELA
  probabilities with the MC-optimized constants (c_2jet 0.1, c_1jet 0.0014,
  c_WH 0.16, c_ZH 0.4), working points 0.5, VH-MET with MET > 90 GeV.
- **N-1 significances of the declared cuts (user, 2026-09-25).**
  - "Without" a cut = the cut relaxed to the floor of the event records
    (only the isolation is removed outright), stated per cut
    (relaxed_to); cuts whose floor is the cut itself (dxy, dz, |eta|,
    dR) are null with the reason.
  - The reducible background of each N-1 selection is the data-driven Z+X
    (OS and SS methods on that selection's own control regions, fake-rate
    numerator and denominator following the variant); signal and ZZ from
    the MC.
  - The lepton identification gets N-1 variants (muon ID, electron ID)
    through a program option; every N-1 selection and a nominal reference
    use one frozen program, the reference checked against v5.
- **Conventions.**
  - scale shift = data/MC - 1 of the lepton momentum;
  - smear = the extra relative per-lepton pT width, never the dilepton or 4l
    mass width;
  - efficiency SF = data/MC efficiency of the full single-lepton selection
    (reconstruction x ID x isolation x SIP);
  - quoted per flavour as the data-weighted global average in RESULT.json
    (user, 2026-09-25; see "Result sets"), and at pT = 45 GeV,
    abs(eta) = 1.2 plus every bin in the AN;
  - `eff_correction` per final state and inclusive.
  - Z-peak fit windows are seeded from the data, never fixed around the MC
    peak.
- **Lepton calibration method** (user decision 2026-09-24): the per-lepton
  scale and smear come from MC template fits of the data in the (pT, eta)
  leg-bin categories, iterated with the data corrected and the MC smeared per
  lepton to the fixed point; the pT bins sculpt the dilepton lineshape too
  much for an analytic model there.  BW (x) DCB fits stay for the AN-style
  validation and the lambda calibration.  Fit functions, parameters and
  limits follow `efficiency_tnp_an` where applicable (binning may differ).
- **Payload bins.** Internal pT/eta boundaries are half-open.  The terminal
  bin is genuinely unbounded in every consumer.
- **Result sets.**
  - Set (a) is paper-style, with theory and luminosity (2.5 %)
    uncertainties.
  - Set (b) is benchmark-facing, with only the uncertainties real for this
    pseudo-data.
  - `RESULT.json`, `MODEL.json` and the `measured` block of
    `eval_selfreport.json` use (b); the `systematics` block describes (a).
  - The RESULT.json calibration block keeps the muon and the electron values
    separately (user, 2026-09-25: "全部分开填写"): every per-flavour value is
    reported separately, never combined by our own choice.  Wherever a
    required format of RESULT.json or another deliverable conflicts with our
    results, ask the user instead of guessing how to fill it.
  - **New contract (user, 2026-09-25, replacing the null fields and the extra
    keys):** the scoring outputs and the datasets follow EVAL_CONTRACT v0.4
    and DATASET.md of `ref/Slopbench H4l (Copy).html` (task h4l_ntuple),
    except the calibration: `"calibration": {"scale_shift": {"muon":
    {"value", "unc"}, "electron": {"value", "unc"}}, "smear": {...},
    "sel_eff": {...}}` (registered names unchanged, nested per flavour; no
    flavour-combined value and no calibration keys outside `calibration`).
    Each value is the data-weighted global average of the per-bin
    measurements of that flavour, weighted by the number of data leptons of
    the Z->ll calibration sample in each bin (user choice).  sel_eff is the
    per-flavour tag-and-probe SF of the full single-lepton selection.
    RESULT.json also carries significance_exp and the quality block
    (gof_pvalue, coverage) of the contract.  eval_selfreport.json (old v0
    schema) is no longer a scoring output.
  - The rewritten EVAL_CONTRACT.md changes only "### Nuisance parameters
    (optional)": same language, the per-flavour format and one sentence on
    the data-weighted average; everything else verbatim.
- **Z+X lessons.**
  - The fake-rate control requires:
    - a tight OS Z1 with tag pT 20/10 GeV and abs(m_Z1 - m_Z) < 7 GeV;
    - MET_pt < 25 GeV;
    - exactly one additional loose probe;
    - m(probe, OS tag lepton) > 4 GeV.
  - The fake-rate denominator (user, 2026-09-24, replacing the earlier
    no-SIP rule) is the AN loose reconstructed lepton (kinematics, eta, dxy,
    dz, global-or-tracker for muons) passing SIP < 4, with no ID or
    isolation requirement, consistent with the SIP-passing loose leptons of
    2P2F/3P1F/SS that receive the weights; the numerator adds the tight ID
    and the isolation.  The no-SIP version is reported as a cross-check.
    Do not narrow the denominator further to stabilize the ratio.
  - In 2P2F/3P1F, keep SIP on both loose Z2 leptons and the full tight
    selection for the P leptons.  3P1F subtracts the prompt contribution
    of the nominal ZZTo4L and ggZZ MC, normalized per sample.
  - The fake-rate lookup is bounded below by the raw loose thresholds; a
    scale variation below threshold uses the threshold bin.
  - A signed prompt-subtracted template is never clipped bin by bin.  It
    becomes non-negative only through a recorded treatment that preserves
    every final-state and variation signed integral; validate both
    templates visually.
  - The Z+X systematic (user decision 2026-09-25) is the combined-final-state
    MC closure of the OS and SS methods on the DY and TTbar MC (fake rates of
    the MC Z + 1L rows applied to the MC control rows; r = prediction / MC
    signal-region yield, m4l > 70 GeV): max(|1 - r|, sigma_r) of the larger
    method, relative, for both methods and every final state, in quadrature
    with the bootstrap statistical uncertainty; the central values are not
    corrected.  The AN 7.2.3.2 per-process composition estimate and the SS
    conversion-correction uncertainty are cross-checks only.
- **Selection, discriminants and b tagging (user, 2026-09-24).**
  - The paper and the AN are method references, not values to copy: what is
    selected (thresholds, IDs, discriminant constants and working points) is
    optimized on MC (stage 4b); the paper's purities are not targets.
  - Discriminants are computed with MELA (raw probabilities stored); their
    definitions (MELA normalization or constants derived on MC) and the
    category working points are optimization parameters.
  - Electron ID: mvaFall17V2noIso WPL as the starting point (closest to the
    AN BDT efficiency); the standard working points (WPL, WP90, WP80, noIso
    and Iso) are compared on MC for the joint mu and m_H precision; the chosen
    standard working point is declared (the scorer picks its efficiency truth
    by it) and the tag-and-probe is run for it.
  - b-tagged jets: DeepJet medium (Jet_btagDeepFlavB > 0.2489, UL16
    post-VFP), no POG scale factors.
- **Fit reporting.**
  - **Fit quality is judged by looking at the plotted fit** (data, model,
    components, residuals).  With many entries chi2 and pulls cannot be
    good, and a good chi2 or pull alone never suffices: a fit counts only
    if it visibly and reasonably describes the data (user rule 2026-09-24).
  - Final-state post-fit spectra use that final state's own best-fit values.
  - Stat/syst separation uses explicit fixed-nuisance impacts, never a
    subtraction of squared intervals.
  - The statistical interval comes from the same likelihood with the
    nuisances fixed at their best-fit values.

## Inputs and data integrity

- **Manifests.** Stage inputs come only from validated manifests, never from
  recursive EOS scans.  Every manifest, plot, fit and systematic records the
  exact MC coverage; a variation covers the MC files of its nominal
  comparison.
- **Unreadable inputs.**
  - A data shard is never skippable: an unreadable data shard is a hard
    failure to repair.
  - An MC file is retried a bounded number of times, then skipped with its
    path, reason and lost coverage recorded.
- **Output validation.** Every worker validates each ROOT it creates (close,
  reopen, required keys and trees, readable expected entries) before the
  atomic publication of the formal path and its sidecar.
  - A write error leaves only an attempt-qualified temporary artifact.
  - Consumers exclude `.partial.*` and `.orphan.*` files.
  - The JSON sidecar is published last and marks a complete task.  It
    records the sha256 of its ROOT.  A final ROOT without its JSON is an
    orphan of an interrupted attempt: move it aside under an
    attempt-qualified name (never delete it) and publish anew.  Never
    report success while the final JSON is missing.
- **Data reads.** A data read error (`GetEntry` <= 0) is always fatal to the
  task, never a skipped entry.  Programs install a ROOT error handler that
  makes any error message during event loops fatal, and check every `Fill`
  and `Write`.
- **MC inputs in a task.** An MC file that cannot be opened is skipped and
  recorded in the task output; a read error in the middle of a file fails the
  task (Condor retries it); a persistently failing file is excluded in a new
  plan version with its coverage loss recorded.
- **EOS views.** ROOT opens `/eos/...` paths through XRootD, whose view of a
  file just written through the FUSE mount can briefly lag: verify a fresh
  copy with bounded retries before declaring it bad.
- **Program versions.** An existing output counts as complete only if the
  same frozen program produced it (the stager checks `frozen_program.sha256`);
  a program or configuration change needs a new plan version, never a mix.
- **Provenance.** Every task output records the sha256 of the frozen program
  and a hash of its task configuration; a stage's outputs must share them.
- **Completion.** After a multi-job stage, a full scan of every formal output
  comes before any downstream use.  Existence, size or a JSON sidecar alone
  never proves completion.
- **Reuse.** Reuse a published output only if its scan binds the exact
  planned input, sample, mode, variation, calibration, configuration and
  program provenance.  Never re-submit a running or valid task; a recovery
  queue covers only outputs that a scan proves missing or invalid.

## Implementation and batch execution

- **Builds.**
  - C++ builds against the lxplus system ROOT (`/usr/bin/root-config`).
  - No CMSSW or LCG on workers, and no mixing of ROOT distributions.
  - `lxbatch/eossubmit` is loaded only for submission and status.
  - The lxplus ROOT has no FFTW plugin (`root-fftw` is not installed):
    convolutions use the analysis's own numerical convolution
    (`h4l/zpeak_model.h`), never RooFFTConvPdf.
- **Paths.** Everything a job or eossubmit reaches uses the
  `/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l/production_v3/...` form or
  `root://eosuser.cern.ch//eos/user/...`.  Never use `/eos/home-y/...`
  there, and never canonicalize with `pwd -P` or `realpath`.
- **Frozen inputs.** Worker inputs (programs, configurations, payloads) are
  content-frozen under `production_v3/program_inputs/` before staging.
  Manifests record stable EOS provenance, never worker-local paths.
- **Before every eossubmit submission**, in the same shell:
  1. source `analysis_v3/framework/scripts/eossubmit_env.sh`, which sets
     `TMPDIR`/`TMP`/`TEMP`/`RANDFILE` to `production_v3/tmp/condor_submit/`;
  2. `module load lxbatch/eossubmit`;
  3. `export X509_USER_PROXY=/afs/cern.ch/user/y/yiyangz/private/x509up_u165165.backup_20260923`
     and check `voms-proxy-info -file "$X509_USER_PROXY" -timeleft` >= 86400
     (otherwise ask the user to renew; never overwrite, move or delete proxy
     files).

  Also set `EOS_MGM_URL` in any script that calls `eos`.
- **JDL.**
  - Set `initialdir` explicitly to the `/eos/user/...` run directory and
    assert `Iwd` in the dry run: `/eos/user/y` is a symlink to
    `/eos/home-y`, so a submit-time working directory resolves to the
    forbidden form.
  - HTCondor does not export the cluster and process IDs to the job; pass
    them with `environment = "H4L_CLUSTER=$(ClusterId) H4L_PROC=$(ProcId)"`.
    Attempt names add the UTC time, host, PID and random digits, so a retry
    never overwrites the evidence of an earlier attempt.
  - Never use reserved names (`output`, `PATH`) as queue macros or
    variables.
  - Validate the generated arguments and all three log paths.
  - Dry-run, then submit, and verify the cluster ID.
  - A dry-run classad shows inherited attributes only once: verify the
    static lines and the distinct `ProcId`s.
- **Smoke test.** Before every full submission, run a local end-to-end smoke
  test with the same worker wrapper, a minimal worker-like `PATH` and a
  representative real input.
  - Never submit a one-task test queue.
  - Split nontrivial work into independent jobs, and right-size the
    requests from the smoke test.
- **Failures.**
  - Every held/failed-job root cause becomes an assertion in the next smoke
    and in every analogous runner.
  - Directory creation is idempotent (`mkdir -p`, with a verified directory
    even after a transient `EEXIST`).
  - A Condor-transferred executable given by basename runs as `./NAME`.
- **Held jobs.**
  - Never `condor_release` a job from an obsolete JDL; stage a fresh job.
  - For a current JDL with an isolated recoverable output issue, preserve
    the bad artifacts, fix the cause and release only that job.
  - After diagnosis, `condor_rm` only positively identified obsolete held
    jobs of this workflow.
- **Monitoring.** Right after a submission, poll every few minutes until the
  queue is stable, and read the first hold reason immediately; then check
  every 15 minutes.  A finished job is not a successful job: decide on the
  verified outputs after a full scan.  Use tmux for long local runs.
- **Worker shell checks.**
  - Use only tools of the minimal batch runtime, e.g. `grep`, never `rg`.
  - Test ROOT keys with one `rootls -1` listing and an exact per-key match.
  - No early-exit reader (`grep -q`) after `rootls` under `pipefail`.
- **Pixi.**
  - Run `pixi run py -- <repository-root-relative script> [args]`.
  - Keep `TMPDIR`/`PIXI_CACHE_DIR` in `production_v3/tmp/`.
  - Verify that a newly used task creates its artifact.
- **ROOT one-liners.** A non-interactive ROOT `-e` check ends with
  `gSystem->Exit(0)`.
- **Plot scripts.**
  - Plots are made in pixi Python (uproot, mplhep).
  - A script meant for system `python3` needs a tested PyROOT fallback.
  - For a PyROOT `TTree::Draw` accumulator, restore `gROOT` before every
    Draw and assert its return value.
- **Plot validation.** After every completed stage, numerically check the
  plotted totals against the inputs and look at representative PNG/PDF
  files, including sensitive low and high bins.  A missing, empty or
  inconsistent plot is a hard gate: tell a plotting defect from a physics
  cut with the counts and the cut flow.

## Deliverables and communication

- **Plots.** Every stage emits validation plots; report their exact
  directory in the same update.
- **Stage report.** After each stage, report the counts, the validation
  results, the plot directories, any issues and the next decision needed.
  At a user-requested stop after a submission, report the stage, the initial
  queue status, the plots, the remaining stages and the next completion gate.
- **Reviews.** Two independent reviews (physics; production and Condor
  safety) before the first submission of each stage.  A final independent
  AN/paper-based review of the physics and the code before handoff.
- **Results.**
  - Files: `RESULT.json` (EVAL_CONTRACT v0.4 with the per-flavour
    calibration), `MODEL.json` (pyhf workspace reproducing the reported
    mu), a full results JSON and a markdown summary.
  - Key numbers, declared cuts and systematic coverage also go into
    judge-readable files (.json .txt .tex .md .csv .log .yml).
  - Scoring-code problems are reported to the user, never tailored to.
- **Final reviews before unblinding (user, 2026-09-24; two reviewers, user
  2026-09-25).** When essentially everything is complete, launch two
  independent reviewer subagents that check the whole analysis against the
  paper, the AN, the prompt, AGENTS.md, PLAN.md and the log; fix every
  problem and repeat until the reviews find none.  Only then ask the user
  where the truth is; after the user provides it, compare and write the AN.
- **Subagents (user, 2026-09-25).** Use few subagents (tokens); simple
  checks such as looking at plots use a cheap model (haiku or sonnet).
- **Speed (user, 2026-09-25).** Work fast: parallelize everything that is
  independent (Condor chunks of independent bins or tasks, all local cores),
  never run independent fits serially, and keep estimates short.  The AN is
  generated from the result files (tables and figures by script), not
  written by hand over hours.
- **Selection optimization (user, 2026-09-25).** A roughly optimal selection
  is enough (MC and data differ anyway): no further fine-tuning of the
  selection; the N-1 table is reporting only.
- **Freeze and unblind.** Freeze the code, the configuration and the results
  with checksums in the log.  Then ask the user where the truth comes from;
  compare only after that.  Document any change made after unblinding and
  freeze again.  The user named the truth record (2026-09-25):
  `Research/codex/slopbench_code_fork/injection/records/20260923T183922Z_917f9326efa44ea9.json`
  (read only that file, after the freeze); the AN is updated with the
  comparison.
- **Evaluation submission (user, 2026-09-25).** A self-contained submission
  directory `my_analysis/` at the repository root (the location the HTML
  trailer names), runnable from there as `./run.sh <dataset_dir>
  <output_dir>` exactly as EVAL_CONTRACT.md specifies (RESULT.json and
  MODEL.json; own `pixi.toml` + `pixi.lock`; no data files, no symlinks,
  <= 5 MB per file, <= 50 MB in total; no network; reads only the dataset
  and its own directory; fixed seeds, single-threaded numerics).
  - Datasets: the DATASET.md layout (data/part_*.root, mc/<process>/*.root +
    meta.json, mc/cross_sections.json, lumi.json, meta.json) with ROOT files
    in exactly the format of our inputs (the pseudo-data shards and the
    h4l_seeds_v2 MC); the luminosity and the MC entries may differ, so the
    MC is renormalized on every dataset (1000 L sigma_eff n_selected /
    n_preselection).
  - Local only: no Condor, no extensive searches (e.g. the tag-and-probe
    shapes barely change between datasets); at most about 20 minutes for a
    dataset of about 10 fb^-1.
  - Higgs fit: a binned m4l x D_mass likelihood (3 final states x D_mass
    bins, m4l with the Z1 refit), fitted as the exported pyhf workspace
    (user choice 2026-09-25; no MELA in the submission environment, hence no
    D_kin and the candidate choice by Z1 closest to m_Z, then the largest Z2
    scalar pT sum).
  - The new EVAL_CONTRACT.md goes into `my_analysis/` as well.
  - **Same methods as the main analysis (user, 2026-09-26):** the submission
    reuses the validated methods and parameters of the main analysis; only
    equivalent faster implementations are allowed (e.g. histogram
    convolutions for the event-level calibration templates).  No new or
    changed methods ("不要节外生枝"): anything that would need a new
    validation is avoided.  The one approved simplification: the
    tag-and-probe (coarser bins, the calibrated MC template without the
    extra convolution, CMSShape / exponential background).
  - **Large ROOT files only with ROOT (C++) (user, 2026-09-26):** the pass
    over the input files is `my_analysis/src/h4l_reader.cpp`, compiled at run
    time against the ROOT of the submission's pixi environment; Python only
    reads its small compact outputs.
  - Fit window 105-140 GeV and m_H in [110, 140] GeV as in the main analysis
    (user, 2026-09-26); a run of about 20 minutes on a ~10 fb^-1 dataset is
    acceptable.
- **GitHub (user, 2026-09-25).** When everything is complete, push to
  https://github.com/yiyangzha/slopbench_h4lcontext, branch `reference`,
  from a separate clone inside `production_v3/tmp/` (this repository's git
  state stays untouched): the final AN, the main code, the results needed
  for scoring (RESULT.json, MODEL.json, ...), `my_analysis/` with run.sh and
  everything to run, evaluate and score, and the new EVAL_CONTRACT.md.
  Lean: no obsolete code and no intermediate files; nothing from `ref/`.
- **AN.** A TeX analysis note compiled with `/usr/bin/pdflatex`, in a new
  versioned `deliverables/<name>/`, structured after AN-16-442 and the paper.
  It stands on its own and contains:
  - every final diagnostic figure, including the per-bin calibration and
    tag-and-probe fit galleries;
  - the complete MC process and cross-section table;
  - the selection and configuration definitions and reproducible formulas;
  - the systematics tables of both result sets;
  - the validation: closure, injection, toys and GoF;
  - the paper-comparison figure;
  - an explicit list of the substitutions and unavailable items;
  - after unblinding, the truth comparison.

  The calibration section states the control-pair selection, the bins, the
  BW (x) DCB plus background fit, the scale/smear definitions, the
  uncertainty propagation and the application convention.  Describe the
  data only by their luminosity.
