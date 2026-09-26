# H→4ℓ 2017 reanalysis experiment log

This is the live execution and planning record.  Update it with every plan
change, stage start/submission, job-status check, completion check, and stage
completion.  `AGENTS.md` remains the durable rule set.

## Overall plan

1. Build a certified 10/fb pseudo-data NanoAOD from modified MC only.
2. Build data and unmodified-MC manifests; make loose-v1 AN-like Z-peak lepton
   calibration (60–120 GeV `BW ⊗ double-CB +` smooth background).
3. Run AN-like loose-v1 H→4ℓ reconstruction with muon/electron thresholds 3/5 GeV;
   propagate lepton scale/smear variations.
4. Build reducible-background inputs and fit inclusive/three-final-state
   signal strengths, simplified peak-shift `mH`, and feasible fiducial or
   differential/category results, with validation plots and uncertainties.

## Current execution state

- 2026-07-30: Finalized the loose-v1 calibration pT binning after the user's
  statistical-precision refinement.  It uses five pT regions per flavor
  (3/5--20, 20--35, 35--50, 50--65, 65+ GeV) in every eta region: one fewer
  than the Figure-38 AN display but substantially finer than the old three
  regions.  The final exact minimal-PATH sidecar-wrapper smoke reopened all
  20 required histograms, had zero unreadable input, and its single real data
  calibration shard already had 227 entries in its least-populated bin.
  Rebuilt the sidecar, syntax-checked the analogous runner/stager, verified
  the certified 94/94 data and 75/75 DY input status, frozen selection/bins
  provenance, group-EOS outputs, personal-EOS logs, absence of `/eos/home-y`,
  6-GB classads, and four Condor dry-run classad files.  Submitted the four
  independent loose-v1 sidecars: `972124` data-muon, `972125` data-electron,
  `972126` unmodified-DY-muon, and `972127` unmodified-DY-electron.  Initial
  status is one idle process per cluster, zero held, each normalized to
  2 CPUs/6000 MB.  Wait for all four jobs, then fully reopen/validate all four
  outputs and check the full per-bin entry distributions before the Z fits.

- 2026-07-30: All four loose-v1 sidecar jobs ended normally (exit code zero;
  no hold) and their stderr logs are empty.  The independent full validator
  reopened every sidecar, required all 20 60--120-GeV spectra, and wrote
  `calibration_validation/calibration_mass_sidecars_loose_v1_20260730T1640Z.json`.
  The minimum populated spectrum is 20,515 data-muon tag/probe pairs, 25,387
  data-electron pairs, 56,762 DY-muon pairs, and 69,519 DY-electron pairs;
  the near-AN binning is therefore statistically sound.  The next active
  step is the local portable BW-convolved-double-Crystal-Ball plus smooth
  background fit with 3/5-GeV payload eligibility and per-bin diagnostic
  plots.

- 2026-07-30: Added the user's durable visual-validation requirement: every
  future completion check must inspect both formal ROOT/JSON outputs and its
  corresponding rendered figures.  For the completed loose-v1 calibration,
  all 40 per-bin PNG/PDF fit figures and four summary PDFs exist; the
  low-pT muon and high-pT/endcap electron PNGs were visually inspected and
  show the expected shifted pseudo-data and nominal-MC Z peaks with fitted
  curves, not an empty-plot defect.

- 2026-07-30: The loose-v1 calibration fit is payload-ready with all 20 muon
  and 20 electron bins valid and no failed payload bin.  Its payload and fit
  summary are respectively `calibration/lepton_calibration_loose_v1_20260730T1645Z.json`
  and `calibration/lepton_calibration_loose_v1_fit_summary_20260730T1645Z.json`.
  A worker-like, minimal-PATH reconstruction smoke used the frozen generic
  runner/program/configuration and a real 9.4/fb-manifest pseudo-data shard:
  it atomically published/reopened a valid ROOT with four candidates (three
  above 105 GeV), and its inspected cutflow figure is at
  `h4l_reco/validation/loose_v1_worker_smoke_20260730T1647Z/`.
  After static routing/provenance checks and ten JDL dry-runs, submitted the
  loose-v1 full reconstruction: cluster `973752` for the remaining 93 data
  files and clusters `973754`, `973763`, `973771`, `973779`, `973781`,
  `973784`, `973786`, `973788`, and `973790` for the 404 original-MC files in
  nominal, mu-scale up/down, electron-scale up/down, muon-smear up/down, and
  electron-smear up/down respectively.  Initial status is idle with no held
  job in every submitted cluster.  The old-selection reconstruction outputs
  remain historical and are excluded by the new calibration/configuration
  provenance.  Per the unattended-run instruction, check these queues every
  30 minutes; after completion run the formal ROOT scan, numerical/visual
  validation plots, then continue automatically.

- 2026-07-30: First unattended 30-minute reconstruction check: no held jobs.
  The data and the first seven MC-variation clusters have left the active
  queue; `973786` has 11 running, `973788` has 55 running and six normal
  completed processes, and the last variation `973790` has 28 running plus
  376 idle.  Continue the 30-minute polling cadence; no recovery action is
  warranted.

- 2026-07-30: Second unattended reconstruction check found all ten queues
  absent from `condor_q`; no held job required recovery.  The full C++ ROOT
  validation re-opened and provenance-checked all loose-v1 formal products:
  data 94/94 and every nominal/scale/smear MC variation 404/404 are valid,
  with zero missing or invalid files.  The reports are
  `h4l_reco/validation/reco_output_loose_v1_*_20260730T1720Z.json`.
  The nominal data cutflow has 51,131,174 input events, 20,729,573 HLT
  events, 1,837 four-lepton events, 481 H4l candidates, and 344 with
  m4l>105 GeV.  Numerical cutflow totals and the complete non-empty
  reconstruction validation figures were produced at
  `h4l_reco/validation/reconstruction_loose_v1_20260730T1728Z/`; both its
  cutflow and m4l PNG were visually inspected.  The m4l spectrum has the
  expected low-mass Z-rich region and the retained 125-GeV signal region, not
  an empty-plot failure.  Next is to implement/validate the approved AN-like
  MET fake-rate requirement and prompt-ZZ subtraction before staging the new
  loose-v1 reducible controls.

- 2026-07-30: The user approved retaining the AN/paper selection and pairing
  logic while making a controlled, versioned loose-v1 numerical adaptation for
  the low-statistics 9.4/fb pseudo-data: muon/electron object pT 3/5 GeV,
  leading/subleading 15/7 GeV, isolation 0.40, SIP3D 5, Z1/Z2 lower bounds
  30/8 GeV, and OS-pair mass 3 GeV.  Eta, dxy/dz, upper mass, FSR, smart-cut,
  and m4l requirements are unchanged.  The pT calibration bins retain
  near-Figure-38 AN-like granularity with a 3/5--20 first bin, 20--35,
  35--50, 50--65, and tail bins in each eta region.  Created immutable loose-v1
  analysis, fake-rate, and calibration-bin configurations.  The sidecar now
  freezes and records the exact analysis configuration, so the new 3/5-GeV
  payload cannot be confused with the old 5/7-GeV payload.  Existing 94/94
  data and 75/75 unmodified-DY calibration skims remain fully certified and
  reusable because they already retain both flavors above 3 GeV.  A
  minimal-PATH single-data-input worker-wrapper smoke passed with the rebuilt
  sidecar, reopening all loose-v1 mass spectra and reporting zero unreadable
  inputs.  Next action is to stage four new versioned data/MC x muon/electron
  mass-sidecars; no old-selection sidecar, fit, reconstruction, or reducible
  output will be reused as a loose-v1 product.

- 2026-07-30: Prepared (but did not submit) the next reducible-background
  `controls` queue.  The fresh plan has 94 certified 9.4/fb pseudo-data task
  identities; its exact minimal-PATH wrapper/binary smoke processed shard 0,
  reopened its required six ROOT objects, and published 544,897 read events,
  101,770 Z tags, and 209 fake-rate probes.  A new staging run correctly
  excludes that valid smoke output, leaving 93 disjoint pending tasks; its
  JDL/dry-run classads, personal-EOS logs, group-EOS outputs, and worker
  provenance have been checked.  An initial unsubmitted plan was preserved
  after its JDL audit found `/eos/home-y` configuration-source paths; the
  stager now publishes immutable `/eos/user/.../stable_sources/` copies before
  freezing the worker copies, and its plan hash includes that stable-source
  identity so old task JSON is never overwritten.
- 2026-07-30: Submission is deliberately paused at a physics-method gate.
  AN Sec. 7.2.1.1 and paper Sec. 7.2 specify the fake-rate control selection's
  `|m(ll)-mZ|<7 GeV` plus `MET<25 GeV` requirement; the current scanner has
  the mass condition but no MET cut.  Paper Sec. 7.2.1 also subtracts the
  unmodified-MC `ZZ` contribution from the OS 3P1F component, while the
  current config explicitly has `irreducible_control_subtraction.enabled:
  false`.  Both alter the reducible central value rather than merely adding a
  systematic.  No data/MC stage has been rerun and no controls Condor job has
  been submitted pending the user's method decision.
- 2026-07-30: The user approved adding both AN-like reducible corrections and
  requested a looser final selection, reusing existing code and submitting
  only genuinely necessary new work.  The corrections are feasible with the
  current inputs: every tested pseudo-data and original ZZ/ggZZ NanoAOD has
  `MET_pt`, and the four allowed prompt-irreducible samples total 188 public
  MC files with `genWeight` and `Runs:genEventSumw`.  The exact current
  94-file reconstruction has only 21 candidates in the 105--140-GeV fit
  window (4/4/13 in 4mu/4e/2e2mu); that is expected to be statistically small
  at 9.4/fb.  Selection changes are paused only for a necessary threshold
  definition: the payload is calibrated down to, but not below, muon/electron
  pT 5/7 GeV, so lowering those object minima would require a new calibration.
  Once a calibrated loose profile is specified, create a new versioned
  reconstruction target (data plus MC nominal/eight variations), never reuse
  the old-selection outputs, and retain the old outputs/plans unchanged.

- 2026-07-30: Completed the code-level implementation of the AN Eq. (7.1)
  OS reducible-background correction before any new controls were submitted.
  The fake-rate scanner now enforces the required `MET_pt <25 GeV` Z+probe
  control selection.  A separate immutable `prompt_templates` mode admits
  only the four unmodified public `ZZTo4L`, `GGZZ2E2Mu`, `GGZZ4Mu`, and
  `GGZZ4E` samples, validates each task against the frozen nominal-MC
  manifest, records file-level `Runs:genEventSumw`, and applies the same
  `f/(1-f)` 3P1F factor as the data term.  The merger normalizes each prompt
  sample with `1000*9.4*xsec/sum(genEventSumw)` and forms the final templates
  as data `3P1F-2P2F` minus the prompt 3P1F component, retaining both audit
  histograms.  It requires all four samples.  New collectors reopen every
  ROOT, check its task-bound provenance and every required key/binning before
  merging; data remains 94/94 mandatory while prompt MC may retain only the
  documented sub-10% shortfall.  The corresponding plotter produces PNG+PDF
  fake-rate maps, all systematic/final-state spectra, and a prompt-subtraction
  audit plot with numerical totals.  Both C++ binaries rebuild with no
  warnings and the worker/stagers pass bash/Python syntax checks.  Next:
  stage a fresh versioned 94-way data `controls` plan, execute one exact
  minimal-PATH local wrapper smoke, then submit only its remaining tasks.

- 2026-07-30: Staged and locally exercised the new versioned loose-v1
  reducible `controls` plan.  The exact minimal-PATH runner smoke processed
  pseudo-data shard 00000 in 4.11 s with 336 MB peak RSS, atomically reopened
  all six required ROOT objects, and recorded 544,897 input events, 49,554 Z
  tags, 258 probes, and the required 25-GeV MET threshold.  Its formal ROOT
  provenance now correctly names the final rather than temporary atomic path.
  The refreshed plan at
  `.../runs/reducible_controls/20260730T182522Z_2596457` safely reused that
  certified smoke output and has 93 remaining disjoint data shards.  The
  93-job final JDL was dry-run parsed with `RequestCpus=1`, `RequestMemory=6
  GB`, EL9, private submission temporary storage, personal-EOS logs, and no
  AFS or `/eos/home-y` task reference.  Submitted it as cluster `974840`;
  its immediate status is 93 idle, 0 held, normalized by eossubmit to two
  CPUs/6 GB.  Data completion remains mandatory: after the queue ends,
  fully reopen/provenance-check all 94 controls ROOTs, merge fake rates,
  render and inspect fake-rate figures, then stage the data templates and
  public prompt-ZZ templates together.

- 2026-07-30: The completed v1 controls ROOT scan itself passed 94/94, but
  visual inspection of the required muon fake-rate plot exposed a physics
  defect rather than a plotting defect: nearly every numerator/denominator
  ratio was unity.  The plot was correctly non-empty; its linear 3--1e6-GeV
  x-axis also made the finite-pT bins unreadable, so the validation plotter
  now uses a log pT axis.  Comparison with AN Secs. 3.1.1, 3.2.1, and
  7.2.1.1 found the cause: v1 had incorrectly narrowed the *loose* fake-rate
  denominator with ID/isolation/SIP requirements.  The AN loose denominator
  has only reconstructed-lepton kinematics, eta, dxy/dz (and global-or-
  tracker muons); SIP remains required only when applying 3P1F/2P2F Z2
  leptons, while P leptons remain fully tight.  This direct AN-correctness
  fix is versioned as `reducible_2017_loose_v2.json`; old v1 artifacts remain
  preserved and excluded.  The scanner and all related local checks rebuilt
  cleanly.  Next: execute a new v2 local worker smoke, inspect its fake-rate
  map, then submit the new 94-file data controls queue.

- 2026-07-30: The v2 exact minimal-PATH smoke reopened every formal object in
  3.87 s/335 MB and now has 3,689 loose probes (rather than 258) in the first
  pseudo-data shard.  Its independently merged test has 1,520 muon and 2,169
  electron denominator probes, with 164/39 tight numerators; its inspected
  log-pT maps at
  `reducible_background/validation/controls_loose_v2_smoke_20260730T2105Z/plots/`
  show physical O(0.01--0.3) populated-bin rates, rather than the v1
  near-unity artifact.  The plotter's map axis is now logarithmic so the
  actual 3--100-GeV bins are visible despite the overflow-tail bin.  A fresh
  93-way plan reused that certified smoke output, passed JDL and dry-run
  checks (EL9, one CPU/6 GB request, no invalid paths), and was submitted as
  Condor cluster `974844`.  Immediate status: 93 idle, zero held.  On
  completion, require the v2 94/94 ROOT/provenance collector, merge/inspect
  the complete fake-rate maps, then stage the two template queues.

- 2026-07-30: Diagnosed the apparently empty reconstruction validation figures
  from `reconstruction_20260730T155019Z/`.  They were not caused by a strict
  H→4ℓ cut: the certified 94-file, 9.4/fb pseudo-data reconstruction has
  51,131,174 input events, 20,729,573 HLT events, 990 four-loose-lepton
  events, 309 H4l candidates, and 259 with `m4l > 105 GeV`.  A minimal PyROOT
  reproduction showed that opening each input TFile changed the current ROOT
  directory, so `TTree::Draw` filled a temporary same-named histogram while
  the retained plot histogram stayed zero.  Restoring `gROOT` before every
  Draw and asserting its return fixed both data and MC histogram accumulation.
  The fresh, inspected non-empty validation set (including the data cutflow)
  is at
  `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/h4l_2017_nanoaod_analysis/h4l_reco/validation/reconstruction_20260730T155019Z_fixed_20260730T1623Z/`.
  The cutflow plot was additionally made compatible with system Matplotlib
  3.4, and the analogous future reducible-background plot now has a
  system-Python PyROOT/plain-Matplotlib fallback that passed import/CLI smoke;
  it still requires its mandatory real-input plot smoke before publication.
  No reconstruction data,
  MC, calibration, or physics selection was changed or resubmitted.

- 2026-07-27: Diagnosed the original mixer as sampling the shard-front entries
  repeatedly.  Replaced it with a global deterministic affine permutation:
  no replacement in each complete source pass.
- 2026-07-27: Frozen `10fb_global_v3` plan at
  `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/h4l_2017_nanoaod_analysis/fake_data_nanoaod/10fb_global_v3/metadata/mix_plan_global_v3.json`.
  It uses the already audited 404/404 modified source files, has 100% source
  coverage, and targets 54,497,625 events.  Modified MC remains mixer-only.
- 2026-07-27: Submitted and completed shard-0 preflight cluster `967753`.
  It exited 0 and produced a validated data-only ROOT/JSON pair with 544,897
  events.  ROOT/XRootD warnings were retried successfully; no input was
  omitted.
- 2026-07-28: Submitted 32-core all-shards cluster `968102`; it remained idle
  because eossubmit normalized the request to 32 CPUs/96 GB.  User approved
  replacement; it was removed before execution.
- 2026-07-28: Replaced the oversized job with 20 disjoint batch jobs, five
  shards each and two worker threads per job.  Range parsing, range-resume on
  completed shard 0, C++ compilation, shell syntax, frozen personal-EOS
  inputs, and the 20 exact non-overlapping ranges `[0,100)` were checked.
  Each job requests 2 CPUs and 6 GB; no data/MC event selection was relaxed.
- 2026-07-28: Submitted the 20 batch jobs as Condor cluster `968104`.  Initial
  status: all 20 are idle, none held; every classad has the intended
  `RequestCpus=2`, `RequestMemory=6000`, and its assigned disjoint range.
- 2026-07-28: Cluster `968104` exceeded its 6-GB cgroup limit (one recorded
  peak was 10.577 GB) and was removed in full by explicit user instruction.
  The replacement plan is 50 disjoint two-shard batches, one worker per batch
  and an explicit 16-GB request; it retains full `[0,100)` pseudo-data
  coverage and all frozen MC-source use.
- 2026-07-28: The removed jobs left 85 atomically published shard ROOT files:
  41 have JSON sidecars and 44 need only their embedded-summary sidecar
  restored.  The new 50-range submission includes every range explicitly;
  completed shards are resumed without overwrite and remaining shards are
  produced normally.
- 2026-07-28: Submitted the 50 one-worker/16-GB ranges as Condor cluster
  `968107`.  eossubmit normalized each classad to 6 CPUs and 18 GB to satisfy
  its memory-per-CPU policy; the mixer itself retains one worker for the
  two-shard range.  Initial status: all 50 idle, none held.
- 2026-07-28: Held job `968107.5` was an isolated output-resume failure, not a
  memory issue: old `shard_00010.root` was unreadable while shard 11 resumed.
  The unreadable ROOT and its JSON were preserved with
  `unreadable_resume.968107_5.20260728T071308Z` suffixes, then this one
  current-JDL job was released.  Its range will regenerate shard 10 and
  resume shard 11; no job was removed or resubmitted.
- 2026-07-28: Diagnosed 27 held jobs in cluster `968107`.  This is neither a
  Condor memory hold nor a general batch-routing fault: each fails on the same
  modified TTBar source `40000/6AF3BAC7-B767-044A-8801-4A81E3A80BF4.root`.
  Independent direct tests show unstable EOS/XRootD readability (zombie or
  `kRecovered` handles interspersed with successful opens); one worker reached
  the resulting ROOT basket abort (exit 134).  The existing frozen worker
  accepted recovered handles, so it is incompatible with a safe release.
- 2026-07-28: Started a non-destructive v4 recovery: preserve all v3 products,
  explicitly exclude this one unreliable modified-MC mixer-only input with its
  reason recorded, strengthen the input probe with representative Events reads,
  and reject/retry recovered ROOT handles in the fresh worker.  A new complete
  v4 plan/output will be used; no v3 product will be overwritten.
- 2026-07-28: The first v4 audit used `TTree::GetEntry` with all branches
  enabled, unnecessarily loading generator baskets while checking data-only
  mixer inputs.  Before it wrote any manifest, it was replaced by the same
  three-point integrity check on the required `run`, `luminosityBlock`, and
  `event` data branches only.  This keeps the demonstrated `run`-basket
  failure detectable while making the 404-file control step practical.
- 2026-07-28: Resource audit of completed `968107` one-worker/two-shard jobs
  found only 0.96--0.99 GB peak `MemoryUsage`; their 18-GB/6-CPU eossubmit
  allocation was excessive and not the cause of the holds.  The replacement
  will use 50 one-worker jobs with a 6-GB request (normalized to 2 CPUs/6 GB),
  retaining a large measured margin while improving scheduling availability.
- 2026-07-28: Added a source-audit hard gate for automatic unreadable-file
  exclusion: no more than 5% unavailable files per sample and 2% globally.
  Each allowed omission is recorded in the v4 manifest/plan; a larger failure
  cannot silently produce a pseudo-data plan.
- 2026-07-28: v4 audit completed with 399/404 readable modified sources.  The
  five exclusions are recorded with reasons; the per-sample maximum is 2/75
  (DY) and the global loss is 5/404 = 1.24%, both below the hard gates.  v4
  still targets all 54,497,625 Poisson-drawn 10/fb events from its complete
  frozen *available* source catalogue.  The old builder process predated the
  hard-gate code, so the current staging script independently enforces those
  same limits before it can submit this v4 plan.
- 2026-07-28: Staged and inspected the v4 replacement JDL at
  `/eos/user/y/yiyangz/h4l_2017_nanoaod_analysis/condor/fake_data_mixing/runs/20260728T092748Z_full/submit.sub`.
  It has exactly 50 non-overlapping ranges `[0,2) ... [98,100)`, one mixer
  worker per range, 6-GB memory request, personal-EOS logs, and newly frozen
  v4 plan/worker inputs.  The stager's independent unreadable-source limits
  passed before the JDL was written.
- 2026-07-28: Submitted the inspected v4 JDL as Condor cluster `968142`.
  Initial check: all 50 ranges are idle, none held; eossubmit normalized the
  resource classads exactly as intended to 2 CPUs/6 GB per range, while each
  wrapper runs one mixer worker.  v3 output remains preserved and excluded
  from all downstream inputs.
- 2026-07-28: Diagnosed 12 held v4 ranges.  This is a common source-catalogue
  failure, not an isolated range issue: repeated DY failures use
  `DYJetsToLL/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F`, with corrupted
  basket aborts also involving
  `DYJetsToLL/280000/D43DBBA7-B019-7F4F-9632-07A0817D785B` and
  `TTBar/40000/1DA1E53F-DCA9-BE49-BAE2-7AA2D6FDA7B4`.  Two further ranges
  (`39`, `47`) independently exceeded the 6-GB cgroup limit (10.17--11.01
  GB).  At diagnosis, only range `12` remained running; 99 ROOT and 86 JSON
  products existed, but cannot be combined with a revised source plan.
- 2026-07-28: The remaining original v4 job finished its published shard, so
  99 ROOT and 87 JSON products now exist.  To preserve those completed v4
  outputs and their global affine selection, the recovery does **not** create
  a v5 plan.  Added the frozen
  `fake_data_nanoaod/config/source_skip_global_v4.json` registry for the
  three newly diagnosed unreadable modified inputs.  The mixer now skips only
  a listed (or ordinary read-failing) file's already-selected entries without
  replacement, recording requested/written/skipped counts and path/reason;
  it never shifts any completed deterministic draw or replaces a data-like
  event.  Future plan-time probes now read every non-MC Events branch in an
  isolated process, rather than three ID points, to find comparable damaged
  baskets before mixing.
- 2026-07-28: Compiled the mixer/finalizer/validator successfully, checked
  shell syntax, registry-to-v4-manifest membership, and generated an
  inspected fresh recovery JDL.  It selects only the twelve old held ranges
  (`4--6`, `8--10`, ..., `98--100`); each wrapper resumes any existing shard
  and only creates missing ROOT/JSON artifacts.  The JDL transfers the
  content-frozen skip registry, writes logs to personal EOS, requests one
  worker and 12 GB, and keeps scientific output in group EOS.
- 2026-07-28: Removed exactly the twelve confirmed obsolete held processes
  from `968142` after retaining their diagnostics.  Submitted the inspected
  recovery JDL as cluster `968520`.  Initial status is 12 idle / 0 held;
  eossubmit normalized the 12-GB request to `RequestCpus=4`,
  `RequestMemory=12000` per range.  No completed shard was requeued and no
  scientific ROOT/JSON artifact was removed.
- 2026-07-29: Re-reviewed the AN and paper before the next gate.  The pending
  pseudo-data finalization does not alter their lepton calibration or 4-lepton
  method: the subsequent analysis remains anchored to 7/5-GeV e/µ thresholds,
  20/10-GeV candidate leading/subleading requirements, Z1/Z2 construction,
  three final states/categories, and the documented signal-strength, simplified
  mass, and feasible fiducial outputs.
- 2026-07-29: Verified 100/100 non-empty *formal* v4 ROOT+JSON shard pairs;
  the additional twelve `.partial.*.root` files are retained interrupted
  attempts and cannot match the finalizer's exact shard names.  Aggregated
  accounting is requested=54,497,625, written=54,402,728,
  skipped-unreadable=94,897, with all recovery stderr files clean.  The
  recorded skips are deliberately not replacement draws.
- 2026-07-29: Before staging the single finalization/validation job, fixed the
  finalization and analogous revalidation stagers to content-freeze every
  worker binary, wrapper, plan/manifest in personal EOS.  Both now request
  one worker with a 6-GB memory request; finalization additionally rejects a
  missing/empty formal pair before a job can be submitted.
- 2026-07-29: Re-read the durable execution rules and AN/paper method before
  the v4 completion gate.  Recompiled the finalizer and syntax-checked every
  analogous finalization/revalidation wrapper.  The finalizer now compares its
  transferred frozen plan byte-for-JSON with the stable group-EOS plan and
  publishes that stable EOS path, rather than an ephemeral worker-local path,
  in the downstream dataset manifest.  The formal-pair preflight again found
  exactly 100 non-empty ROOT+JSON pairs; the plan is complete over its frozen
  available source catalogue.
- 2026-07-29: Staged, inspected, and submitted the single v4
  finalization/validation task as Condor cluster `968763`.  Its JDL contains
  only personal-EOS frozen runner/binary/plan inputs; its explicit stable
  group-EOS plan provenance and all three personal-EOS log paths were checked
  before submission.  Initial status is `idle` (0 held); eossubmit has
  normalized the submitted 1-CPU/6-GB request to `RequestCpus=2`,
  `RequestMemory=6000`, while the wrapper itself executes one finalizer and
  one validator serially.
- 2026-07-29: Cluster `968763.0` is held with exit code 1, not for memory or
  transfer: its finalizer observed v4 shard 71 as a ROOT zombie and correctly
  refused to publish any finalization artifact.  A subsequent identical
  full-read preflight instead reached shard 57 and found it zombie, while
  `rootls` can intermittently list shard keys.  Historical mixer stderr has
  `TFile::Flush ... No such device` for formal shards 52, 57, 71, and 87.
  The old mixer then renamed the temporary ROOT and wrote its sidecar without
  a close/reopen check, so these four formal pairs are suspect even though
  their files are non-empty.  No final manifest, completion certificate, or
  validation ROOT/JSON was created.  Keep the finalization job held; do not
  release or reprocess until the user approves the four-shard repair.
- 2026-07-29: Added the user's mandatory integrity gates to `AGENTS.md`:
  certify each new non-official input before its first submitted use (while
  reusing immutable prior certifications and exempting official Open Data),
  have each worker close/reopen/check its ROOT before publication, and run a
  full formal-output ROOT scan after every multi-job stage before submitting a
  consumer.  Workers must bounded-retry and then record/skip an isolated
  unreadable independent input rather than hold an entire queue.  The explicit
  exception is a required pseudo-data shard, which must be regenerated from
  its frozen plan rather than silently omitted from the required 10/fb data.
- 2026-07-29: User approved the explicit v4 exception: do not regenerate
  suspect formal shards 52, 57, 71, and 87.  Exclude them through a frozen
  output-skip registry and set the pseudo-data analysis luminosity solely by
  file coverage, `10 fb^-1 * 96/100 = 9.6 fb^-1`; do not derive a finer
  event-count luminosity.  This supersedes the prior v4-specific
  regeneration gate but not the general output/integrity requirements.
- 2026-07-29: Implemented the frozen `output_skip_global_v4.json` registry,
  finalizer/validator propagation of formal-file coverage and effective
  luminosity, and the matching README command.  The finalizer now retries a
  formal ROOT three times and fully reads every Events entry of every retained
  shard before publication; an independently unreadable formal input is
  explicitly recorded as another runtime output skip instead of holding the
  queue.  The obsolete held JDL `968763.0` was removed after its diagnostics
  were retained.  The fresh pre-submission input check passed all 96 retained
  ROOT headers/required keys, skipped exactly 52/57/71/87 by registry, and
  found no additional warning.
- 2026-07-29: Submitted the inspected replacement finalization/validation
  task as Condor cluster `968769`.  It uses only frozen personal-EOS
  runner/binary/plan/output-skip inputs, logs exclusively to personal EOS,
  and its arguments retain stable group-EOS plan provenance.  Initial status
  is `idle` with 0 held; eossubmit normalized the 1-CPU/6-GB request to
  `RequestCpus=2`, `RequestMemory=6000`.  The worker serially scans all Events
  entries of the 96 retained shards before it can publish any manifest.
- 2026-07-29: Completion verification passed: `968769` exited 0 after a full
  retained-shard scan, used 866 MB peak memory, and published mutually
  consistent manifest/completion/validation artifacts for 96/100 files,
  9.6/fb, and 52,220,289 events.  The validation ROOT reopens with all six
  required objects.  Recorded the durable Pixi repository-root script-path
  rule in `AGENTS.md`; the final validation plotting command created all five
  expected PDFs below
  `10fb_global_v4/validation/plots/` (the first immediate EOS directory check
  raced the asynchronous EOS directory visibility, but later inspection
  verified all five non-empty PDFs).
- 2026-07-29: Fixed pseudo-data manifest construction to use only the exact
  `mixed_data_completion.json:data_files` list, rather than recursively
  enumerating a shard directory that also contains retained partial attempts
  and four excluded formal ROOTs.  Published and checked the new central
  `manifests/data_9p6fb.json`: 96/96 certified files, 9.6/fb, and exact
  completion-certificate binding.  The existing immutable `nominal_mc.json`
  was inspected rather than overwritten; it matches the current 11-sample
  public-MC catalogue and records 404/404 available official files.
- 2026-07-29: Strengthened the next calibration-skim stage before it is first
  submitted: each C++ output now closes/reopens and fully reads its retained
  tree before publication; the worker retries a classified unreadable input
  and records a skip with exit 0 after three failures; resume checks validate
  existing output keys; and the stager freezes runner/binary/task JSON in
  personal EOS.  Compiled and shell/Python syntax-checked all analogous
  calibration scripts.  A dry DY calibration test JDL was staged at
  `.../calibration_skim/20260729T073614Z_3553672/submit_test.sub`: one
  unmodified official DY file, personal-EOS frozen worker inputs/logs, group
  EOS output, and a 1-CPU/4-GB request.
- 2026-07-29: Submitted that inspected DY calibration test as Condor cluster
  `968797`.  Initial status is idle with 0 held; eossubmit normalized the
  request to `RequestCpus=2`, `RequestMemory=6000`.  It is a single official
  unmodified DY input and cannot alter pseudo-data or any completed output.
- 2026-07-29: Cluster `968797.0` exited normally in four seconds with no
  memory/resource issue, but it did not execute the skim binary: the frozen
  worker used unavailable `rg` for ROOT-key checking.  A second audit found
  that its ordered multi-key regexp was also intrinsically invalid because
  ROOT does not guarantee key order.  The worker and stager are now repaired
  to use minimal-runtime `grep`, `rootls -1`, and independent exact key tests.
  The worker had moved the pre-existing formal DY skim aside because of that
  false negative; its retained copy has all required keys and a PyROOT full
  scan read all 638,344 `CalibrationEvents` entries.  No replacement job has
  been submitted, no old artifact has been deleted, and the test's
  `skipped_unreadable_input` summary is known to be a checker false negative.
- 2026-07-29: User replaced Condor one-task tests with a mandatory local
  end-to-end smoke gate before every full queue.  All calibration, mixer, and
  reducible-background staging interfaces/docs no longer create a test queue;
  the calibration smoke runner uses the actual worker wrapper under a minimal
  `/usr/bin:/bin` path, writes only an isolated `calibration_smoke/` output,
  and fully rereads it.  The root-key logic now lists each ROOT once and makes
  independent exact checks, avoiding both unavailable `rg`, key-order
  assumptions, and `pipefail`/`grep -q` SIGPIPE false negatives.
- 2026-07-29: Restored the wrongly preserved formal DY output only after its
  full 638,344-entry scan passed.  The corrected local DY smoke test completed
  with the same official input, 638,344 fully reread output entries, complete
  schema, correct `DYJetsToLL`/`mc` provenance, and no retry.  The fresh full
  DY staging run
  `.../calibration_skim/20260729T081619Z_3816032` was syntax/JDL-routed
  checked and has `planned_outputs=75`, `pending_outputs=0`: all 75 formal
  DY skims are present with the required top-level keys, so there is no
  missing task that Condor can or should submit.  No existing scientific ROOT
  was overwritten or reprocessed.

## Next gate

Wait for the current 94-way pseudo-data nominal reconstruction queue to finish.
Require 94/94 completed records and a full reopen/key/event-read scan of every
formal reconstruction ROOT before generating the data cutflow/validation plots
or staging a background consumer.  Then perform the same local-smoke and
coverage-aware procedure for the unmodified nominal MC plus its eight matching
lepton scale/smear variations.

## Active scan gate

- 2026-07-29: Implemented a C++ full calibration-skim scanner which forks one
  child per ROOT, checks all required objects, and calls `GetEntry` on every
  `CalibrationEvents` entry with all branches enabled.  A crashing/unreadable
  child is recorded as an isolated invalid MC file, rather than holding the
  aggregate job.  The real-worker local smoke passed on one DY skim
  (638,344/638,344 readable entries); a second local smoke deliberately added
  one nonexistent file and correctly produced valid=1/2, invalid=1 without
  stopping the scan.  The fresh frozen, personal-EOS staged 75-file scan JDL
  is `.../calibration_output_validation/20260729T082645Z_3853140/submit.sub`;
  its group-EOS report target is
  `calibration_validation/skim_full_scan_DYJetsToLL_20260729T082645Z_3853140.json`.
  The JDL has one queue item, personal-EOS executable/inputs/logs, group-EOS
  report, 1 CPU/4 GB request, and no AFS worker reference.
- 2026-07-29: Submitted the inspected full 75-file DY scan as Condor cluster
  `968808`.  Initial state is idle with 0 held; eossubmit normalized the
  request to `RequestCpus=2`, `RequestMemory=6000`.  It is an integrity gate
  only and cannot overwrite any skim: it publishes a new timestamped JSON
  report after the isolated per-file scans finish.
- 2026-07-29: User made parallelism a durable execution rule: after a local
  smoke gate, split every nontrivial safely independent workload across
  multiple Condor jobs; reserve a one-job submission for genuinely small or
  non-partitionable work.
- 2026-07-29: Cluster `968808` completed normally in 16m57s with 231 MB peak
  memory and no stderr.  Its formal report passed the full gate: 75/75 listed
  DY skims valid, zero invalid, 100% coverage, and 31,559,821 fully readable
  `CalibrationEvents`.  This scan was already in flight before the new
  parallelism rule and is not repeated.  The next coverage aggregation is
  being redesigned as parallel file batches plus a merge, rather than the old
  single aggregate worker.
- 2026-07-29: Replaced the next coverage execution with 75 independent
  one-skim C++ coverage jobs.  Each worker uses the fully scanned DY report as
  its frozen one-file input, retries locally, validates its own partial ROOT,
  and records a failed input as a task-summary skip rather than holding the
  queue.  The local minimal-PATH worker smoke produced the expected seven
  ROOT objects and a 1/1 coverage JSON.  The inspected fresh run
  `.../calibration_coverage_parallel/20260729T093441Z_141954` has 75 disjoint
  frozen inputs, personal-EOS program/log paths, group-EOS partial outputs,
  one CPU/4 GB per task, and no Condor test queue.  A small merge will be
  staged only after these independent parts are verified.
- 2026-07-29: Submitted the inspected parallel DY coverage queue as Condor
  cluster `969544`.  Initial status is 75 idle / 0 held; eossubmit normalized
  every task to `RequestCpus=2`, `RequestMemory=6000`.  Each task has exactly
  one distinct fully scanned calibration ROOT and timestamped group-EOS
  partial outputs, so no task can overwrite another.
- 2026-07-29: Cluster `969544` completed with 75/75 `completed` task
  summaries and no nonempty stderr.  The full post-stage check found exactly
  75 ROOT + 75 JSON partial artifacts, every ROOT has all seven required
  coverage objects, and summed JSON counters are readable=75,
  unreadable=0, events=31,559,821.  The worker's minimal-PATH local merge
  smoke correctly combined two parts.  The inspected merge staging run
  `.../calibration_coverage_merge/20260729T141550Z_1503126` freezes its
  binary/runner/two lists under personal EOS and uses all 75 validated pairs;
  its one-job merge is deliberately small and non-partitionable.
- 2026-07-29: Submitted that inspected small merge as Condor cluster `969598`.
  Initial state is idle with 0 held and eossubmit-normalized
  `RequestCpus=2`, `RequestMemory=6000`.  Its timestamped final coverage ROOT
  and JSON paths are recorded in the merge staging summary; after it succeeds,
  they will be independently validated and plotted.
- 2026-07-29: Cluster `969598` completed normally with empty stderr.  Its
  final coverage JSON reports 75/75 readable DY skim parts, zero unreadable,
  100% coverage, 31,559,821 calibration events, 14,912,000 muon Z pairs, and
  12,885,737 electron Z pairs.  Independent ROOT rereading confirmed all
  seven required objects.  Four one-page PDF checks (both Z peaks and both
  probe `pT`--`|eta|` maps) were generated locally with system Python/PyROOT,
  byte-for-byte copied and reread from group EOS at
  `calibration_validation/plots_DYJetsToLL_20260729T141550Z_1503126/`.
  Direct Pixi plotting was not used because its interpreter is on EOS and
  blocked in filesystem I/O; the plot script's system-Python PyROOT fallback
  was verified instead.
- 2026-07-29: Began the certified pseudo-data calibration-skim stage.  The
  immutable `data_9p6fb.json` input passed its prior full-completion gate
  (96/96 files at 9.6/fb).  The actual minimal-PATH calibration worker smoke
  on shard 00000 selected and fully reread 165,419 calibration events from
  544,897 NanoAOD events.  The inspected full queue at
  `.../calibration_skim/20260729T143110Z_1569978` contains 96 disjoint
  pseudo-data shards, 96 unique group-EOS outputs, frozen personal-EOS
  program/task inputs and logs, and a 1-CPU/4-GB `el9` JDL.  It is ready for
  the formal parallel submission; data completion remains a hard 96/96 gate.
- 2026-07-29: Submitted that inspected 96-way pseudo-data calibration-skim
  queue as Condor cluster `969599`.  Initial query is 96 idle, 0 held; as
  expected eossubmit normalized each request to 2 CPUs and 6,000 MB.  Logs
  are exclusively under the frozen personal-EOS run directory.  Wait for all
  96 jobs to finish before performing the required full output scan or any
  pseudo-data calibration consumer.
- 2026-07-29: The user explicitly approved proceeding without the one
  outstanding calibration-skim task, `969599.54`, which maps uniquely to
  source/data shard `fake_data_10fb_shard_00055.root`.  The 95 already
  published formal calibration ROOTs have empty stderr logs.  This is a
  global analysis-selection change, not merely a plotting-label change: the
  immutable v4 source certificate remains 96/100 files at 9.6/fb, while all
  downstream data consumers use the newly published, complete-relative-to-
  selection `manifests/data_9p5fb.json` with exactly 95 files, shard 55
  excluded, and luminosity 9.5/fb.
- 2026-07-29: Extended the manifest builder with a guarded certified-subset
  mode.  It accepts only a strict subset of the source completion list, needs
  an explicit reason, and enforces luminosity equal to source luminosity times
  selected/certified file count.  A local 95/96 manifest smoke test and the
  published 9.5-fb manifest both passed.  `AGENTS.md` now makes this selection
  and its luminosity mandatory across calibration, reconstruction, background,
  templates, fits, and plot labels.
- 2026-07-29: Replaced the pending pseudo-data full-skim scan with a parallel
  one-ROOT scan/merge implementation.  A local minimal-PATH worker test fully
  read data calibration shard 00000 (165,419 entries), and a local one-report
  merge test passed.  The inspected fresh 95-way staging run
  `.../calibration_output_validation_parallel/20260729T152324Z_1760747`
  selects exactly the new manifest's 95 outputs, excludes shard 55, freezes
  program/list inputs on personal EOS, gives each task a unique group-EOS
  report and personal-EOS log, and records that the later merge must require
  95/95 valid outputs.
- 2026-07-29: Submitted the inspected 95-way pseudo-data calibration full
  scan as Condor cluster `969604`.  Its initial state is 95 idle and 0 held;
  eossubmit normalized every 1-CPU/4-GB request to 2 CPUs and 6,000 MB.  The
  deliberately ignored old calibration task `969599.54` is still running, but
  it is not a member of the frozen 9.5-fb manifest or the new validation
  queue, so it cannot enter any downstream input list.  Wait for all 95 scan
  reports, then run the verified merge and require 95/95 valid before the
  next calibration consumer.
- 2026-07-29: Per the user's file-count rule, removed the sole remaining
  unfinished/held validation task `969604.1`; it was the one-file scan for
  calibration shard 00001 and had not produced a report.  The other 94
  reports are completed.  Published the immutable certified-subset manifest
  `manifests/data_9p4fb.json`, excluding shards 00001 and 00055 from the
  96-file/9.6-fb source certificate.  Its enforced luminosity is
  `9.6 * 94/96 = 9.4 fb^-1`; it is now the only downstream data input.
- 2026-07-29: Added the user's held-prevention rule to `AGENTS.md`: each
  observed hold root cause must be tested locally under the worker-like PATH
  and corrected in every analogous runner/stager before another submission.
  The current one-file validation runner already has the required idempotent
  `mkdir -p` plus existing-directory check; the held task used the older
  frozen wrapper.
- 2026-07-29: Re-reviewed the calibration portions of the AN/paper before
  staging the Z-mass sidecars: retain the 5/7-GeV muon/electron object cuts,
  Z tag/probe range, muon FSR recovery with isolation subtraction, no
  electron-FSR approximation, and the subsequent `BW ⊗ double-CB +` smooth
  background fit.  The newly published data availability certificate is
  94/94 calibration inputs at the user-approved 9.4/fb; the unmodified-DY
  certificate is 75/75 inputs.  Both input lists were checked before staging.
  A clean local end-to-end sidecar smoke used the exact Condor basename
  transfer convention, `/usr/bin:/bin` PATH, and already-existing group-EOS
  output/summary directory; it reopened the ROOT metadata and all 12 mass
  histograms and required JSON `unreadable_outputs=0`.  Four distinct frozen
  JDLs were staged and parser dry-run checked, with unique group-EOS outputs,
  personal-EOS logs/frozen inputs, no AFS or `/eos/home-y` worker reference,
  one CPU/6-GB request, and `on_exit_hold` for a true nonzero worker exit.
- 2026-07-29: Submitted the four inspected calibration-mass sidecars as
  independent Condor clusters: `969611` (data muon), `969612` (data electron),
  `969613` (unmodified-DY MC muon), and `969614` (unmodified-DY MC electron).
  The immediate eossubmit status check reports all four as idle, none held;
  each normalized to `RequestCpus=2`, `RequestMemory=6000`.  No fit or
  downstream reconstruction consumer has been submitted; wait for all four
  sidecars to finish before the mandatory full output validation.
- 2026-07-30: All four calibration-mass sidecars completed with empty Condor
  stderr.  The formal full sidecar validator re-opened each ROOT, required its
  metadata and all 12 mass histograms, checked the 60--120-GeV/120-bin axis,
  and confirmed positive entries and `unreadable_outputs=0` for data 94/94
  and unmodified DY 75/75.  Its certificate is
  `calibration_validation/calibration_mass_sidecars_20260730T091240Z_479060.json`.
  The old schema-v1 human-readable JSON formatted large histogram entry counts
  to six significant digits; ROOT is authoritative and all physical histogram
  contents were verified.  The sidecar writer now records future counts as
  exact integers; no completed ROOT or physics input needed reprocessing.
- 2026-07-30: Per the user's approval, upgraded the portable AN-like Z fit to
  float both double-Crystal-Ball tail exponents and to initialize at the
  observed Z peak.  A 24-bin diagnostic used the certified four sidecars and
  showed the improved high-pT behavior; the user approved proceeding.  The
  formal `BW ⊗ double-CB +` exponential-background fit then atomically
  published the schema-validated central calibration payload
  `calibration/lepton_calibration_20260730T092122Z_data9p4_dy75_free_tail.json`
  and matching fit summary.  It records exactly the 94-file/9.4-fb data and
  75-file unmodified-DY sidecars, has 12 valid bins for each flavor, and no
  failed payload bin.  Its 24 per-bin fit PDFs/PNGs and four scale/smear PDFs
  are at `calibration/plots_20260730T092122Z_data9p4_dy75_free_tail/`.
- 2026-07-30: Re-reviewed the AN/paper and current rules before reconstruction.
  The AN-like 5/7-GeV object thresholds, 20/10-GeV candidate requirement,
  muon-only NanoAOD FSR recovery/isolation subtraction, no electron FSR, Z1/Z2
  pairing, smart cut, and 70-GeV candidate threshold are retained.  Before
  staging, the reconstruction driver was tightened to freeze the validated
  calibration payload under personal EOS while retaining its stable group-EOS
  provenance in every output; the selection configuration and future cutflow
  label use 9.4/fb.  The C++ build and Python syntax checks passed.  The exact
  minimal-PATH generic worker smoke processed a real pseudo-data NanoAOD shard,
  reopened all formal ROOT objects, checked summary/provenance, and verified
  idempotent task-summary routing.
- 2026-07-30: Staged and dry-run validated the full pseudo-data nominal
  reconstruction plan at
  `.../runs/20260730T135527Z_3487963_h4l_reco_pseudo-data/`: 94 distinct
  certified 9.4/fb inputs, 94 unique group-EOS outputs, personal-EOS frozen
  program/config/calibration inputs and logs, and a 6-GB one-worker request.
  Submitted it as cluster `972097`.  Initial status is 92 running, 2 idle, and
  0 held; eossubmit normalized every task to `RequestCpus=2`,
  `RequestMemory=6000`.  No downstream data consumer has been submitted.
- 2026-07-30: During the pseudo-data reconstruction queue, two accidental
  submissions of the same 94-task frozen plan overlapped.  Their formal output
  paths are identical and workers preserve an already published ROOT.  The
  generic runner now makes a post-publication auxiliary attempt-record write
  non-fatal, so a record-directory race cannot turn a valid retained or
  atomically published ROOT into a held task.  Per the user's current rule,
  do not stop or duplicate any completed/running task: after the live queues
  finish, fully re-open every formal ROOT and stage only scan-proven missing
  or invalid shards.
- 2026-07-30: The final pseudo-data reconstruction scan passed 94/94 current
  versioned outputs: zero missing/invalid ROOTs, 51,131,174 input events and
  309 selected candidates.  The report is
  `h4l_reco/validation/reco_output_snapshot_20260730T143808Z.json`.  The
  obsolete held tasks and redundant remaining duplicate queue work were
  removed only after this proof.  Before MC staging, a version audit found
  402--451 historical ROOTs per variation in the old layout, all with the old
  unversioned `calibration/lepton_calibration.json`; none can enter the
  current 9.4/fb calibration analysis.  The current
  `calib_fc4303851caa1060` target has no existing output, so all 404 official
  MC inputs are correctly pending in each of nominal and eight scale/smear
  variations.
- 2026-07-30: The pre-submit MC minimal-PATH worker smoke caught an MC-only
  bookkeeping error before any formal MC task was submitted: the output check
  compared the all-event `GenFiducial` tree to the subset of reconstructed
  candidates that are also fiducial.  The reconstruction now writes a distinct
  `Metadata:nGenFiducialEvents`, retains
  `nSelectedRecoGenFiducial` for the subset, and checks the tree against the
  former.  The initial temporary artifact is retained.  After rebuild, the
  same real unmodified GGH NanoAOD smoke passed with 603 reconstructed
  candidates and 1,027 particle-level fiducial events.  The full validator
  likewise now requires this MC consistency field.  The nine earlier staged
  MC plans froze the pre-fix binary and remain unsubmitted historical plans;
  fresh plans are required before any MC Condor submission.
- 2026-07-30: Re-froze the corrected reconstruction binary (SHA256
  `0f55db4535ea12b4...`) for all 404 unmodified official MC files and all nine
  calibration versions.  Every JDL passed a fresh 404-classad Condor dry-run,
  uses only personal-EOS frozen worker inputs/logs and the current versioned
  `fc4303851caa1060` output subtree.  Submitted the queues as: nominal
  `972099`, mu-scale up/down `972100`/`972101`, electron-scale up/down
  `972102`/`972103`, muon-smear up/down `972104`/`972105`, and
  electron-smear up/down `972106`/`972107`.  The initial status is 404 idle,
  zero held in every cluster.  Stop here for execution: once the jobs have
  run, fully scan formal ROOT outputs before any validation plot/template or
  fit.  MC may retain a recorded below-10% shortfall, but the accepted nominal
  files must have a valid matching ROOT in every one of the eight variations.

- 2026-07-30: The loose-v2 reducible-controls cluster `974844` completed 92
  queued shards plus its prior certified smoke; the only missing formal shard
  is `fake_data_10fb_shard_00011`.  Its original worker encountered a real
  1.53-TeV muon and aborted because the sentinel `65--999999` calibration
  tail had been treated as bounded.  The formal output scan finds exactly
  93/94 ROOTs, so all successful siblings are retained.  The terminal-bin
  semantics were corrected consistently in the reducible scanner,
  reconstruction, and calibration-skim consumers.  An exact worker-like
  minimal-PATH smoke with shard 00011 and the new scanner processed all
  545,750 events (3,703 fake-rate probes), atomically published/reopened all
  six required ROOT objects, and therefore exercises the former failure.
  A single-task recovery stager was added; it keeps the original immutable
  task/output identity, stages the corrected binary, dry-runs its JDL before
  submission, and refuses to overwrite a formal artifact.  Next: submit only
  this scan-proven missing data shard, then require the full 94/94 collector,
  merged fake-rate ROOT/payload, numerical totals, and visual inspection of
  the full low/high-pT rate maps before template work.

- 2026-07-30: Reinforced the permanent completed-stage gate: ROOT/JSON
  validation alone is insufficient.  Every completed/recovered stage must
  numerically compare plotted totals with certified inputs and visually inspect
  representative rendered figures; absent, empty, or inconsistent plots stop
  downstream progression.

- 2026-07-30: After syntax checking the recovery stager, confirming its
  immutable current scanner/runner hashes, and inspecting its one-row JDL,
  the recovery dry-run verified EL9, 1 CPU, 6000 MB, personal-EOS logs, the
  current scanner binary, and exactly the missing shard-00011 formal target.
  Submitted only that recovery as cluster `974848`; no completed sibling was
  resubmitted.  The next completion gate is formal 94/94 controls validation,
  followed by full fake-rate merge and numerical/visual rate-map inspection.

- 2026-07-30: The recovery completed normally (empty stdout/stderr) and the
  original controls plan now passes the formal 94/94 collector.  Its v2
  merge/payload and visual rate maps were completed and numerically checked:
  the merged numerator/denominator maps exactly equal all 94 certified inputs
  (muon 13,670/137,498; electron 3,316/202,124).  The inspected map directory
  is `reducible_background/validation/controls_plots_loose_v2_20260730T2205Z/`.
  The four gray 3--5-GeV electron bins are the intentionally unreachable
  region below the versioned 5-GeV electron object threshold, not a plotting
  defect; every reachable bin is measured.  The v2 products are preserved as
  historical after an additional AN/paper audit found a central-value defect:
  the fake-rate scanner had counted every extra loose lepton as a probe and
  omitted the AN's explicit 20/10-GeV tag and opposite-sign tag-probe
  mass-above-4-GeV conditions.  Paper Sec. 7.2.1 and AN Sec. 7.2.1.1 require
  exactly one loose probe.  Under the user's standing approval for direct
  AN-correctness changes, this was fixed and frozen as
  `reducible_2017_loose_v3.json`; the scanner rebuild and shell/config syntax
  checks pass.  Fresh v3 controls, merge, and payload are mandatory before
  template submission; no templates have been submitted.

- 2026-07-30: The v3 controls plan was staged with all 94 current 9.4/fb data
  shards.  Its exact minimal-PATH wrapper smoke reopened all six formal ROOT
  objects and processed 544,897 events; the equality of its 3,048 selected
  tags, 3,048 probes, and merged denominator entries explicitly tests the
  exactly-one-probe contract.  The 94-way JDL passed dry-run/routing checks
  (EL9, 1 CPU/6000 MB requested, current frozen scanner, personal-EOS logs,
  no modified NanoAOD), then was submitted as cluster `974858`.  Next gate:
  wait for the full data queue, remove only a diagnosed obsolete held task if
  one occurs, validate 94/94 ROOTs, merge the v3 payload, and inspect/numerically
  check the v3 fake-rate plots before staging data/prompt templates.

- 2026-07-30: Immediate `974858` status is 94 idle and zero held; eossubmit
  normalized the one-CPU/6000-MB request to 2 CPUs/6000 MB.  The dry-run
  audit procedure now records that inherited Condor classad values appear once
  rather than once per process, avoiding a false pre-submit failure.

- 2026-07-30: The AN-correct v3 controls queue completed with 94/94 valid
  ROOT/provenance products.  The v3 merge and rate-map check are complete at
  `reducible_background/validation/controls_plots_loose_v3_20260730T2236Z/`;
  the merged maps exactly equal all 94 certified inputs (muon 11,831/113,323;
  electron 2,756/166,744), and the visible maps are non-empty and physical.
  Before any template submission, a real minimal-PATH data-template smoke
  passed (87 control candidates and all required templates).  The required
  prompt-ZZ smoke exposed an MC scale-down edge case before Condor: an
  electron accepted at the 5-GeV raw loose threshold could shift infinitesimally
  into the deliberately unmeasured 3--5-GeV electron fake-rate bin.  The
  correction clamps fake-rate lookup to the flavor's raw loose threshold;
  this preserves the physical selected-object domain and affects only this
  boundary systematic evaluation.  The scanner has rebuilt; the unsubmitted
  pre-fix template/prompt plans remain historical and must be freshly staged,
  then both real worker-like smokes rerun before submission.

- 2026-07-30: Fresh plans froze the boundary-safe scanner and the certified v3
  fake-rate payload.  Their exact minimal-PATH smokes both pass: the data
  template smoke has 87 Z+X candidates and all 24 required ROOT objects; the
  unmodified public `ZZTo4L` prompt smoke has 5,335 3P1F candidates, a valid
  `genEventSumw=874272.3248`, and all 24 prompt-template objects.  The two
  JDLs pass count/routing/classad checks (94 complete pseudo-data tasks and
  exactly 188 allowed public prompt-ZZ/ggZZ tasks).  Submitted data templates
  as `974860` and prompt ZZ subtraction templates as `974861`.  Both initially
  have only idle jobs, zero held; eossubmit uses 2 CPUs/6000 MB for each.
  Next gate after both queues finish: fully validate data 94/94 and prompt
  MC coverage/per-sample ROOTs, merge the normalized prompt subtraction with
  the data template, then numerically and visually validate every final Z+X
  template figure before signal fitting.

- 2026-07-30: The data-template queue completed; the prompt queue produced
  166 formal outputs and 22 held jobs.  All 22 failures have the same
  non-physics EOS directory-creation race: concurrent first writes to a
  sample directory returned `EEXIST` even though the directory was usable.
  The runner now accepts a directory only after explicitly verifying both
  output and summary paths are directories, and treats this EOS `EEXIST` as
  idempotent.  A real public-MC smoke wrote two distinct prompt templates into
  the same already-created sample directory, reopening all 24 keys each time.
  The old cluster then contained exactly its 22 diagnosed held jobs (no running
  work), so it was removed after log capture.  A 22-task recovery queue keeps
  the old immutable task/output paths and reuses all 166 existing roots; its
  JDL passed dry-run/routing checks and was submitted as `975125`.  Next:
  wait for this recovery, then validate 94 data plus 188 prompt outputs and
  make/inspect the final Z+X templates.

- 2026-07-30: Recovery cluster `975125` is initially 22 idle, zero held,
  normalized by eossubmit to 2 CPUs/6000 MB.

- 2026-07-30: Prompt recovery completed without a held residual.  The formal
  collectors reopen/provenance-check all template products: data 94/94 and
  prompt MC 188/188 (ZZTo4L 109, GGZZ2E2Mu 22, GGZZ4Mu 31, GGZZ4E 26).
  The normalized AN Eq. 7.1 merge completed with 9.4/fb prompt normalization;
  the final nominal Z+X integral is 17.3331, exactly equal to the retained
  data `3P1F-2P2F` integral 18.1899 minus the audit prompt-ZZ subtraction
  0.8569.  All rate/template/final-state/audit figures were created and
  visually inspected at
  `reducible_background/validation/templates_plots_loose_v3_20260730T2214Z/`;
  they are nonempty and show the expected sparse-data bin fluctuations rather
  than a plotting failure.  Reducible background is ready for the final
  signal/background template and likelihood-fit stage.

- 2026-07-31: Re-read the AN/paper, the current workflow rules, and the
  latest certified reconstruction/reducible records before starting final
  inference.  The only permitted final inputs are the 94-file/9.4-fb
  pseudo-data reconstruction at `calib_8a0bbf9902cfa9c2`, the matching 404
  unmodified-public-MC files in nominal plus every scale/smear variation, and
  the AN-style loose-v3 data-driven Z+X merge.  A fresh non-overwriting final
  template configuration was frozen as
  `inference/config/template_config_loose_v3_20260731T005409Z.json`.
  Template building will enforce identical nominal/variation inventories,
  record the supplied zero cross-section uncertainty rather than invent one,
  atomically publish its NPZ/JSON products, and produce pre-fit m4l figures
  before any likelihood fit.

- 2026-07-31: The final deliverable was extended to include a compiled TeX
  analysis note in an orderly in-repository archive.  It will follow the
  reference AN/paper presentation, include all validated intermediate/final
  figures, give a dedicated principal-result comparison with the references,
  and document the implemented substitutions and unavailable selections or
  systematics.  This requirement is now a persistent project rule.

- 2026-07-31: The current final-template build completed and formally reopened
  its non-overwriting products at
  `results/templates_loose_v3_20260731T005409Z/`.  It contains all 94
  pseudo-data candidate files (481 candidates), all 404 nominal public-MC
  outputs and an exactly matching 404-file inventory in every eight lepton
  variation.  The 105--140-GeV pre-fit totals are data 67.0 versus 55.628
  expected (signal 8.532, qqZZ 27.195, ggZZ 0.875, data-driven Z+X 19.026).
  Every scale/smear and reducible variation has a finite positive applicable
  template.  The complete pre-fit figures were numerically checked and
  visually inspected at
  `results/plots_prefit_loose_v3_20260731T010147Z/`; they are nonempty in the
  inclusive and each final-state channel.  Two retained attempt-qualified JSON
  files document concurrent diagnostic invocations that correctly stopped at
  the no-overwrite gate; the verified formal NPZ/JSON are unaffected.  The
  next stage is the inclusive/final-state profile likelihood and simplified
  fixed-shape signal-translation mass fit.

- 2026-07-31: The inclusive, three-final-state, and simplified mass profile
  fits completed on the certified templates at
  `results/fit_loose_v3_20260731T010318Z/fit_h4l.json`.  The active correlated
  nuisance sources are muon/electron scale, muon/electron smearing, fake-rate,
  and reducible lepton-scale variations; the provided effective-cross-section
  uncertainty is explicitly zero.  The fit JSON was re-read for finite POIs
  and all 41 finite mass-profile points.  Final figures, including per-final-
  state m4l spectra, signal strengths, and the mass profile, were numerically
  checked and visually inspected at
  `results/plots_fit_loose_v3_20260731T010430Z/`.  The initial mplhep smoke
  exposed an unavailable cosmetic `mpl_magic` helper; it was removed before
  the formal figures were made, with no change to any template or fit value.
  Next: determine which particle-level fiducial/differential observables are
  actually represented in the certified NanoAOD-derived trees, then run only
  the sound optional result(s).

- 2026-07-31: During the required independent AN/paper-based final code audit,
  a central-value defect was confirmed in the final-inference reader: it
  clipped each signed prompt-subtracted Z+X bin to zero without preserving the
  AN Eq. 7.1 total.  The certified merge integral is 17.333066 (4mu 9.751805,
  4e 0.239266, 2e2mu 7.341995), whereas the preliminary NPZ used 19.025655
  after clipping.  The preliminary final fit/plots/fiducial products must not
  be treated as final.  Per the standing reprocessing rule, no replacement
  template or fit has been run: await the user's decision on a minimal
  non-negative smoothing/coarsening remedy that preserves every final-state
  and variation integral.  This requires only local template/fitting work,
  never NanoAOD or Condor reruns, and all existing artifacts are retained.

- 2026-07-31: The user approved the local final-inference correction.  A new
  versioned template configuration freezes the same certified 94-file/9.4-fb
  data, public-MC inventories, and v3 Z+X source, but replaces per-bin
  clipping with the narrowest non-negative reflected-Gaussian template that
  preserves each signed source integral exactly.  The replacement fit figures
  will use the final-state likelihood's own best-fit signal strengths and all
  nuisance parameters; fixed-nuisance one-sigma envelopes will be reported
  explicitly rather than inferred by subtraction of two profile intervals.
  The current fiducial helper is retained only as a labelled SM-extrapolated
  loose-v1 proxy and is excluded from the formal final results because it
  lacks the AN/paper response, efficiency, and nonfiducial-signal treatment.

- 2026-07-31: The corrected local template build completed with the same
  certified 94 data files and 404 nominal files in each of eight matching
  variations.  All 15 signed Z+X source templates were checked non-negative
  after smoothing and equal to their raw signed source integrals to
  1e-10.  The nominal Z+X yield is 17.333066 (9.751805 in 4mu, 0.239266 in
  4e, and 7.341995 in 2e2mu), replacing the invalid historical clipped yield
  19.025655.  The checked template, smoothing figures, and pre-fit figures
  are at results/templates_loose_v3_smoothed_20260731T030108Z/,
  results/plots_reducible_smoothing_loose_v3_20260731T030108Z_v3/, and
  results/plots_prefit_loose_v3_smoothed_20260731T030108Z/.

- 2026-07-31: The corrected profile fits, final-state post-fit plots, and
  limited reference comparison completed locally.  The fit uses its matched
  final-state best-fit vector for each post-fit spectrum, and all final
  figures were numerically checked and visually inspected.  The formal
  results are in results/fit_loose_v3_smoothed_20260731T030108Z/ and
  results/plots_fit_loose_v3_smoothed_20260731T030108Z/.  The simplified
  mass result is 128.180 GeV with a direct profile interval
  +1.320/-1.413 GeV and an explicitly recorded lepton scale/smear envelope
  +0.221/-0.058 GeV.  The figure archive and compiled six-page analysis note
  are in deliverables/final_reviewed_20260731T030108Z/; its archive validation
  confirmed all referenced figures and excluded historical preliminary
  010759Z/010846Z inference artifacts.

- 2026-07-31: A detailed, self-contained successor AN is being prepared before
  the requested central-value audit.  It will retain the existing reviewed
  note, include all final-version diagnostic figures (including all 40
  calibration-bin fits), document the full calibration/selection/background/
  likelihood workflow and MC catalogue, and describe the data solely as
  9.4 fb^-1.  No physics result or input is changed by this documentation work.

## 2026-09-23: UL16 PFNano v3 renewal (prompt_renew.md)

The task is now the renewal described in `prompt_renew.md`: the AN-16-442 /
JHEP 11 (2017) 047 method on the UL16 PFNano v3 benchmark pseudo-data
(20 fb^-1) and MC.  Everything above this heading is the obsolete 2017
Open Data record; it is kept for provenance only.  `AGENTS.md` and `README.md`
are untracked and are not edited until their dated copies are preserved in a
user-approved location (prompt §2.4, §10 step 3); until then new durable user
requirements are logged here.

- 2026-09-23T20:27Z, user requirement (durable; move to `AGENTS.md` after the
  preservation step): do not be constrained by the existing code.  The
  priorities are correctness, completeness and consistency with the paper and
  the AN, not reuse.  The existing code is only a reference.
- 2026-09-23T21:01Z, user requirement (durable; move to `AGENTS.md` after the
  preservation step): read large ROOT files with C++ wherever possible,
  because it is faster; switch to Python for plotting and similar work.  This
  also applies to ad-hoc inspection of ROOT inputs.
- 2026-09-23T21:01Z, user: `pixi install` had not been run earlier; it has now
  been run.  Re-check the pixi environment before relying on it.
- 2026-09-23T21:19Z, user requirement (durable; move to `AGENTS.md` after the
  preservation step): whenever Condor is used, submit through eossubmit and
  set its parameters correctly; every output (scientific products and Condor
  artifacts alike) goes into a folder inside `production_v3/`.
- 2026-09-23T21:32Z, user requirements (durable; move to `AGENTS.md` after the
  preservation step):
  - If the pixi environment lacks a package, add it to `pixi.toml` and run
    `pixi install` without asking (class-level permission for pixi packages).
  - Record every user requirement promptly in this log or in `AGENTS.md`.
  - Record and update the log and the plan promptly, and review them
    regularly (at least before every stage, submission and report).
  - Keep the old (2017) code and the new (v3) code in clearly separate
    folders, so that the two are never mixed up.
- 2026-09-23T21:32Z, user requirement (durable; move to `AGENTS.md` after the
  preservation step): follow the paper as closely as possible; where a
  limitation of the inputs makes an identical implementation impossible, do
  not force it: a reasonable and accurate substitute is enough.  Whenever
  several methods are possible, ask the user promptly which one to use.
- 2026-09-23T21:35Z, pixi (user permission of 21:32Z): uncommented
  `pyhf = ">=0.7"` and `iminuit = ">=2.25"` in `pixi.toml` and ran
  `pixi install` with `TMPDIR` and `PIXI_CACHE_DIR` inside
  `production_v3/tmp/`.  Installed pyhf 0.7.6 and iminuit 2.33.0; a two-bin
  pyhf/minuit smoke fit (`production_v3/work_reading/pyhf_smoke.py`) returns
  the injected mu = 1.00 +- 0.88.  Purpose: build and check `MODEL.json`.

## 2026-09-23T21:39Z: section 10 step 1 -- reading summary, inventory, gap analysis

### 1. Sources read

- `prompt_renew.md`, all 13 sections (read-only task specification).
- Repository: `AGENTS.md`, `README.md`, this log (2017 record), `prompt.md`,
  `.gitignore`, `pixi.toml`, `.analysis_config`, `retrieval_log.md`; every
  stage directory with its code, configuration and READMEs
  (`production_framework/`, `calibration_skim/`, `h4l_reconstruction/`,
  `reducible_background/`, `inference/`, `fake_data_nanoaod/`, `ntuplizer/`);
  `deliverables/final_reviewed_20260731T030108Z/` and the superseded
  `deliverables/figures/`; `work/`, including the AN and paper text dumps.
- `ref/AN2016_442_v8.pdf` page by page (printed pages 1-15 through
  `pdftotext` after the PDF-image budget ran out, pages 16-93 as images) and
  `ref/JHEP11(2017)047.pdf` (through `pdftotext`, checked against
  `work/paper_text_20260730.txt`).
- `ref/Slopbench H4l v0 (Copy).html`: metric definitions, `RESULT.json`
  contract, closure history and lessons, model-validation probes.  Its render
  script is extracted to `production_v3/work_reading/scoreboard_render.js`.
  Its sealed truths are not used.
- `slopbench_code_fork`: `AGENTS.md`, `README.md`, `seed_v3.md`, `seed_v2.md`,
  `seed.md`, `injection/README.md`, the header of `injection/pfnano_inject.cpp`,
  `.claude/skills/slopbench-eval/SKILL.md`, `benchmark/production/PFNANO_V3.md`,
  `PFNANO_V2.md`, `PFNANO_EFFICIENCY.md`, `PFNANO_TNP.md`, the plans in
  `benchmark/configs/`, `benchmark/tnp/*.py`,
  `results/h4l_seeds/v3/an/v3_validation.pdf`, the YR4 spreadsheet and the
  legacy `contract/`.  Of `h4l_pfnano_seed_plan_v3.json` only the sample,
  cross-section, selection, branch-policy and trigger blocks are used.  The
  response and efficiency definitions are read only to understand the
  benchmark; they are never analysis inputs.
- `slopbench_code_main`: `eval/reference_hzz4l.json`,
  `eval/model_emits_block.md`, `eval/cross_sections_effective.json`,
  `eval/run_eval.py`, `eval/harness_analysis.py`, `eval/gof_combine.py`,
  `eval/truth_from_record.py`, `eval/gen_truth.py`, `harness/judge.py`.
- jfc framework: `agents/`, `conventions/`, `methodology/`, jfc `README.md`.
- Not opened (section 4.3): `injection/records/`, `h4l/ref_v3/injection_tmp/`,
  `results/h4l_seeds/` except `v3/an/` and the YR4 spreadsheet, any
  `eval/truth/`.

### 2. Verified facts about the inputs

Single-file inspections with C++ (and, before the 21:01Z rule, PyROOT) in
`production_v3/work_reading/`.  No data entry count, shard size or provenance
field other than schema and luminosity was looked at.

Pseudo-data (`production_v3/pseudo_data/917f9326efa44ea9/`, 83 shards):
- keys `Events`, `Runs` (0 entries), `LuminosityBlocks` (0 entries),
  `PFnanoFakeDataProvenance` (schema `pfnano_fake_data/v3`, `lumi_fb` 20);
- 1,280 branches; no `run`/`luminosityBlock`/`event`, generator, `Pileup_*`,
  `genWeight` or PF-candidate branches; every data branch exists in the MC;
- every HLT path of AN Table 2 exists in data and MC (NanoAOD names);
- electrons: `mvaFall17V2noIso`/`mvaFall17V2Iso` raw scores and
  `_WPL/_WP90/_WP80` flags, `cutBased`, `cutBased_HEEP`, `mvaTTH`, `lostHits`,
  `convVeto`, `deltaEtaSC`, `energyErr`, `pfRelIso03_all/_chg` (rho x A_eff
  corrected, no components), `sip3d`, `dxy`, `dz`;
- muons: `looseId` (= PF and (global or tracker)), `mediumId`, `tightId`,
  `softId`, `highPtId`, `isPFcand`, `isGlobal`, `isTracker`, `isStandalone`,
  `nStations`, `nTrackerLayers`, `ptErr`, `pfRelIso03_all/_chg`,
  `pfRelIso04_all`, `sip3d`, `dxy`, `dz`, `fsrPhotonIdx`, `tunepRelPt`;
- `FsrPhoton_{pt,eta,phi,relIso03,dROverEt2,muonIdx}`: muons only; `Photon`
  starts at about 10 GeV;
- `IsoTrack` excludes lepton tracks (it matches 0.4-1.3 % of prompt MC
  muons), so it cannot probe the muon reconstruction efficiency;
- `TrigObj_{id,filterBits,pt,eta,phi,l1pt,l1iso,l2pt}`,
  `PV_{npvs,npvsGood,ndof,x,y,z}`, `MET_*` (re-clustered), `RawMET_*`,
  `ChsMET_*`, `Flag_*`, `fixedGridRho*`;
- jets are re-clustered from PF candidates identically in data and MC:
  `Jet_jetId` 2 (tight) or 6 (tight + lepton veto), no loose bit in UL16;
  `Jet_puId`; `Jet_btagDeepFlavB` defined for 99.7 % of jets with pT > 30 GeV,
  while `Jet_btagCSVV2`/`Jet_btagDeepB` are -1 for 17-35 % of them;
  `Jet_rawFactor` median 0.048 (energy corrections applied); MC reco/gen
  response median 1.00-1.05.

MC (`h4l_seeds_v2/mc/`): DYJetsToLL 3,637 files, 190 GB, 4 blocks; ZZTo4L
1,854, 124 GB, 5 blocks; TTBar 488, 58 GB, 6 blocks; GGZZ4Mu 48, 6.5 GB;
GGZZ4E 3, 0.48 GB; GGZZ2E2Mu 1, 0.04 GB; GluGluToHToZZ_M125 99, 0.85 GB;
VBF_HToZZ_M125 98, 0.90 GB; VHToZZ_M125 99, 0.88 GB (flat).  Total 6,327 files,
381 GB.
- `Runs`: `genEventCount`, `genEventSumw`, `genEventSumw2` and LHE sums;
  `PFnanoFinalMCProvenance.entries_before_selection` equals `genEventCount` in
  every inspected file.
- `ZZTo4L` carries `LHEScaleWeight`, `LHEPdfWeight`, `LHEReweightingWeight`,
  `PSWeight` and `LHEPart`.  The signal carries `PSWeight` only, plus
  `HTXS_stage_0`, `HTXS_Higgs_{pt,y}`, `HTXS_njets30`, `GenPart` with
  `statusFlags`, `GenDressedLepton` (NanoAOD dR 0.1 dressing), `GenCands`
  (stable particles), `GenJet` (without neutrinos), `Pileup_*`.
- All nine effective cross sections equal the prompt table and
  `samples[].effective_xsec_pb` of `h4l_pfnano_seed_plan_v3.json` (TTBar
  87.58 pb).  YR4 spreadsheet sha256 `47b55325...809ce49c` verified.
- Count- vs genWeight-based normalization, every 37th file (ZZTo4L 50 files,
  TTBar 14 files): ratio of the preselected fractions (weight/count) 1.0008
  and 1.0010; negative weights 0.46 % and 0.36 % of the selected events; Kish
  n_eff/n 0.982 and 0.986.  To be recomputed on the full manifests.
- Lepton ID efficiencies, gen-matched prompt leptons that pass the AN loose
  cuts (ZZTo4L, 6 files; ggH agrees within 1-2 %):
  - muons, 5-15 / > 15 GeV: `looseId` 0.998 / 0.998, `mediumId` 0.982 / 0.991,
    `tightId` 0.958 / 0.974; presence in the collection 0.997 / 0.999;
  - electrons, 7-15 / > 15 GeV: noIso WPL 0.912 / 0.989, noIso WP90
    0.812 / 0.865, noIso WP80 0.616 / 0.728, Iso WPL 0.925 / 0.989, cutBased
    >= veto 0.756 / 0.919; presence 0.917 / 0.964; a `Photon` (supercluster)
    exists for 0.46 / 0.98 of the prompt electrons.
- Conditions: RunIISummer20UL16, post-VFP (`seed_v3.md`, `PFNANO_V2.md`).

Tools: pixi (Python 3.14.6, numpy 2.5.1, scipy 1.18.0, uproot 5.7.5, awkward
2.11.0, hist 2.10.1, boost-histogram 1.7.2, matplotlib 3.11.1, mplhep 1.3.2,
vector, particle, rich, pyhf 0.7.6, iminuit 2.33.0, pandoc 3.9 +
pandoc-crossref; no ROOT); system ROOT 6.40.04 with RooFit, RooStats,
HistFactory, RooFitHS3, Minuit2 (OpenMP), TMVA, MathMore and PyROOT for system
python3 3.9; g++ 11, gfortran 11.5, cmake, make, pdflatex (TeX Live 2020),
latexmk, apptainer, the combine container on CVMFS, `lxbatch/eossubmit`.  MELA,
JHUGen and MCFM are in neither pixi nor the CVMFS CMSSW externals.  Proxy:
596,873 s left on 2026-09-23.

### 3. The reference method, condensed (AN-16-442 v8, HIG-16-041)

- Triggers: OR of the 2016 paths of AN Table 2; efficiency from the 4l
  tag-and-probe method of AN 2.1.2 (tag matched to a single-lepton trigger
  object; the three probes must rebuild any analysis trigger).
- Leptons (AN 3): electrons pT > 7, |eta| < 2.5, dxy < 0.5, dz < 1; tight =
  BDT (Table 7) and RelPFIso(0.3, rho A_eff) < 0.35.  Muons pT > 5,
  |eta| < 2.4, dxy < 0.5, dz < 1, global or arbitrated tracker, ghost
  cleaning; tight = PF (pT < 200) or PF/tracker-high-pT (pT > 200) and
  RelPFIso(0.3, delta-beta) < 0.35.  SIP3D < 4; electrons within dR < 0.05 of
  a selected muon are removed.
- FSR (AN 3.3): PF photons pT > 2, |eta| < 2.4, relIso < 1.8, supercluster
  veto, attached to the closest loose+SIP lepton, dR < 0.5, dR/ET^2 < 0.012,
  the lowest dR/ET^2 per lepton, removed from the isolation sums.
- Jets (AN 3.4): AK4 PF, loose ID, pT > 30, |eta| < 4.7, dR > 0.4 from tight
  leptons and FSR photons; b tag CSVv2M > 0.8484.
- ZZ candidates (AN 4): 12 < m(ll(gamma)) < 120, Z1 closest to m_Z,
  m_Z1 > 40, pT 20/10, dR > 0.02, every OS pair m > 4 (no FSR), smart cut,
  m4l > 70; best candidate by the highest D_bkg^kin (for the same four
  leptons, Z1 closest to m_Z).
- Discriminants (AN 5): MELA D_bkg^kin, D_2jet, D_1jet, D_WH, D_ZH; per-event
  mass uncertainty with lambda corrections and closure (AN 5.3); Z1-mass
  constrained refit (Eq. 16).
- Categories (AN 6): seven, defined with the MELA-only working points
  (D_2jet > 1.043 - 460/(m4l + 634), D_WH > 0.951 or D_ZH > 0.9937,
  D_1jet > 0.697), jets, b tags, additional leptons and MET.
- Backgrounds (AN 7): qqZZ with NNLO/NLO and EW k-factors; ggZZ with a
  k-factor and 10 %; Z+X from OS (2P2F, 3P1F with prompt subtraction) and SS
  (conversion correction, r_OS/SS) methods, combined; Landau shapes.
- Signal (AN 8): DCB with parameters linear in m_H per final state;
  non-resonant WH/ZH/ttH as a Landau.
- Statistics: unbinned 2D L(m4l) L(D_bkg^kin | m4l) in 21 subsamples; mass
  fits 1D/2D/3D with and without the Z1 refit; on-shell width with
  interference; STXS stage 0; fiducial volume of paper Table 4 (dressing
  dR < 0.3, isolation from stable particles except e, mu, nu, jets pT > 30
  and |eta| < 2.5, Z1 closest to m_Z then the highest-sum-pT Z2 at generator
  and reconstruction level, fit at fixed m_H without categories, final-state
  fractions free); systematics of AN Tables 10, 11, 21, 22.

### 4. Scoring requirements, condensed

- Scoreboard: recovery (exit 0, valid `RESULT.json`); per POI the pull,
  sensitivity gate |slope - 1| <= max(0.15, 3 err), bias
  |mean pull| <= max(3/sqrt N, 0.05), width |w - 1| <= max(0.25, 3/sqrt(2(N-1))),
  precision median min(1, sigma_ref/sigma); S_phys = recovery x mean over
  POIs (sens x calib x precision); scored if the expected Z >= 3; reference
  closure also needs max|pull| <= 3.5; `MODEL.json` MLE within 0.25 sigma + 0.05.
- v3 eval: procedural = mean(systematic groups covered/12, N-1 cut ablation
  on the harness's own cuts, combine saturated GoF with 0.05 < p < 0.95);
  accuracy = pulls of `mu_incl`, `muScaleShift`, `elScaleShift`, `resSmear`
  (efficiency alongside).  The judge reads only .json .txt .tex .md .csv .log
  .yml .yaml files, up to 20,000 characters each.

### 5. Scoring-code issues (to report; the analysis is not tailored to them)

1. `harness_analysis.py` and `cross_sections_effective.json` still use the old
   five signal names and TTBar 52.70 pb; the v3 directories `*_M125` and
   `VHToZZ_M125` are not found by its `mc/<proc>` lookup.
2. `gof_combine.py` applies the per-lepton smear directly to m4l; for four
   similar legs the m4l smear is about r/2.
3. The efficiency truth is one tier for both flavours (default
   Muon_mediumId / Electron WP90), whatever working points the analysis uses.
4. The `model_emits` `gof` is read from the self-report, whose schema has no
   `gof` field.
5. The harness takes the reducible background from DY/TT MC, and its muon
   isolation uses the 0.4 cone (the AN uses 0.3 for both flavours).
6. TTBar 52.70 pb also survives in `contract/TASK_CARD.md`, the fork's
   `seed.md`, `gen_truth.py` and the old `prompt.md`.

### 6. Component inventory

Verdicts: *rewrite* = new implementation in the v3 code tree, the old code
being only a reference (user rule of 20:27Z); *reference* = read for its
patterns, not executed; *discard* = not needed for v3; *keep* = left untouched.

| component | verdict | reason |
|---|---|---|
| `production_framework/scripts/build_manifest.py` | rewrite | 2017 catalogue, opendata prefix and mixing-certificate gates; no ROOT-level validation |
| `production_framework/scripts/stage_condor.py`, `monitor_plan.py` | rewrite | group-EOS output prefix, 2017 work root; keep its sha-frozen program inputs and JDL/log checks as patterns |
| `production_framework/condor/run_root_task.sh` | rewrite | exit 75 turns a data-input failure into a silent `skipped_unreadable_input` |
| `production_framework/scripts/eossubmit_env.sh` | reference | re-created in the v3 tree with the approved TMP location |
| `samples_2017.json` (three copies), `CMakeLists.txt`, `program_options.example.json` | discard | 2017 catalogue; Makefiles are used |
| `calibration_skim/calibration_skim.cpp` | rewrite | requires event IDs (silent skip), 2017 objects |
| `calibration_skim/calibration_mass_sidecar.cpp` | discard | unweighted MC, no trigger, fixed binning |
| calibration coverage/validation programs and scripts | reference | integrity and coverage patterns |
| `calibration_skim/scripts/fit_lepton_calibration.py` | rewrite | SciPy fit, 5/7 GeV defaults, dilepton-level results |
| calibration bin files, payload schema, identity payload | reference | half-open bins, unbounded tail bin |
| `h4l_reconstruction/src/h4l_reconstruct.cpp` | rewrite | event IDs required and used for seeds, 2017 triggers and WPs, no discriminants, crude categories, fiducial proxy |
| `h4l_reco_output_validation.cpp`, `stage_reconstruction.py`, plot scripts | rewrite | new schema |
| `analysis_2017*.json`, `calibration_identity.json`, `tests/` | discard | 2017 configuration |
| `reducible_background/src/reducible_scanner.cpp` | rewrite | OS only, event IDs, 9.4 fb^-1 gate; its OS lessons are kept |
| `reducible_background/src/merge_reducible_background.cpp` | rewrite | no SS method, no combination |
| reducible Condor stagers and `reducible_2017*.json` | discard | 94-file / 9.4 fb^-1 gates |
| `inference/src/build_templates.py` | rewrite | keep the integral-preserving Z+X smoothing idea |
| `inference/src/fit_h4l.py` | discard | SciPy binned fit, shift-only mass, stat from a separate likelihood (stat > total), bounded mu |
| `inference/src/derive_fiducial_cross_section.py` | discard | non-formal proxy |
| other `inference/src/*.py` and `inference/config/*` | discard / rewrite | old results and conventions |
| `fake_data_nanoaod/` | keep | historical mixing code, not needed (data delivered); `src/atomic_output.h` is a reference |
| `ntuplizer/` | keep | tracked, unused |
| `deliverables/final_reviewed_20260731T030108Z/`, `deliverables/figures/` | keep / reference | only record of the old results; never modified |
| `work/` | reference | text dumps, checked against the PDFs |
| `AGENTS.md`, `README.md` | preserve, then update | section 12 |
| `.analysis_config` | update | points to deleted 2017 paths |
| `pixi.toml` | updated 21:35Z | pyhf and iminuit; the tectonic `build-pdf` task is obsolete |
| `.gitignore` | propose | stale entries, unignored binary, contradictory comment |
| `task.local.*.json` (repository root) | keep | old smoke-test task files |

### 7. Gap analysis against sections 6 and 9

| requirement | old analysis | needed |
|---|---|---|
| 6.1.1 yields and distributions | m4l per final state only | full/low-mass m4l per final state and category, m_Z1/m_Z2 (2D in 118-130), D_bkg^kin vs m4l with D_mass, category discriminants, Tables 1 and 2 |
| 6.1.2 signal strength | SciPy fit, stat > total, mu_4mu on the boundary | observed/expected mu with fixed-nuisance stat/syst; per final state, category and production mode; 2D contours |
| 6.1.3 STXS stage 0 | none | HTXS stage 0, abs(y_H) < 2.5, normalized to the SM |
| 6.1.4 fiducial | non-formal proxy | Table 4 volume, A_fid, epsilon, f_nonfid, per final state, differential pT(H), N(jets), pT(jet1) with response-matrix unfolding |
| 6.1.5 m_H | shift-only template | 1D/2D/3D with and without the Z1 refit, per final state, stat/syst, pre/post-fit expected, D_mass calibration and validation, Z->4l m_Z |
| 6.1.6 width | none | on-shell Gamma_H limit, observed and expected |
| 6.1.7 backgrounds, systematics | OS Z+X only, partial systematics | OS + SS + combination; full systematics table with impacts |
| 6.1.8 paper comparison | old figure | new figure with the new results |
| 6.2 m_H, mu with the SM expectation at m_H | fixed-mass normalization | YR4 sigma x BR(m_H) and A x eps(m_H) with systematics |
| 6.2 calibration | Z-peak scale/smear (unweighted MC, no trigger); no tag-and-probe | per-lepton scale and smear at (45 GeV, 1.2) in every bin; tag-and-probe SFs; eff_correction |
| 6.2 other scored items | none | significance and discovery flag, saturated GoF, declared cuts with N-1 Z, systematic coverage, pyhf workspace |
| 9.1 contract and closure | none | `RESULT.json` and generality (no dataset-specific code), MC pseudo-experiments with slopes, pulls, widths, false discoveries, paired injection, crosstalk, GoF uniformity, toy bias, MODEL consistency |
| 9.2 v3 eval | none | `eval_selfreport.json`, harness-recognizable cut names, judge-readable text files |
| 9.3 conventions | smear at dilepton level | per-lepton scale/smear/SF at (45 GeV, 1.2), flavour combination documented, per-flavour extras |

Conflict recorded: section 2.6 of the prompt puts new files where their
function already lives, while the user's newer rule (21:32Z) asks to separate
old and new code folders.  The newer instruction wins; the concrete layout is
question Q1 of `PLAN.md`.

## 2026-09-23T21:44Z: section 10 step 2 -- plan v1 and questions

- The complete plan (stages, commands, outputs, validation, Condor resources,
  ordering, volume) is the new living file `PLAN.md`, version 1
  (sha256 23f8207cc16e9c5a...).  The plan is kept there and updated when it changes; this
  log records each version and decision.
- Questions Q1-Q22 (PLAN.md section 6) cover every decision point of prompt
  section 11 plus: the code-tree separation, the submission TMP location, the
  preservation location, the `.gitignore` changes, the muon reconstruction
  efficiency (no usable probe object), the trigger-efficiency method, where
  the fits run, the fiducial isolation/dressing ambiguity and the fixed or
  profiled m_H of the fiducial fit, and the stage-2 skim design.
- Asked the user on 2026-09-23T21:44Z.  No code, configuration or Condor
  change before the answers and the approval (the only change so far is the
  user-authorized pixi addition of 21:35Z).
- 2026-09-23T21:48Z, user requirement (durable; supersedes the approval
  gates of prompt section 2.2 and the ask-before-install rule): ask the user
  questions directly (interactively).  Installing what is needed, changing
  code and configuration, and submitting Condor jobs may be done directly,
  without asking.  The remaining limits: change only files inside this
  repository directory, and never delete files (ask if a deletion is really
  needed).  The 21:32Z rule still holds: when several methods are possible,
  ask which one to use.
- 2026-09-23T21:51Z, answers, round 1:
  - Selection profile: do not fix it now.  Start from the AN values, but the
    data are smaller than the AN's, so the optimum may differ: optimize the
    thresholds for the actual situation.  (Plan: optimize on MC expected
    sensitivity at 20 fb^-1, never on the data.)
  - Electron ID: likewise, study the actual thresholds; training a BDT, NN or
    similar is allowed if needed.
  - Muon reconstruction efficiency: the user asks whether tag-and-probe can
    be used; if it cannot, explain why and ask.  (Follow-up asked in round 2.)
  - Mass morphing: rest-frame scaling of the 4l system is accepted, with a
    corresponding systematic uncertainty.
- 2026-09-23T21:55Z, answers, round 2:
  - Muon reconstruction efficiency (a standard tag-and-probe has no probe
    independent of the muon collection): measure it with all three methods,
    the Z yield ratio, a tag-and-probe with HLT trigger-object probes, and
    the MC value.  Choose the most suitable one in the end: if all agree,
    the one with the smallest statistical uncertainty.  The differences
    between the methods become a systematic uncertainty.
  - mu is unbounded, and the m_H range is wider (initial choice
    110-140 GeV, with the m4l window widened accordingly; final values from
    the toy study, with no boundary artifacts).  Discovery stays at a local
    significance of at least 5 sigma (default, not contradicted).
  - MC normalization: genWeight-based (sum of genWeight over the selected
    events / genEventSumw over the files present).  The count-based value
    remains a documented cross-check.
  - Trigger efficiency: the AN 4l method plus tag-and-probe, but check
    whether the data suffice.  If they do not, a Z-peak tag-and-probe per
    lepton in pT and eta bins may be used instead.
- 2026-09-23T21:59Z, answers, round 3 (all as recommended):
  - Fiducial isolation: the photons used for dressing (dR < 0.3) are not
    counted in the isolation sum (consistent with the reconstruction-level
    FSR removal).
  - Fiducial fit: m_H profiled; a variant at 125.09 GeV for the paper
    comparison.
  - MODEL.json: binned pyhf templates at the fitted m_H (final state x
    category channels, 1 GeV m4l bins, set-(b) nuisances as
    normsys/histosys + staterror), checked against the reported mu.
  - Result sets: the prompt's proposal (RESULT.json, MODEL.json and the
    `measured` block of eval_selfreport.json use set (b); the `systematics`
    block describes the full model (a); luminosity 2.5 % in set (a)).
- 2026-09-23T21:59Z, answers, round 4:
  - Optimization objective for the selection thresholds and the electron
    ID (on MC only, never on the data): the joint expected precision of mu
    and m_H, i.e. minimize (sigma_mu/sigma_mu,AN)^2 + (sigma_mH/sigma_mH,AN)^2
    from Asimov data (mu = 1, m_H = 125 GeV, 20 fb^-1).
  - Unblinding: after the freeze, ask the user where the truth comes from;
    the user will say.
  - Flavour combination in the RESULT.json calibration block: keep both
    (the muon and the electron values) for now; the user will revise
    RESULT.json later.  Ask again when RESULT.json is written.
  - muEffScale/elEffScale: the SF of the full single-lepton selection
    (reconstruction x ID x isolation x SIP) at (45 GeV, 1.2); the
    components as extras.
- 2026-09-23T21:59Z, decisions taken under the 21:48Z permission (no
  alternative physics method involved; each is reversible):
  - New code in a new top-level tree `analysis_v3/` (PLAN.md 2.1); the 2017
    folders stay untouched.  Products per PLAN.md 2.2; Condor artifacts in
    `production_v3/condor/<stage>/<UTC stamp>/`; submission-time TMP in
    `production_v3/tmp/condor_submit/` (to be verified by an eossubmit dry
    run before first use).
  - AGENTS.md/README.md dated copies in `deliverables/pre_v3_rules_20260923/`.
  - `.gitignore`: the four changes of PLAN.md Q4.
  - Electron FSR: muon FSR only (`FsrPhoton` with the AN criteria), since
    NanoAOD has no electron FSR photons and `Photon` starts at 10 GeV.
  - Kinematic discriminants: download and build JHUGen-MELA into
    `production_v3/external/` (the AN method); if it cannot be built, ask.
  - Jets: pT > 30, abs(eta) < 4.7, `jetId >= 2`, dR > 0.4 cleaning, DeepJet
    medium (UL16 post-VFP, 0.2489) as the only usable b tagger; AN
    categories.
  - Width: BW(Gamma_H) (x) per-event resolution without interference (the
    interference model is unavailable).
  - Fits in compiled C++ RooFit/Minuit2; Python only for plots,
    orchestration and the pyhf export.
  - AN written directly in LaTeX and compiled with pdflatex; plot label
    "CMS Private work", data called pseudo-data.
  - Stage 2 as one pass writing the dilepton tables, the multilepton slim
    NanoAOD and the signal generator table.
  - Muon ID: the AN PF definition (`looseId`, `highPtId` above 200 GeV) as
    the baseline; muon WP and thresholds are part of the optimization study.
- 2026-09-23T22:05Z, log correction: the times of the entries appended in
  this session (21:01Z-21:59Z) had first been estimated, some later than the
  clock; they were corrected in place to the real append times taken from the
  session record.  The uncorrected file is kept as
  `production_v3/tmp/experiment_log.before_time_fix.md`.

## 2026-09-23T22:05Z: stage 0 done (preservation and rules)

- PLAN.md version 2 (sha256 2b89a573ebb3b030...): the decisions table replaces the
  questions; stages 2-4 and 7 updated for the user's answers (MC-only
  selection and electron-ID optimization as stage 4b, three-way muon
  reconstruction efficiency, genWeight normalization, wider m_H range).
- Preserved the 2017 rules: `deliverables/pre_v3_rules_20260923/AGENTS.md`
  (sha256 e8f41cae83cb98be...) and `README.md` (sha256 d70e5e545819031c...),
  identical to the originals.
- `AGENTS.md` rewritten for v3 (sha256 b42bec08a9ca4a18...): prompt section 12 applied
  (loose-v1, 9.4 fb^-1/94-file, mixing, skip-registry and old-path rules
  removed; "do not ask" and unattended-monitoring rules replaced); every user
  requirement of 2026-09-23 added; the integrity, Condor, plotting, Z+X and
  fit-reporting rules kept.
- `README.md` rewritten as the v3 workflow skeleton (sha256 381b18f88c09f232...); the
  old notes stay in the preserved copy.
- `.analysis_config`: data_dir = the v3 pseudo-data, allow = the v3 MC.
- `.gitignore`: comment reworded; stale `calibration_skim/calibration_skim`
  and `calibration_skim/calibration_coverage` lines removed; new rules
  `h4l_reconstruction/h4l_reco_output_validation` and `/analysis_v3/*/bin/`.
  `git ls-files | git check-ignore --stdin --no-index --verbose` prints
  nothing, and sample source paths under `analysis_v3/` are not ignored.

## 2026-09-23T22:25Z: stage 1 (infrastructure) -- code, plan and local smoke

- New code (all under `analysis_v3/`):
  - `common/include/h4l/{io.h,hash.h,root_io.h}`: EOS-safe directories,
    atomic text publication, FNV-1a/splitmix seeds (file key, entry, object,
    variation), ROOT open-with-retries and close/reopen/tree validation;
    `common/python/h4l_style.py`: mplhep CMS style, "Private work" label,
    PDF+PNG output.
  - `framework/`: catalogue `config/samples_ul16_v3.json` and branch
    requirements `config/input_requirements_v3.json`; C++ `validate_inputs`
    and `root_check`; `condor/run_task.sh` (sandbox run, copy to an
    attempt-qualified EOS name, size and root_check verification, `mv -n`
    publication, JSON last, attempt records); `scripts/stage_condor.py`
    (content-frozen program inputs, run directory, JDL, local smoke with
    PATH=/usr/bin:/bin, dry-run check, submit and cluster-ID record);
    `scripts/build_manifests.py` (plan: the only directory listing;
    assemble: full output scan and manifests); `scripts/eossubmit_env.sh`
    (TMP in `production_v3/tmp/condor_submit/`, eossubmit pool, proxy
    check >= 86400 s); `scripts/plot_manifests.py`.
- Plan `production_v3/manifests/v1/plan.json`: 84 tasks (pseudo-data
  83 shards in 17 tasks of 5; MC 6,327 files in 67 tasks of up to 100).
- Local minimal-PATH smoke through the frozen wrapper: validate_pseudo_data_0000
  (5 shards ok, 2 s) and validate_ZZTo4L_0000 (100 files ok, 58 s); both
  published their final JSON; attempt records written.  eossubmit
  environment verified (pool eossubmit, proxy 591,582 s left).
- Two independent reviews (production/Condor safety; physics bookkeeping)
  started before the first submission.
- 2026-09-23T22:35Z, finding for stage 3 (trigger efficiency): in the UL16
  NanoAOD (run2_HLTconditions_2016, CMSSW_10_6_30 triggerObjects_cff.py) the
  muon TrigObj quality bits are only 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter
  PFTau, 8 = IsoTkMu, 1024 = Mu50, and muon objects passing none of them are
  not stored.  Measured in one ZZTo4L file: muon objects carry bits 0, 1, 2,
  3, 10 only; electron objects carry the 2e, 3e, 1e-1mu, 2e-1mu and 1e-2mu
  leg bits.  Hence the muon legs of HLT_TripleMu_12_10_5,
  HLT_DiMu9_Ele9_CaloIdL_TrackIdL and HLT_Mu8_DiEle12_CaloIdL_TrackIdL cannot
  be rebuilt from trigger objects in the AN 2.1.2 4l method; its impact is to
  be quantified in stage 3 (user rule: fall back to the Z-peak per-leg
  tag-and-probe if the 4l method is insufficient).
- 2026-09-23T22:40Z, findings for stages 2-3 (CMSSW_10_6_30 NanoAOD sources
  on CVMFS):
  - FSR: `muonFSRProducer` keeps its default `muonPtMin = 20 GeV`, so
    `FsrPhoton` exists only around muons above 20 GeV; its preselection is
    looser than the AN (photon pT > 2, abs(eta) < 2.5 outside 1.4-1.6,
    1e-4 < dR < 0.5, closest muon, electron-footprint veto, relIso03 < 2.0,
    dR/ET^2 < 0.05).  The AN cuts (relIso < 1.8, dR/ET^2 < 0.012,
    abs(eta) < 2.4, closest selected lepton) are applied on top; FSR of
    muons below 20 GeV (and of electrons) cannot be recovered: a documented
    limitation.
  - Isolation with FSR removal is exact from the NanoAOD sums: iso =
    chg + max(0, (all - chg) - sum pT(FSR in cone)/pT), because
    `pfRelIso03_all` = (chg + max(0, neutral - PU))/pT.
  - Muon collection (NanoAODv9 default): pT > 15 GeV kept whatever the ID;
    3-15 GeV kept only if it passes CutBasedIdLoose, SoftCutBasedId,
    SoftMvaId or a high-pT ID.  A reconstruction inefficiency invisible to a
    collection-based probe can therefore occur only below 15 GeV; above
    15 GeV standalone-only muons can serve as reconstruction probes.
- 2026-09-23T22:47Z, stage 1 reviews and fixes.  Production/Condor review:
  two A findings, fixed (explicit `initialdir` with an Iwd assertion, since
  `/eos/user/y` resolves to `/eos/home-y`; cluster/process IDs passed with
  `environment`, and attempt names made unique with UTC time and random
  digits), plus B/C fixes (records with message and program-output tails,
  sha256 check of every EOS copy, orphan-ROOT handling, root_output parsing
  checked, JSON `publication` block with the ROOT sha256, program and config
  provenance in every report, full read of every required data branch,
  transient errors not published, restartable all-or-nothing assemble, an
  in-flight guard, a passed-smoke requirement for --submit, and Arguments,
  Iwd and Environment checks in the dry run).  Physics review: no A finding;
  its B/C points were adopted (hard genEventCount/entries_before_selection
  check, Runs LHE sums recorded, consumers re-sum genEventSumw over the files
  they process, collection-tagged seeds with one deviate per lepton,
  preselection bias of probe-object reconstruction tag-and-probe, record
  coverage from seed_v2.md, genPartIdx/genPartFlav and pfRelIso03_chg
  required, data shard sizes dropped from data.json).  The rules are in
  AGENTS.md.  The physics reviewer wrote a scratch file outside the
  repository (session scratchpad `.../scratchpad/report.txt`); it is left
  in place (no deletion) and was reported to the user.
- 2026-09-23T22:47Z, plan v2 (`production_v3/manifests/v2/plan.json`,
  84 tasks; v1 is superseded and unused).  Smoke through the frozen wrapper
  with PATH=/usr/bin:/bin: validate_pseudo_data_0000 (5 shards, full read,
  peak 345 MB) and validate_ZZTo4L_0000 (100 files, peak 381 MB).  eossubmit
  dry run ok (82 jobs; Iwd, Arguments, Environment, logs verified).
- 2026-09-23T22:47Z, SUBMITTED stage 1: cluster 1152509, 82 jobs,
  run directory production_v3/condor/manifests_v2/20260923T224646Z,
  1 CPU, 1000 MB, 1000 MB disk, espresso.  Completion gate: all 84 outputs
  in production_v3/manifests/v2/validation/, then `assemble`.
- 2026-09-23T22:56Z, STAGE 1 DONE.  Cluster 1152509: 82/82 attempt
  records "completed" (no failure, no hold, no retry needed); with the two
  smoke outputs all 84 task outputs exist.  `assemble --version v2` (the
  full output scan: every report parses, has schema /2, one frozen program
  sha256 910b4a20d366..., one task configuration per sample, exact input
  coverage) published production_v3/manifests/v2/{data.json, mc_*.json,
  summary.json (sha256 2421fe1851313294...), summary.md}:
  - pseudo-data: 83/83 shards valid with every required branch read in full;
    schema pfnano_fake_data/v3 and lumi_fb 20 in every shard; one branch list.
  - MC: all 6,327 files valid.  Coverage of the documented records: ZZTo4L
    0.8042, GGZZ4Mu 1.0000, GGZZ4E 0.0680, GGZZ2E2Mu 0.0120, DYJetsToLL 1.0000,
    TTBar 0.2456 (as in seed_v2.md); signal preselected events ggH 89,245, VBF
    90,076, VH 87,900 (as in seed_v3.md).
  - genWeight/count preselected-fraction ratio on all files: ZZTo4L 1.00094,
    TTBar 1.00065, others 1.00000-1.00001; negative weights 0.46 % (ZZ) and
    0.37 % (TTBar) of the preselected events.
- Plots: production_v3/manifests/v2/plots/ (first version) and plots_r2/
  (corrected: the log-scale file-count plot hid the single GGZZ2E2Mu file;
  bars now annotated and the axis starts at 0.5).  Checked numerically
  against summary.json and inspected (normalization_ratio,
  events_per_file_ZZTo4L, files_per_sample).
- 2026-09-23T23:25Z, stage 2 development: `analysis_v3/common/config/analysis_ul16_v3.json`
  (version skim_v1: trigger list and mask order, skim floors muon pT > 3 /
  electron pT > 5 GeV with abs(eta), dxy < 0.5, dz < 1; pair window 40-140 GeV;
  probe definitions; AN FSR cuts; slim branch patterns; paper Table 4
  fiducial volume) and `analysis_v3/skims/src/skim_v3.cpp`.  Local tests on
  real inputs: ggH (2 files, 3.3 s), DY (2 files, 4.6 s, 8 % of events in the
  slim), ZZTo4L (1 file, 5.1 s, 52 % in the slim), one pseudo-data shard
  (70 s, 856 MB peak; only structural checks of the data output).  Expected
  total skim volume about 50 GB.
- 2026-09-23T23:25Z, fiducial cross-check (ggH test file): the full Table 4
  selection passes 53.3 % of the preselected events and an independent bare-
  lepton check 55.0 %, consistent.  Relative to generated events A_fid =
  0.478, above the paper's 0.398 (POWHEG+JHUGen): in this private Pythia8
  sample only 3.5 % of the kinematically accepted events have m_Z2 < 12 GeV,
  consistent with Pythia's default Z mass floor (23:mMin = 10 GeV) and the
  absence of gamma*; the pseudo-data signal comes from the same generator, so
  the analysis is internally consistent.  To be stated in the AN's paper
  comparison.
- 2026-09-23T23:45Z, stage 2 physics review (before any submission): four A
  findings, all adopted.  (1) Trigger-object matching kept only the closest
  object: UL16 muons carry near-duplicate HLT objects with different bits
  (Iso lost for 61 % of the muons carrying it), and HLT electron objects sit
  at the supercluster position (measured in DY MC: eta = eta_SC to 0.0005,
  q*dphi = -0.750/pT barrel and -0.537/pT endcap, as expected from
  0.3*B*r/2 with B = 3.8 T and r = 1.29 m or z/sinh(eta)).  Fix: merge muon
  objects within dR < 0.05 (OR of the bits), OR the bits of all objects in the
  cone, and match electrons at (eta_SC, phi - q*0.3*B*r/(2 pT)).  (2) The
  TChain::CopyTree slim copy ignored GetEntry errors.  Fix: explicit
  CloneTree(0) loop with fatal GetEntry <= 0, file key and entry written into
  Events, and a ROOT error handler that makes any I/O error during the event
  loops fatal.  (3) Pairs lacked the inputs of a trained electron ID.  Fix:
  store the electron (sieie, hoe, r9, eInvMinusPInv, scEtOverPt, dr03 sums,
  miniPFRelIso, vidNestedWPBitmap, tightCharge, ip3d, dxy/dz errors, eCorr)
  and muon (nTrackerLayers, tunepRelPt, segmentComp, tkRelIso, mvaLowPt,
  softMva, ip3d, mediumPromptId, highPurity) ID inputs.  (4) HTXS is empty in
  the private Pythia8 signal (stage_0 = 0, Higgs_y = NaN).  Fix: STXS stage 0
  and the VH decay class from GenPart and the sample mode.  B/C findings
  adopted: merged TrigObj probes, dressed leptons and candidate indices in
  GenTable, fromHardProcess for the H-decay leptons, AN FSR cuts before the
  lowest dR/ET^2 choice, electron floor abs(eta) OR abs(eta_SC) < 2.5,
  portable Box-Muller deviates.  Analysis configuration now skim_v2.
- 2026-09-24T00:05Z, stage 2 revision checks on MC test inputs (all pass):
  - VH sample content: the private Pythia8 generation forces every Z to
    ee/mumu (seed_v3.md), the associated Z of ZH included, so ZH appears only
    as Z->ll: vh_class W->lnu 878, Z->ll 101, W->qq 1701, Z->qq 0, Z->nunu 0
    of 2680 events (ZH 3.8 %, as expected from sigma(ZH) BR(Z->ee,mumu) /
    sigma(WH)); W decays are not forced (W->qq 66 %).  To be stated in the AN.
  - STXS stage 0 from GenPart: ggH 11 (98 %) / 10; VBF 21 (99 %); VH 23 / 31 /
    41 as above, 1.8 % forward.  H-decay leptons: exactly 4 in 1795/1795 ggH
    events with the fromHardProcess requirement.
  - Trigger matching: 89 % of the matched ZZ muon legs now carry the Iso bit.
    Electron matching at the predicted supercluster position with dR < 0.1
    reproduces the wide uncorrected dR < 0.3 cone (7-10 GeV 0.434 vs 0.441,
    10-15 GeV 0.779 vs 0.782) while the uncorrected dR < 0.1 loses 17 % at
    7-10 GeV.  Merged muon objects: 1.92 raw objects per TrigObj probe in DY.
- 2026-09-24T00:15Z, stage 2 production review (before any submission): two A
  findings.  (1) Slim copy: fixed by the explicit CloneTree(0) loop, now with
  per-file span checks (LoadTree/GetTreeNumber), fatal GetEntry, Fill and
  Write checks, the error handler, and a full read of every output tree
  before publication.  (2) The v1 plan and its two smoke outputs came from
  the reviewed binary: superseded (kept, unused); new plan
  production_v3/skims/v2/plan.json (config skim_v2, sample mode and manifest
  entry counts in every task), and the stager now refuses existing outputs of
  another program version.  B/C fixes: -lGenVector dropped (the shared-library
  list now equals the proven validate_inputs one), unreadable MC files at
  open are skipped and recorded, manifest entry counts asserted, flavour
  microcentury.
- 2026-09-24T00:15Z, smoke of skims v2 through the wrapper (PATH=/usr/bin:/bin):
  VHToZZ_M125 (99 files, 137 s, 456 MB), DYJetsToLL_0000 (39 files, 97 s,
  653 MB), pseudo_data part_00000 (71 s, 955 MB).  The first data smoke
  failed at root_check on the fresh EOS copy although the sha256 matched; a
  manual re-run passed: ROOT reads /eos through XRootD, whose view of a file
  just written through FUSE lagged.  The wrapper now retries root_check
  (6 attempts, growing delay); rule added to AGENTS.md.
- 2026-09-23T23:53Z, SUBMITTED stage 2: cluster 1152510, 280 jobs (283 tasks,
  3 complete from the smoke), run directory
  production_v3/condor/skims_v2/20260923T235319Z, 1 CPU, 2000 MB, 3000 MB
  disk, microcentury, frozen programs 8ede04ddc00f.  Completion gate: all 283
  outputs in production_v3/skims/v2/, then validate_skims.py (full scan,
  coverage, histograms, plots).
- 2026-09-24T00:05Z, stage 2 monitoring: 66 of the first 250 attempts failed
  with exit 73 "cannot create the output directory" (ZZTo4L 58, TTBar 28,
  GGZZ4Mu 1): the directories of samples without a smoke output were created
  by the first job of each sample, and the concurrent jobs of that sample did
  not see them within the wrapper's 15 s window through their EOS FUSE
  mounts.  No hold; Condor retries them (max_retries = 2) now that the
  directories exist.  Root-cause fix for the next submissions: the stager
  creates every output directory before submission; the wrapper waits up to
  8 attempts with growing delays.
- 2026-09-24T00:27Z, stage 2: job 1152510.193 (skim_TTBar_0008, node
  b9p29p4279) ran 27 min against about 3 min for its peers; its ROOT copy
  was complete on EOS at 00:13Z (389 MB, opens cleanly) but the job did not
  finish, most likely stalled in an XRootD read of the root_check step.
  Removed with condor_rm (positively identified job of this workflow); the
  task is re-run locally through the same frozen wrapper; the partial copy is
  kept under its attempt-qualified name.
- 2026-09-24T00:45Z, STAGE 2 DONE.  All 283 skim outputs of production_v3/skims/v2/
  exist (280 from cluster 1152510; 91 first attempts had failed on the EOS
  directory race and were retried by Condor successfully; job 193 stalled
  and was removed, its task and the three smoke tasks ran locally through
  the same frozen program).  validate_skims.py: full scan 283/283 valid (one
  program sha256, one configuration, ROOT publication records, root_check on
  every output), coverage complete (every manifest file processed, no MC
  skip).  Histograms and plots: production_v3/skims/v2/validation/ and
  production_v3/skims/v2/plots_r2/ (plots/ is the first version with a fixed
  0.8-1.2 ratio range that hid the Z-peak ratio; plots_r2 has an adaptive
  range).  Numerical checks: data/MC pair totals mumu 0.996, ee 0.993; in
  81-101 GeV 0.986 and 0.975.  The pseudo-data Z peaks sit visibly below the
  MC (mumu about 89.3 against 91.0 GeV), a lepton momentum-scale difference
  for stage 3 to measure; N_PV data/MC flat at about 0.99 (no pileup
  difference: no reweighting); muon ID flags in the Z window: data/MC about
  0.98, mediumId/tightId/mediumPromptId about 0.97.
- 2026-09-24T00:57Z, stage 3a development (lepton momentum scale and
  resolution): code in `analysis_v3/calibration/` (`calib_extract`,
  `zpeak_histograms`, `fit_zpeaks`; shared `analysis_v3/common/include/h4l/calibration.h`
  with the bin lookup (lower edges, half-open, unbounded terminal bins), the
  named lepton IDs and the payload application; configuration
  `analysis_v3/calibration/config/calibration_ul16_v1.json`).  Extraction
  plan `production_v3/calibration/v1/extract/plan.json`: 43 tasks over the
  283 skims v2 outputs, analysis trigger-OR mask 0x1fffff, raw mass 50-130
  GeV, both legs iso03 < 1 and SIP < 10.  Local smoke of
  calx_pseudo_data_0000 through the frozen wrapper (PATH=/usr/bin:/bin):
  5 shards, 2,087,440 pairs read, 1,875,410 kept, 24 s, 478 MB peak;
  published.  Two independent reviews (physics; production) started before
  the Condor submission.
- 2026-09-24T01:10Z, finding: the lxplus system ROOT 6.40.04 has no FFTW
  plugin (`root-fftw` is not installed; "TVirtualFFT::FFT plugin not
  found"), so RooFFTConvPdf cannot run here or on workers.  Substitute: the
  analysis's own numerical BW (x) DCB convolution on a 0.05 GeV grid with
  the Breit-Wigner tabulated once (`h4l/zpeak_model.h`) and a Minuit2 binned
  extended Poisson fitter with HESSE errors and the sandwich covariance for
  weighted histograms (`h4l/binned_fit.h`); about 1-2 s per category.
- 2026-09-24T01:25Z, finding (test fits on calx_pseudo_data_0000, about 6 %
  of the data, against 1/11 of the DY MC): in the (pT, eta) x (pT, eta)
  leg-bin categories the pT bin edges sculpt the dilepton lineshape (both
  muons 40-50 GeV and abs(eta) < 0.9: the far tails fall far below any
  Breit-Wigner).  BW (x) DCB + exponential gives chi2/ndf about 180 in the
  window [mode - 25, mode + 20] and about 23 in +-8 GeV, with degenerate
  tails.  First look at the peaks: data/MC about 0.983 in both flavours
  (e.g. mumu 40-50 GeV central: 89.55 against 91.11 GeV).  Question to the
  user: which method for the per-lepton scale/smear extraction (MC template
  fit; AN-literal iterative BW (x) DCB; BW (x) DCB in eta-only categories).
- 2026-09-24T13:00Z, USER ANSWER (verbatim): "可以用MC模版拟合。但是我也要提醒你，
  对于拟合质量，要通过看图，而不是chi2或者pull，因为entry很多就很难做到chi2、pull很好。
  重点是看起来拟合上了，合理即可，并且即使chi2、pull比较好，也还是要看起来合理拟合上了才行。
  另外对于拟合应该怎么拟合，用什么函数、参数、限制等，可以参考/afs/cern.ch/user/y/yiyangz/Research/codex/slopbench_code_fork/results/h4l_seeds/efficiency_tnp_an里的参数，
  但是可能分bin不完全一样".  Decisions: (a) the per-lepton scale and smear
  come from MC template fits in the leg-bin categories, iterated with the
  data corrected and the MC smeared per lepton; BW (x) DCB fits stay for the
  AN-style validation and the lambda calibration; (b) fit quality is judged
  by looking at the plotted fits, never by chi2/pulls alone (durable rule in
  AGENTS.md); (c) the fit functions, parameters and limits of
  `efficiency_tnp_an` may be consulted (a user-authorized exception to the
  h4l_seeds blinding rule, limited to the fit configuration).
- 2026-09-24T13:10Z, read the fit configuration of `efficiency_tnp_an`
  (`tnp_validation.tex` sections "Tag and probe" and "The fit", the record
  structure of `fits/mc_dy.json`, and the `mass`/`binning`/`fit` sections of
  its configuration `slopbench_code/benchmark/configs/h4l_tnp.json`).
  Conventions taken over where applicable (binning may differ): window
  60-120 GeV in 0.5 GeV storage bins fitted at 1 GeV (coarsened up to 3 GeV
  for at least 40 entries per bin); pass and fail fitted simultaneously with
  the efficiency as a parameter; nominal BW (PDG m_Z, Gamma_Z fixed) (x)
  Crystal Ball with mass shift, width, alpha and n shared by pass and fail,
  exponential backgrounds per category; alternatives BW (x) two Gaussians
  sharing one mean and a monotonically falling Bernstein-2 background;
  starting width base (1 + a <|eta|>)(1 + b <pT>/100 GeV) (1.21 GeV muons,
  1.89 GeV electrons inclusive), bounded by a factor 3 and clipped to
  [0.3, 8] GeV, a width on its bound refitted with the bound opened; staged
  freezing by effective probes (> 3000 all free; < 3000 fix n; < 1000 also
  alpha; < 300 also shift and width, at the inclusive values); rejection of a
  total curve with more than one maximum or a background peaking within
  12 GeV of m_Z; HESSE failure repeated with the resolution fixed;
  convergence judged on the EDM; every fit drawn.  One printed DY record
  incidentally showed a field of predicted ratios for the injected profiles
  of the older seeds; it concerns another dataset family and is not used
  anywhere.  The seed results and paired closures of that directory were not
  read.
- 2026-09-24T13:15Z, stage 3a reviews (both before the first Condor
  submission of stage 3).  Production review: no A finding for the
  submission.  To do before any use: a full scan of the extraction outputs
  (A1: `analysis_v3/calibration/scripts/scan_extract.py`, written) and
  per-input row checks in the histogram jobs (A2, done).  B findings
  adopted: a new configuration version instead of in-place edits
  (calibration_ul16_v2.json, calib_v2; the v1 file, embedded in the v1 plan
  with sha256 2a36888b / fnv1a64 46bd5763 and edited twice in place, is
  labelled superseded, sha256 now d5a24e35; the v2 plan embeds only the
  extract section and the report records its hash); input sha256 in the
  extraction report; bounded retries of the output check and an exclusive
  `.lock` reservation (exit 17) in zpeak_histograms; frozen program identity
  for the local runs; the closure-injection crash below the first pT edge
  (removed by the interpolated payload); `timeout` around the wrapper's EOS
  copy (1800 s), checksum (900 s) and root_check (600 s), the assertion for
  the stage-2 stall of job 1152510.193.  Physics review: extraction sound;
  adopted: more CalibPairs fields (pair pT/y, FSR mass and photons, iso03_chg
  and the FSR-subtracted isolation used for both the preselection and the
  tight selection, dxy/dz, r9, gen_flav, leg trigger bits, npv); the lowest
  pT bins split at 10 and 15 GeV; the scale and smear applied by bilinear
  interpolation between the bin mean positions (a piecewise-constant scale
  would turn a within-bin scale slope into a spurious smear); windows kept
  inside 60-120 GeV (DY generator cut at 50 GeV); the template design:
  event-level MC templates with frozen pair deviates and linear bin sharing
  (no smoothing of the binned MC, which biases the smear), MC statistics in
  the likelihood (Barlow-Beeston lite), a common known data and MC smear
  delta so that negative residuals are measurable, iteration with per-lepton
  re-categorization, per-category kinematics plots, a second calibration
  pass with the tag-and-probe SFs on the MC, and closure on the same events
  and on independent halves with an injection grid.
- 2026-09-24T13:21Z, SUBMITTED stage 3a extraction v2: cluster 1153910,
  41 jobs (43 tasks; calx_TTBar_0000 and calx_pseudo_data_0000 complete from
  the smokes: 34 s, 474 and 540 MB), run directory
  production_v3/condor/calib_extract_v2/20260924T132145Z, frozen programs
  f6839018f744 (wrapper with timeouts), 1 CPU, 2000 MB, 1000 MB disk,
  microcentury.  Plan production_v3/calibration/v2/extract/plan.json (the v1
  plan and its single smoke output are superseded and unused).  Completion
  gate: all 43 outputs, then scan_extract.py --version v2 --manifests v2.
- 2026-09-24T13:40Z, stage 3a extraction v2 DONE: cluster 1153910, 41/41
  attempt records "completed", no failure, no hold, no partial or orphan file.
  Full scan (`scan_extract.py --version v2 --manifests v2`): 43/43 outputs
  bound to their plan tasks (inputs with sha256, trigger mask 0x1fffff,
  extraction-config hash 5dffadcedd316b9a), one program (sha256 1e3a7173...),
  every ROOT sha256 recomputed equal to its publication record, root_check
  with the JSON tree counts, CalibPairs = kept = kept_mm + kept_ee, and the
  per-sample original-file sets equal to the v2 manifests (complete
  coverage).  Pairs: pseudo_data 30,794,553; DY 23,484,358; ZZTo4L
  11,364,393; GGZZ4Mu 2,269,261; TTBar 1,444,886; GGZZ4E 135,629; signals
  84,914 / 88,268 / 91,474; GGZZ2E2Mu 7,909; 6.5 GB.  Scan record:
  production_v3/calibration/v2/extract/scan.json.
- 2026-09-24T14:20Z, stage 3a template-fit development (tests on about 12 %
  of the data against 27 % of the DY MC, then the full run):
  - Linear bin sharing of the event-level template gave a wiggly likelihood
    (MC-noise structure of about one NLL unit, false local minima in D,
    HESSE errors from the curvature of a wiggle; seen in likelihood slices).
    Fix: every MC pair is spread by a fixed Gaussian kernel of relative width
    0.3 % (smooth likelihood; the smoothing does not vary with D, so it cannot
    favour a larger smear); the reported E = D + kernel^2.  The common data
    smear delta raised to 0.5 % so that D reaches -2.5e-5.  Slices now smooth
    parabolas; fitted templates visually describe the data in every tested
    category.  Empty template bins in the tails met by data events (infinite
    likelihood) fixed by likelihood groups of consecutive 0.1 GeV bins with at
    least 20 effective MC entries (fixed per fit, the template stays fine).
  - Kinematic-edge categories: in e.g. mumu 17_21 (30-40 x 40-50 GeV, both
    0.9-1.5) the pT-bin boundaries give a sharp mass edge fixed in
    reconstructed pT, so the pair-level scale does not respond to the lepton
    scale (fit at ln k about 0, D at its bound).  Response passes added (MC
    playing the data with a known -0.5 % scale or 0.8 % smear, same events as
    the MC role): the good diagonal categories respond with about 0.9 (mumu
    20_20: 0.89 scale, 0.89 smear), but low-statistics and edge categories
    give responses from about -3 to 13 (edges and the noise of a 0.5 %
    probe); the cut R >= 0.5 excluded more than half of the categories.
- 2026-09-24T14:20Z, nominal calibration run production_v3/calibration/v2/nominal
  (tmux session calib_nominal; frozen programs under
  production_v3/program_inputs/calibration_v2/): iteration 0 scale residuals
  -1.4 to -2.4 % per lepton (e.g. muons 40-50 GeV |eta| < 0.9 -1.71 %,
  electrons 40-50 GeV |eta_SC| < 1.0 -1.65 %), smear residual about +8e-5 per
  lepton centrally (r about 0.9 %); worst residual 70.6, 16.3, 6.7, 4.2 sigma
  in iterations 0-3 (with the errors inflated by chi2/ndf).  FINDING: the
  per-lepton least squares has chi2/ndf about 15-40 even after the first
  corrections; the largest pulls are high-statistics categories pairing a
  central leg (|eta| < 0.9) with one at 0.9-1.5 (mumu 20_21, 651k pairs:
  residual ln k -7.8e-4 +- 0.4e-4), while the diagonal categories are
  consistent.  A model with a per-lepton scale depending only on (pT, |eta|)
  does not describe the data at this precision; charge-, signed-eta- or
  phi-dependent effects (the Rochester-style parameters of the AN muon
  corrections) would average out in |eta| bins but leave exactly such
  pair-level inconsistencies.  Next: diagnostics of the data/MC mass shift
  against each charge's signed eta and phi (diagnose_calibration.py, FullPairs
  with the leg phi added).
- 2026-09-24T14:40Z, DIAGNOSTICS (production_v3/calibration/v2/nominal/diag_iter04,
  payload of iteration 3; data/MC ln k of pairs with a given charge's lepton in
  slices of signed eta and of phi): no charge dependence and no phi
  dependence for either flavour (mu+/mu- and e+/e- agree within errors; flat
  in phi), but a clear |eta| structure symmetric in +-eta inside the
  calibration bins: muons oscillate by about +-0.3e-3 in pair ln k (about
  +-0.06 % per lepton: positive near |eta| 0.2 and 1.8, negative near 1.0 and
  2.2); electrons show a sharp dip of about -0.9e-3 at |eta| 1.4-1.6 (the
  barrel-endcap transition).  After several iterations the high-statistics
  bins of 20-100 GeV are nearly flat in pT (muons |eta| < 0.9 about -1.9 %);
  low-statistics bins oscillate.  The nominal v2 run does not converge
  (worst residual 4.2, 4.6?, 5.2, 9.5 sigma in iterations 3-6): the
  (pT, |eta|) model with 4-5 |eta| bins cannot describe the data.
- 2026-09-24T14:45Z, QUESTION to the user (calibration model: factorized fine
  |eta| + coarse pT; fine 2D grid; keep the bins with a systematic).  USER
  ANSWER: "分解模型：细η+粗pT (Recommended)".  Decision: ln(1 + s) =
  a(|eta|) + b_R(pT), r^2 = c(|eta|) + d_R(pT); a, c on fine |eta| bins (muons
  11 bins of about 0.2 with 1.2 at a bin centre; electrons 13 bins in
  |eta_SC| keeping the 1.4442 / 1.566 edges), measured in pT-inclusive
  (eta_i, eta_j) categories (both legs pT > 20 GeV, no pT-bin sculpting);
  b, d the pT dependence per coarse |eta| region (muons 3, electrons 2),
  measured in (pT, region) x (pT, region) categories, zero at the 40-50 GeV
  reference bin; one joint least squares of both category families with the
  measured leg compositions, iterated.  The v2 nominal run is superseded and
  kept.
- 2026-09-24T15:30Z, factorized-model calibration production_v3/calibration/v3/nominal
  (configuration calibration_ul16_v3.json, calib_v3): iteration 0 per-lepton
  chi2/ndf mumu scale 55208/299 (family A 15114 over 66 categories, B 40094
  over 268, before any response); a family-A-only per-lepton eta fit of the
  same iteration gives chi2/ndf 502/55 (mumu) and 362/73 (ee) and a smooth,
  nearly flavour-independent profile (about -1.88 % at |eta| < 0.3 to
  -2.2 % above 2.1).  Responses (iteration 1, scale -1 %, smear 1.2 %): family
  A 1.00 for the scale (5-95 %: 0.99-1.01) and about 1 for the smear; family
  B broad (pT sculpting).  Iteration 2: family A consistent (mumu scale chi2
  30 over 61 categories, ee 36 over 61); the rest is family B, dominated by
  mumu B_3_16 (10-15 GeV central x 40-50 GeV at 0.9-1.5; pull 34): its data
  mode is about 67 GeV, a kinematic structure of its pT bins, not the Z
  peak.  Fixes (configuration relabelled calib_v3.1): only categories whose
  data mode lies in 80-100 GeV enter; iterative outlier rejection in the
  joint least squares (|pull| > 5, at most 10 % of the categories,
  recorded).  The v3/nominal run was stopped after iteration 2 (kept);
  restarted as production_v3/calibration/v3/nominal_r2.
- 2026-09-24T16:25Z, v3/nominal_r2 (Z-mode selection, outlier rejection):
  per-lepton chi2/ndf for the scale 1.15 (mumu) and 1.4 (ee) from
  iteration 2, smear about 3.4; worst residual 25.5, 7.5, 5.5, 5.8 sigma in
  iterations 1-4.  The core is converged and physical: the |eta| profile
  (muons -1.878 % at |eta| < 0.3 to -2.220 % above 2.1; electrons -1.884 % to
  -2.286 %, gap bin -2.018 %), a consistent pT slope above 20 GeV (the scale
  more negative by about 0.2 % from 30 to 135 GeV), smear c about 1.0-1.1e-4
  per lepton (r about 1.0-1.05 %).  Not converging: the pT terms below 20 GeV
  (and some smear terms), each resting on few strongly sculpted family-B
  categories (Z decays with a soft leg), which jump when the category
  selection changes between iterations (smear terms up to +-4e-3 with equal
  errors).  Fixes (configuration calib_v3.2): a random-walk smoothness prior
  between neighbouring pT bins of each region on the total pT terms (tau
  1e-3 for ln(1+s), 3e-5 for r^2; its strength becomes a systematic) and the
  category set frozen from iteration 2 on.  r2 stopped (kept); restarted as
  production_v3/calibration/v3/nominal_r3.
- 2026-09-24T17:10Z, v3/nominal_r3 (smoothness prior, frozen selection from
  iteration 2): worst residual 162, 35.2, 10.7, 3.4, 4.8 sigma in iterations
  0-4; scale chi2/ndf about 0.8-1.3 (mumu) and 1.4-2.9 (ee), smear 3.9-4.1
  (mumu) and 6-9 (ee).  Muon parameters converge geometrically (factor about
  0.35 per iteration).  The pT terms of endcap electrons below 20 GeV
  oscillate between iterations (e.g. 5-10 GeV: +2.8, -2.7, +3.3, +1.8,
  -4.8 sigma): their few sculpted categories respond discontinuously to the
  data-seeded windows moving between iterations.  Added
  finalize_calibration.py: for a run that reaches max_iterations without
  meeting the 0.2-sigma criterion, the final payload is the mean of the
  payloads applied in the last 4 iterations, and each parameter's iteration
  rms is added in quadrature to its statistical uncertainty (plot and
  closure reports include it).
- 2026-09-24T17:45Z, stage 3a nominal_r3 reached max_iterations (worst
  residual 2.1 sigma at iteration 10; muon parameters within about 1 sigma
  from iteration 5; the endcap-electron pT terms below 20 GeV and above
  100 GeV oscillate).  finalize_calibration.py: payload = mean of the
  payloads applied in iterations 7-10, iteration rms up to 1.4e-4 (muon
  scale), 2.2e-3 (electron scale, 5-10 GeV endcap), 1.4e-5 / 4.3e-5 (smear
  variance).  production_v3/calibration/v3/nominal_r3/payload.json.
- 2026-09-24T17:50Z, stage 3c tag-and-probe extraction: `analysis_v3/tnp/`
  (tnp_extract.cpp, make_tnp_plan.py, config tnp_ul16_v1.json).  Pairs from
  the skim Pairs trees with raw mass 50-130 GeV and at least one tag-candidate
  leg (tag path of the flavour fired: HLT_IsoMu24/IsoTkMu24, filter bits
  Iso|IsoTkMu = mask 10, trigger pT >= 22, raw pT >= 20; or
  HLT_Ele27_WPTight_Gsf/Ele25_eta2p1_WPTight_Gsf, filter bit 1e WPTight =
  mask 2, trigger pT >= 25, raw pT >= 25), no ID/isolation/SIP requirement.
  Substitution: the unprescaled 2016 single-muon tag paths IsoMu24/IsoTkMu24
  instead of the AN IsoMu20/22, and Ele27/Ele25_eta2p1 WPTight (Ele27_eta2p1
  is not in the skim bits).  Independent review (physics + production): no A
  finding; masks verified against the CMSSW 10_6 NanoAOD trigger-object
  definitions; recomputation of the tag rule on the smoke: 0 mismatches.
  Adopted: the electron cut-based ID and MVA inputs and the muon ID inputs
  stored (a stage-4b ID choice needs no new extraction), hence plan version
  v2 (v1 and its two smoke outputs superseded and kept); to do before use: a
  full scan (scan of the tnp outputs), validation plots, and the
  path-consistent final tag (IsoMu24 with bit 2 or IsoTkMu24 with bit 8 and
  trigger pT >= 24; Ele27 with trigger pT >= 27 or Ele25_eta2p1 with
  |eta_SC| < 2.1) with the plateau cut.  Smokes of v2 (tnpx_TTBar_0000 and
  tnpx_pseudo_data_0000): 57 and 72 s, 521 and 671 MB; 1,722,939 of
  2,087,440 data pairs kept.
- 2026-09-24T17:57Z, SUBMITTED stage 3c T&P extraction v2: cluster 1153950,
  41 jobs, frozen programs 14c042679a14, 1 CPU, 1000 MB, 2000 MB disk,
  espresso; plan production_v3/tnp/v2/extract/plan.json.  Completion gate:
  all 43 outputs, then the full output scan.
- 2026-09-24T18:15Z, stage 3c T&P extraction v2 DONE: cluster 1153950, 41/41
  attempt records "completed"; full scan (scan_tnp_extract.py --version v2
  --manifests v2): 43/43 outputs valid (plan binding with input sha256, tag
  masks, extraction-config hash, one program, ROOT sha256 recomputed,
  root_check, row counts), coverage complete.  Pairs: pseudo_data
  28,280,225; DY 21,539,742; ZZTo4L 10,116,915; GGZZ4Mu 2,193,013; TTBar
  1,810,539.  Scan record production_v3/tnp/v2/extract/scan.json.
- 2026-09-24T18:15Z, nominal_r3 smear diagnosis: the electron smear
  inconsistency (chi2/ndf about 10) is localized at the barrel-endcap
  transition: the largest pulls involve the fine bin |eta_SC| 1.566-1.7 (and
  the gap bin), whose apparent extra smear depends on the partner (A_1_8
  E = +3.4e-5 +- 0.7e-5, A_5_8 E = -1.45e-5 +- 0.21e-5), which no per-lepton
  Gaussian smear satisfies.  Covered by the chi2/ndf error inflation (factor
  about 3.2 for the electron smear) and the iteration rms; to be studied in
  the closure runs.
- 2026-09-24T18:40Z, stage 3a nominal RESULT (production_v3/calibration/v3/nominal_r3,
  finalized payload; plots and report.json in plots_r2/, plots/ is an
  interrupted first attempt of the slow kinematics pages): at pT 45 GeV,
  |eta| 1.2: muon scale shift -2.0270 % +- 0.0052 % (stat) +- 0.0002 %
  (iteration), muon smear 1.070 % +- 0.040 %; electron scale shift
  -2.0034 % +- 0.0152 % +- 0.0016 %, electron smear 0.951 % +- 0.142 % +-
  0.046 % (the electron smear error carries the chi2/ndf inflation of the
  transition-region inconsistency).  Profiles: the scale falls smoothly and
  monotonically with |eta| (muons -1.88 % to -2.22 %, electrons -1.88 % to
  -2.28 %), a pT slope above 20 GeV, the smear about 1.0-1.1 % (barrel) with
  lower values in the far endcap for electrons.  Validation: the inclusive
  calibrated Z peaks of data and MC agree over 60-120 GeV (ratio flat at 1),
  against a 1.8 GeV offset and a wider data peak before; the MC
  normalization in 60-120 GeV is x0.982 (mumu) and x0.969 (ee): data/MC
  yield deficits for the tag-and-probe scale factors.  Galleries of every
  template fit (iterations 0 and 10, both families) and per-category leg
  and pair pT data/MC pages are in plots_r2/.  Systematics (closure, fit
  variants, prior strength) pending.
- 2026-09-24T18:50Z, USER (verbatim): "继续啊，你的工作又不是只有刻度，是整个分析啊，还有最后的详细AN。完成全部prompt_new的工作，不要停止啊".
  Decision: work continues without pausing through every stage of
  prompt_renew.md (tag-and-probe, lambda, 4l reconstruction and
  optimization, backgrounds, signal model, statistical model, all section-6
  results, the benchmark files, validation toys, the detailed AN); progress
  is reported briefly; a genuine method choice is asked while other work
  continues.
- 2026-09-24T18:52Z, USER (verbatim): "记得记录log、要求、plan，及时回顾，从而避免忘记了！".
  Decision: every requirement goes into AGENTS.md or this log at once, the
  plan (PLAN.md) is updated whenever the work plan changes, and AGENTS.md,
  PLAN.md and the recent log are re-read before every stage and at least at
  every stage boundary.
- 2026-09-24T18:26Z, TIMESTAMP CORRECTION: the three entries above dated 18:40Z, 18:50Z and
  18:52Z were written at about 18:05Z, 18:15Z and 18:20Z (the times were
  estimated instead of read from the clock).  From now on every log
  timestamp is taken from `date -u` in the same command.
- 2026-09-24T18:26Z, stage 3c tag-and-probe run started: production_v3/tnp/v2/nominal
  (tmux tnp_nominal; run_tnp.py with the calibration payload v3/nominal_r3,
  configuration tnp_ul16_v2.json; programs tnp_histograms and fit_tnp
  frozen under production_v3/program_inputs/tnp_v2/).  MELA build delegated
  to a background agent (production_v3/external/).
- 2026-09-24T19:10Z, stage 3a final-payload diagnostics (diagnose_calibration.py --final,
  production_v3/calibration/v3/nominal_r3/diag_final/plots/charge_eta_phi_{mm,ee}.png):
  residual ln k(data/MC) of pairs with a lepton of given charge in eta or phi
  slices within 1.1e-4 (muons) and 2e-4 (electrons), no charge or phi
  structure.  The data/MC interquartile-range ratio (82-100 GeV) is flat
  within 0.5 % for electrons but about 1.013 for muons with |eta| < 1 (1.002 in
  the endcaps): a possible residual muon width difference in the barrel, to be
  checked against the lambda fits (data/MC sigma per region) and the closure
  runs.
- 2026-09-24T19:10Z, stage 3b lambda v1 (run_lambda.py with lambda_ul16_v1.json,
  production_v3/calibration/v3/nominal_r3/lambda_v1/, SUPERSEDED): independent
  BW (x) DCB fits per bin of the predicted relative mass error and a slope
  through the origin.  Rejected on the fit galleries: the per-bin DCB sigma is
  degenerate with the free tails (e.g. lm_mm_1_11 sigma 0.71 GeV for a peak of
  about 3 GeV width), the high-error electron classes have almost no pairs with
  both legs in the class (<= 600 per bin), and the first/last e bins dropped
  pairs outside the edges.
- 2026-09-24T19:10Z, stage 3b lambda v2 (lambda_ul16_v2.json, new program fit_lambda.cpp,
  lambda_histograms v2): the binned form of the AN 5.3.1 conditional fit.  Pairs
  are classified by the lambda regions (a <= b) of their legs and binned in e =
  0.5 sqrt(d1^2 + d2^2) (first bin from 0, last unbounded); all e bins of a class
  are fitted simultaneously with BW (x) DCB + exponential, shared DCB tails, a
  free shift, yields and slope per bin and sigma_k = (m_Z/2) sqrt(lambda_a^2
  <d_a^2>_k + lambda_b^2 <d_b^2>_k); MC first, the data with the MC tails; a
  per-bin closure refit with sigma free.  Regions as the AN: muons |eta|
  0-0.9-1.8-2.4 (both legs in the region); electrons (|eta_SC|) d < 0.03 in
  0-0.8, 0.8-1.0 and d < 0.07 in 1.0-1.2-1.44-1.57-2.0-2.5 (both legs, the gap
  bin with the reference method), and the high-error classes |eta| < 1, d >
  0.03 and |eta| > 1, d > 0.07 with the AN reference-electron method (one leg in
  |eta| < 0.8, d < 0.03 with its lambda fixed).  Substitution: no ecalDriven flag
  in UL NanoAOD, so no separate tracker-driven corrections.  First run
  (label lambda_v2, kept): all fits converged except the MC fit of |eta_SC|
  1.57-2.0 (Minuit status 1 with n_L at its limit 80), which stopped the run;
  root-cause fix in fit_lambda: a free tail parameter ending at its limit is
  fixed there and the fit repeated once (staged freezing, recorded as
  fixed_at_limit); the reported chi2 now floors the variance at the expectation
  times the mean weight (mixed-weight MC bins gave chi2 ~ 1e10).  Rerun as label
  lambda_v2_r2.  First-run values (data / MC): muons 1.472/1.497, 1.268/1.292,
  1.058/1.087; electrons 1.543/1.585 (0-0.8), 1.473/1.502, 1.548/1.599,
  1.553/1.588, 1.597/1.648 (gap, reference), 1.330/1.348 (2.0-2.5),
  1.187/1.216 (|eta| < 1, d > 0.03), 0.508/0.513 (|eta| > 1, d > 0.07, 3 bins).
  lambda > 1 is expected: the NanoAOD errors do not contain the extra smearing
  of the data (and of the calibrated MC).
- 2026-09-24T19:10Z, stage 4a 4-lepton event records: new program
  analysis_v3/reconstruction/src/h4l_reco.cpp, configuration
  analysis_v3/reconstruction/config/reco_ul16_v1.json (loosest floors: muon /
  electron calibrated pT > 3 / 5 GeV, SIP < 8, Z1 30-120, Z2 4-120, pT 15/7,
  m(OS) > 2, m4l > 70; events with >= 3 loose leptons and an OS SF pair 30-120),
  plan builder make_reco_plan.py.  Design (store all candidates): every loose
  lepton with calibrated and raw pT, deviate, error, FSR-subtracted isolation,
  IDs, FSR photon, HLT match and generator match; jets (tight ID, CSVv2 and
  DeepJet); every ZZ candidate with both orderings of OS-OS pairs
  (c_z1_closer), SS Z2 for the SS region, the smart-cut alternative pairing;
  the electron-muon cross cleaning is left to the final selection through
  l_overlap; theory weights for MC events with a candidate; generator-table
  fields for signal.  Development runs: GGZZ2E2Mu 5,729 slim events -> 4,989
  kept, 3,678 with a candidate; ggH 80,173 -> 73,718 kept, 49,796 with a
  candidate (10 s); a quick AN selection (tight ID, iso < 0.35, SIP < 4, Z1
  40-120, Z2 12-120, pT 20/10, m(OS) > 4, smart cut, trigger, first surviving
  candidate) keeps 26,794 ggH events (27 % of the 99,000 generated): 4mu
  10,700, 2e2mu 12,484, 4e 3,610; m4l median 124.3 GeV.  Plan v1 (74 tasks,
  <= 400k slim events each) and the local smoke through the frozen wrapper
  (reco_ZZTo4L_0000: 307,516 slim events, 236,507 kept, 85,138 with a candidate,
  22.8 s wall, 524 MB peak, 94 MB output) passed; 16 candidates overflowed the
  64-candidate buffer, so kMaxCandidates is now 256 (rebuilt; submission will
  use a new plan version v2; the v1 smoke output stays as a superseded
  artifact).  Two independent reviews (physics; production and Condor safety)
  launched before the first submission.
- 2026-09-24T19:19Z, MELA built (background agent): JHUGenMELA v2.4.4 (commit 10d36ce)
  standalone against the system ROOT 6.40.04 and gcc/gfortran 11.5, MCFM from
  the release's precompiled libmcfm_711.so (el9_amd64_gcc12 URL; the gcc-11 URL
  of setup.sh returns 404; ABI checked); three createIterator loops patched for
  ROOT 6.40 (patches/JHUGenMELA-v2.4.4_ROOT-6.40_createIterator.patch); own
  idempotent build script (setup.sh deletes files and canonicalizes paths);
  MELA's testME_all.py: 503/503 references identical.  Wrapper
  production_v3/external/mela_wrapper/ (MelaWrapper.h, lib/libMelaWrapper.so,
  README.md): h4lmela::init(13, 125), compute(input, output, withProduction);
  timing 0.13 ms (decay only), 2 ms (2 jets), 10 ms (1 jet, eta_J
  integration).  Runtime: writable working directory (MELA creates links
  there) and the /eos/user FUSE mount (compiled-in data path).
- 2026-09-24T19:19Z, QUESTION (AskUserQuestion) and USER ANSWERS:
  (1) discriminant convention, since the CMS-2016 c-constants are not in
  JHUGenMELA and the AN working points refer to the 2016 normalization:
  answer "MELA normalization (Recommended)": MELA's own per-process
  normalization (getConstant / pAvgSmooth) for D_bkg^kin, D_2jet and D_1jet
  with the paper working point 0.5; for D_WH and D_ZH (unbalanced, about 1e-4)
  the constant is derived on our MC with the same average-probability
  normalization; raw probabilities are stored.  (2) b tagger: answer "DeepJet
  medium" (Jet_btagDeepFlavB, UL16 post-VFP medium WP 0.2489 as decided on
  2026-09-23), no POG scale factors.
- 2026-09-24T19:23Z, USER (verbatim): "你注意TnP应该呢你的实际筛选有关，比如对muon的筛选如果是TRG+Soft，那应该用TnP算TRG+Soft的eff，或者是先算TRG eff，再用过了TRG的probe算soft eff（或者反过来），而不是单独算两个乘起来，那样不对哈".
  Rule recorded in AGENTS.md: every efficiency follows the actual selection,
  either measured directly for the combined requirement or as a chain of
  conditional efficiencies whose denominators are the previous numerators;
  never a product of independently measured efficiencies.  Check of the
  current implementation (tnp_histograms.cpp lines 304-307): id | loose, iso |
  id, sip | id and iso (a conditional chain), full = id and iso and sip | loose
  measured directly as the cross-check; the reconstruction step will be the
  first link (its numerator is exactly the loose probe of the id step), and the
  event trigger efficiency is measured on events passing the complete offline
  selection.
- 2026-09-24T19:47Z, stage 3b lambda DECISION: nominal payload
  production_v3/calibration/v3/nominal_r3/lambda_v2_r2/lambda.json (AN window
  60-120 GeV; closure plots closure_{mm,ee}.png: the corrected points lie on
  the diagonal within the AN 20 % band except sparse tail bins; the gallery shows
  well-described cores, while the tails of the lowest-e bins deviate because
  binning in the predicted error sculpts the pT, hence the mass tails).  A
  core-window variant (lambda_ul16_v2b.json, 82-100 GeV, label lambda_v2b)
  converges only partly (free right tails unconstrained in the narrow window);
  where it converges its lambda differs by 1-3 % (e.g. muon barrel 1.460 vs
  1.472): used as the lambda method uncertainty, far inside the paper's 20 %
  resolution uncertainty.
- 2026-09-24T19:47Z, refit inputs v1 (make_refit_inputs.py): FSR photon resolution from 4,926
  signal-MC FSR photons matched to generator photons: sigma/pT = sqrt(a^2/pT +
  b^2), a = 0.232, b = 0.041 (12 % at 3 GeV, 5 % at 30 GeV).  The v1 Z1 line
  shape (histogram smoothed by 0.2 GeV) made the refit likelihood rough in the
  off-shell tail: in a ggH test of h4l_select 10,010 of 63,248 refits did not
  converge.  v2: the generator Z1 mass fitted with a relativistic BW times
  exp(Chebyshev 1-6) (running); h4l_select now records the refit status per
  candidate and accepts a converged MIGRAD with a covariance.
- 2026-09-24T19:47Z, stage 4a h4l_select (new program, analysis_v3/reconstruction/src/h4l_select.cpp,
  config selection_ul16_v1.json = the AN selection): SR and SRZ4l (Z2 > 12 /
  4 GeV) with the best candidate by D_bkg^kin (MELA decay-only; same four
  leptons: Z1 closest to m_Z), CR 2P2F / 3P1F (best by D_bkg^kin) and SS (Z1
  closest, then the largest Z2 pT sum; events of the SR excluded), ZL (Z1 of
  selected leptons within 7 GeV, pT 20/10, exactly one additional loose lepton,
  m(probe, OS tag) > 4); per candidate the lambda-corrected per-event mass
  error, the Z1 refit, MELA probabilities/discriminants with the cleaned jets,
  DeepJet b tags, additional leptons, MET, the MC m4l responses to the lepton
  scale and smearing; the MELA libraries are verified by sha256 before init.
  ggH development run (73,718 records): SR 26,796, SRZ4l 27,687, 2P2F 846, 3P1F
  6,868, SS 1,051, ZL 12,199; 424 s (the 1-jet MELA integration dominates).
- 2026-09-24T19:47Z, stage 4a production review (independent agent) of h4l_reco: A1
  use-after-free (the input TFile was closed while its TTreeReader was alive;
  valgrind: 3 invalid accesses per input) FIXED in h4l_reco and in h4l_select
  (the reader lives in a block that ends before the file is closed).  B/C
  items fixed: the program runs under a per-flavour timeout in run_task.sh
  (program_timeout_s set by the stager; exit 124 recorded) and a SIGTERM leaves
  an attempt record; the stager checks X509_USER_PROXY and >= 86400 s left in
  cluster mode, the stage name against the plan, x509userproxy, TransferInput,
  RequestDisk and JobFlavour in the dry run, writes submission.json whenever a
  cluster number is parsed, and records wall time and output size of the
  smoke; h4l_reco records the triggers hash, the reco-config version, gen_rows
  and the theory-weight branches per input (missing required ones are fatal:
  LHEScaleWeight/LHEPdfWeight for qqZZ, from input_requirements_v3.json),
  checks every event's file key against the plan, writes trees with
  kOverwrite (no stale AutoSave cycles), checks ROOT errors after closing each
  input and after the output validation; root_io.h checks the reopening in the
  full validation; make_reco_plan.py checks the skim configuration hash and
  the payload finalization, keeps the per-input file keys and publishes the
  plan with os.link (never replacing).  Still to do before the submission: a
  valgrind check (0 invalid accesses), the physics review, plan v2, smokes of
  a data and a signal task; a reco output scan before downstream use.
- 2026-09-24T19:47Z, stage 3c tag-and-probe run "nominal" (production_v3/tnp/v2/nominal,
  kept): all 8 fit jobs finished, the SF step stopped on a null efficiency
  error.  Root cause: 8 MC fits near eps = 1 had status 0 but an invalid
  (NaN) covariance and were still flagged usable.  Fix in fit_tnp.cpp: a fit
  is usable only with a finite positive efficiency error (else the bin falls
  back to counting); run_tnp.py guards the SF.  Rerun as "nominal_r2".
- 2026-09-24T20:15Z, stage 4a physics review (independent agent) of h4l_reco: no A item;
  0 mismatches in an independent re-derivation on 15,000 ggH records
  (deviates, calibrated pT, the complete candidate list, masses, smart-cut
  charges, c_z1_closer, iso_fsr, FSR cuts, l_overlap, ghost cleaning).  B
  items FIXED in h4l_reco (config reco_ul16_v2.json, version reco_v2): the
  muon flags now use the stage-2 layout of h4l::muon_id (0 global, 1
  tracker, 3 PF, 8/9 highPt, 10 nStations > 0, ...), so that the shared ID
  helper cannot silently narrow the fake-rate denominator; jets are stored
  from 20 GeV (the 30 GeV cut downstream, for jet-energy variations; 24
  jets); the charged isolation l_iso_chg is stored (other FSR-isolation
  definitions can be rebuilt); the payload eta variables are asserted.  C
  items: 1024 candidates; per-event nlep_all and njet_all; floors below
  every threshold (m4l > 65, Z masses < 130, m_Z2 > 3); l_overlap_mask (all
  loose muons within 0.05); a deterministic tie-break of c_z1_closer; rho
  and Pileup_nTrueInt for diagnostics.  Downstream (h4l_select) items: a
  second ghost cleaning (a non-global non-PF tracker muon within dR < 0.05
  of a same-charge PF muon), the AN FSR association (photons of loose muons
  passing SIP, electron veto |eta_SC| < 1.479 or dR > 0.08) to be evaluated,
  the SS control region with m4l > 100 GeV (AN l.795), the tracker high-pT
  muon alternative to be used identically in the SF and the selection,
  signal PSWeight = [1.0] is a dummy.  Observations: WP90 noIso keeps 86.7 %
  of prompt electrons (4e/4mu SR ratio 0.34 against 0.51 in paper Table 1;
  WPL 0.55) - input to the stage-4b electron-ID optimization; the private ggH
  sample has A_fid = 0.526 (paper 0.398), a property of the sample; the
  calibration smear variance is negative at the lowest pT node in some
  regions (r clamped to 0).  Log correction: the 26,794 selected ggH events
  are 30.0 % of the 89,245 GenTable rows (27 % of the 99,000 generated).
- 2026-09-24T20:15Z, QUESTION (English) answered by USER (verbatim): "以后用中文问我，重新问".
  Rule recorded in AGENTS.md: questions to the user are written in Chinese.
- 2026-09-24T20:15Z, QUESTION (Chinese, AskUserQuestion): the fake-rate denominator, since
  the recorded Z+X lesson (no SIP in the denominator) conflicts with the
  SIP-passing fail leptons of 2P2F/3P1F/SS (Z+X underestimated by about
  P(SIP<4 | loose fake)).  USER ANSWER: "Loose + SIP<4（推荐）": denominator =
  AN loose lepton passing SIP < 4 (no ID, no isolation, not narrowed further),
  numerator = tight ID + isolation, the no-SIP version as a cross-check.
  AGENTS.md Z+X lessons updated (the superseded rule replaced).
- 2026-09-24T20:17Z, USER (verbatim): "什么是对的用什么，原来的要求、规则、代码不一定对，如果有问题你问我，不确定的问我，并且如果确定一定是有问题的就修正，不一定要按照原来的，重点是做对，并且按照html的score能得满分，而不是一致".
  Rule recorded in AGENTS.md ("Correctness over precedent").  Consequences
  taken now: (1) the FSR-subtracted isolation follows AN 3.3 exactly (photons
  attached to loose muons passing SIP < 4, subtracted from the isolation of
  loose leptons passing SIP < 4, veto dR > 0.01 for muons and |eta_SC| <
  1.479 or dR > 0.08 for electrons) in both the selection and the
  tag-and-probe, instead of the earlier approximation kept for consistency;
  (2) the AN tight muon (PF, or tracker high-pT ID above 200 GeV) in both;
  (3) the tag-and-probe is rerun with these definitions (v3) before the SFs
  are used.  Re-read of the scoreboard (ref/Slopbench H4l v0 (Copy).html,
  "How scores work"; truths not read) and prompt section 9: the reported
  efficiency SFs must state standard working points (the scorer chooses its
  truth by them) - input to the stage-4b electron-ID decision; the
  calibration nuisances are scored like POIs (honest errors from the closure
  studies are essential).
- 2026-09-24T20:19Z, stage 4a h4l_reco v2: valgrind on a 2-input DY task: 0 invalid
  reads/writes (the remaining "conditional jump" reports are all inside
  ROOT's streamer/Cling initialization).  Plan production_v3/h4l_reco/v2/plan.json
  (74 tasks, reco_v2 config, required theory branches LHEScaleWeight and
  LHEPdfWeight for qqZZ).  Smokes through the frozen wrapper (program key
  6ce22c40c5c4): reco_pseudo_data_0000 (373,489 slim events, 297,256 kept,
  13,411 with a candidate; 22 s, 546 MB), reco_VHToZZ_M125_0000 (71,418 kept,
  205,990 candidates, no overflow; 20 s, 442 MB), reco_DYJetsToLL_0003
  (largest, 15 inputs; 27 s, 533 MB).  Dry run ok (proxy 512,192 s).
  SUBMITTED cluster 1153957, 71 jobs, microcentury, 2000 MB, 2000 MB disk.
  Completion gate: 74 outputs, then scan_reco.py --version v2 --manifests v2.
- 2026-09-24T20:22Z, stage 4a h4l_reco v2 DONE: cluster 1153957, 71/71 attempt records
  "completed" within about 2 minutes; full scan (scan_reco.py --version v2
  --manifests v2): 74/74 outputs valid (plan binding, config/trigger/payload
  hashes, one program, ROOT sha256 recomputed, root_check, counts), coverage
  complete.  Events4l: pseudo_data 3,051,070; ZZTo4L 7,559,948; DY 2,069,083;
  TTBar 2,364,132; GGZZ4Mu 957,205; GGZZ4E 60,546; GGZZ2E2Mu 4,989; ggH
  73,718; VBF 75,761; VH 71,418.  Overflows: 4 TTBar events with more than 12
  loose leptons (lepton list truncated), no candidate or jet overflow.  Scan
  record production_v3/h4l_reco/v2/scan.json.
- 2026-09-24T20:22Z, h4l_select updated per the user's correctness rule and the physics
  review: the AN FSR-subtracted isolation rebuilt from l_iso/l_iso_chg with
  the photons of loose muons passing SIP (h4l::fsr_in_isolation), the AN tight
  muon (h4l::an_tight_muon: PF, or tracker high-pT above 200 GeV), a second
  ghost cleaning (tracker-only muon within dR < 0.05 of a same-charge PF
  muon), the cross cleaning against any selected muon (l_overlap_mask), the
  SS region with m4l > 100 GeV; selection_ul16_v1.json extended accordingly.
- 2026-09-24T20:29Z, USER (verbatim): "对于eff，你不要管什么MC truth，就是用data TnP的，因为MC truth只是比值，里面还有不是来自Z的muon/ele等".
  Recorded in AGENTS.md.  State of the code: run_tnp.py computes SF =
  eps_data / eps_MC with both efficiencies from the same tag-and-probe fits
  (nominal model; counting only where a fit is unusable); the MC truth
  (generator-prompt probes in 80-100 GeV) never enters an SF and is only drawn
  as a diagnostic of the fit method.
- 2026-09-24T20:29Z, stage 3c tag-and-probe nominal_r2 DONE (sf.json, 20:10Z); plots
  running (plot_tnp.py).  First look (eff_mm_id): muon ID | loose SF rises
  smoothly with pT from about 0.989 (10 GeV) to 0.997 (150 GeV), the same in
  all |eta| bins; the MC fits reproduce the MC truth (method closure).
- 2026-09-24T20:29Z, h4l_select development runs on the v2 records (AN FSR isolation, AN
  tight muon, ghost cleaning): pseudo_data_0000 (1/11 of the data): SR 58,
  SRZ4l 84, 2P2F 1,439, 3P1F 80, SS 1,456, ZL 85,736 (14 s); ZZTo4L_0000: SR
  27,297, SRZ4l 38,191, 2P2F 1,250, 3P1F 8,113, SS 2,190, ZL 37,442 (396 s);
  ggH: SR 26,796, refit failures 1.7 % (398 s).  Reviews (physics,
  production) of h4l_select running.
- 2026-09-24T20:29Z, calibration closure sloped_same (mode same) did not converge in 10
  iterations (chi2/ndf of the muon scale 3,700-15,000 per ~100, oscillating);
  finalized with the mean of the last 4 payloads; closure report at (45 GeV,
  1.2): muon scale -0.922 % vs injected -1.000 % (pull 0.94 with the inflated
  error), muon smear variance 4.0e-5 vs 6.4e-5 (pull -1.05), electron scale
  -1.411 % vs -1.400 %, electron smear variance 1.2e-4 vs 1.0e-4.  Diagnosis
  in progress: the smear model is additive in r^2, while a smear linear in
  pT and |eta| has an eta x pT cross term in r^2 that c(|eta|) + d_R(pT)
  cannot represent; the same mode has no statistical noise between data and
  templates, so any mismatch dominates chi2.  Started: closure_sloped_halves
  (new scenario in closure_ul16_v2.json, independent halves) and
  closure_uniform_halves (run_calibration.py got --closure-config).
- 2026-09-24T20:30Z, USER (verbatim): "不哟追求和truth一样" (read as 不要追求和truth一样): agreement of the tag-and-probe with the MC truth is not a goal; no fit or definition is tuned toward it.  Recorded in AGENTS.md.
- 2026-09-24T20:32Z, USER (verbatim): "还有我不明白，为什么Research/codex/jfc/analyses/ref_h4l/production_v3/tnp/v2/nominal_r2/plots/gallery_ee_sip_data.pdf等里面，data的误差都这么大啊？不应该比如100误差就是大概10吗？".
  Root cause: fit_tnp.cpp wrote the bin variances into the curve fields
  pass_error/fail_error, so the galleries drew variances as error bars (and
  the gallery pulls were too small).  The fits used the variances correctly
  (Poisson/sandwich likelihood), so efficiencies and SFs were unaffected.
  Fixed: the curve errors are sqrt(variance); plot_tnp.py now stops if a data
  bin error is not sqrt(N).  The nominal_r2 galleries are superseded by the
  galleries of the next run.  SF(full) of nominal_r2 at (45 GeV, 1.2): muon
  0.9940 +- 0.0001 (stat) +- 0.0001 (fit model), electron 0.9872 +- 0.0006 +-
  0.0016 (with the approximate isolation and WP90; superseded by tnp v3).
- 2026-09-24T20:44Z, USER (verbatim): "对于TnP，参考Research/codex/slopbench_code_fork/results/h4l_seeds/efficiency_tnp_an里面是怎么分bin的，以及怎么选模型的，注意TnP想要一次性全部拟合好几乎不可能，你需要每次都对于拟合的结果完整检查（看起来不行的），然后一点点修正，有问题的越来越少，不能按照pull或者chi2".
  Read (within the authorized scope): tnp_validation.tex sections "Tag and
  probe", "The fit", "Quality control" and "The fit against the counting",
  and the binning/fit sections of benchmark/configs/h4l_tnp.json (seed
  results and closures not read).  Adopted in tag-and-probe v3
  (tnp_ul16_v3.json, tnp_histograms/fit_tnp rewritten, run_tnp.py
  restructured): the fine binning of the reference (for millions of
  distinct probes: muon pT 5-10-15-20-22.5-...-45-50-60-80-inf, |eta|
  0-0.2-...-1.6-1.85-2.1-2.4; electron pT 7-12-17-20-22.5-...-80-inf,
  |eta_SC| 0-0.3-0.6-0.9-1.2-1.479-1.7-1.9-2.1-2.3-2.5), dR(tag, probe) >
  0.2, one probe per tag and step (closest to m_Z, the step's denominator
  applied before the choice), the chain id | loose -> sip | id -> iso | id,
  sip with the AN FSR isolation and the AN tight muon, starting widths from
  the reference width model W = (1 + 0.76 <|eta|>)(1 + 1.0 <pT>/100) scaled
  from the inclusive fit, width bounds x3 clipped to 0.3-8 GeV with the
  bound relaxation, the HESSE retry with the resolution fixed, convergence on
  the EDM (0.01), rejection of multi-maximum totals and of backgrounds peaking
  within 12 GeV of m_Z, per-bin overrides for the manual fixes, the report
  point by bilinear interpolation between bin centres.  The fit stage is a
  separate label (fits_i1, fits_i2, ...) per inspection round, reusing the
  histograms; the gallery is written as PNG pages (16 fits per page) for the
  complete visual inspection.  Run started: production_v3/tnp/v2/v3 (label
  fits_i1).
- 2026-09-24T20:57Z, USER (verbatim): "也可以参考Research/codex/efficiency_recover/v2_seed18_brief、Research/codex/efficiency_recover/v2_seed25_detailed的TnP设置，但是数据不同，只是借鉴方法、pdf、规律等，不是照抄参数".
  Reading those references for methods (line shapes, sharing rules, pairing,
  binning logic, inspection practice, parametric SF), never their numbers.
- 2026-09-24T20:57Z, USER (verbatim): "你最后基本全都完成之后，要找3个subagent来review，参考paper、AN以及其他要求、plan、log等信息，检查是否有问题。直到确保没问题了，你问我真值的位置，我发给你，然后你撰写AN".
  Recorded in AGENTS.md ("Final reviews before unblinding").
- 2026-09-24T21:01Z, stage 4a physics review of h4l_select (independent agent): A1 the
  refit error (piecewise-linear log density: HESSE sees no curvature; 59 % of
  refits keep the unrefitted error), A2 control-region jets not cleaned from
  the failing candidate leptons (93 % of 2P2F rows have such a jet), A3 stale
  MELA fields in SS rows, A4 no production discriminants for SS rows, A5 the
  Z+1 lepton tree cannot support the SS-method fake rates and conversion
  correction (7 GeV window applied in the program, no tag kinematics); B1
  MELA-normalized D_bkg^kin squeezed near 1 (qqZZ median 0.97) and D_1jet
  above 0.5 for 88.5 % of ggH 1-jet events; B2 electron WP (WPL gives 4e/4mu
  0.55, WP90 0.34); B3 no Z+X control regions for the Z->4l region; C items
  (SS m4l cut in the program, smear response at r = 0, SIP-passing loose
  count, jet cleaning from all selected FSR photons, error closure, signal in
  3P1F, asserts).  Verified correct: SR cuts, candidate choice, lepton roles,
  OS control regions, Z+1 lepton definition; core m4l resolution 1.50 / 2.61 /
  2.00 GeV (4mu/4e/2e2mu) improving to 1.30 / 2.19 / 1.68 with the refit.
- 2026-09-24T21:01Z, QUESTIONS (Chinese) and USER ANSWERS: (1) discriminant constants:
  "我们的目的不是啥都要和paper一样，而是方法参考，具体筛选什么当然要优化啊，怎么可能照搬" -
  the selection, discriminant definitions and working points are optimized
  on MC, the paper is a method reference; (2) electron ID: "WPL 起点，MC
  上比较（推荐）".  AGENTS.md updated (the superseded discriminant rule
  replaced).  Consequence for the plan: stage 4b optimization comes before the
  final tag-and-probe; h4l_select starts from WPL electrons.
- 2026-09-24T21:05Z, h4l_select FIXES after the physics review (config sel_an_v2): A1
  the Z1 line shape is interpolated with a Catmull-Rom cubic in log density
  and the refit covariance is the inverse of the numerical Hessian of the
  likelihood (steps of 0.3 sigma per lepton); a refit without a positive
  Hessian falls back to the unrefitted values (status recorded); A2 jets are
  cleaned per candidate (selected leptons, the candidate's leptons and the
  selected FSR photons, i.e. those of loose muons passing SIP: C4); A3 every
  MELA field is reset per candidate; A4 SS candidates get MELA with the
  second Z2 lepton's charge flipped (mela_charge_flipped = 1; production MEs
  valid, D_bkg^kin approximate); A5 the Z+l tree keeps Z1 in 40-120 GeV
  (windows downstream) and stores m3l, the probe-tag dR, the probe flags and
  the SIP-passing loose count (C3); B3 the control regions are produced for
  both signal regions (region field); C1 SS rows from m4l > 70 GeV (100
  downstream); C2 the smear response is d m4l / d r (absolute); C7 the first
  region is asserted to be SR; the candidate leptons' ID flag words are stored
  (tighter working points downstream); electrons mvaFall17V2noIso WPL (the
  stage-4b starting point).  Development reruns started.
- 2026-09-24T21:25Z, T&P FIT MODEL v4 (after the efficiency_recover lessons, methods only):
  the probe pT binning sculpts the Breit-Wigner tails of the pass and fail
  spectra differently, so the nominal signal becomes a stand-alone
  double-sided Crystal Ball in m - m_Z (zpeak_model.h Kernel::StandaloneDCB)
  with its own shift and width per category and the four tail parameters
  shared; both categories share one resolution when either has fewer than
  300 raw entries (the efficiency_tnp_an objection to free fail shapes at
  high efficiency); alternatives BW (x) one-sided CB (alt_signal) and a
  falling Bernstein-2 background (alt_background); two starting points
  (inclusive tails, default tails), the converged fit with the lower NLL
  kept; storage bins of 0.5 GeV usable directly (the narrowest of 0.5-3 GeV
  with >= 40 entries per bin); held parameters by effective probes (< 1000:
  n; < 300: alpha; < 100: shifts and widths); a total curve with more than
  one maximum per category is rejected; a yield is "at limit" within half an
  event of a bound (the relative criterion flagged a 30 +- 7 background
  yield); the diagnostic chi2 floors the variance at expectation x mean
  weight (weighted MC bins gave chi2/ndf ~ 1e6; never an acceptance
  criterion).  Configuration analysis_v3/tnp/config/tnp_ul16_v4.json
  (tnp_v4): identical to v3 except the electron id step
  mvaFall17V2noIso_WPL (the stage-4b starting point) and the fit block.
  Gallery pages (plot_tnp.py) now show 8 fits per page, pass and fail side by
  side on a linear scale with pulls and the alternative/counting
  efficiencies.  Test on 9 muon isolation bins (production_v3/tmp/tnp_dev/fit4):
  all data fits ok, efficiencies within 0.001 of the alternatives in the
  high-statistics bins; a radiative shoulder at 72-78 GeV in the isolation
  fail spectrum of pT 27.5-30 GeV probes (data and MC alike: an in-cone FSR
  photon not recovered moves the Jacobian peak down and fails the
  isolation) is not described by the shared left tail -- to be counted on the
  full galleries.  The v3 fits_i1 run (the superseded shared BW (x) CB
  model) was stopped at about one third; its logs are kept.  Started
  run_tnp.py --extract v2 --calibration v3/nominal_r3 --run v4 --fit-label
  fits_i1 (tmux tnp_v4).
- 2026-09-24T21:25Z, h4l_select DEVELOPMENT RERUNS (sel_an_v2, production_v3/tmp/reco_dev/run3_*):
  data shard 0: SR 77, SRZ4l 105, 2P2F 1519, 3P1F 151, SS 1660 (m4l > 70),
  Z+l 124,457 rows (Z1 40-120 GeV), refit failures 4 of 7251; ZZTo4L_0000:
  SR 32,839, refit failures 27 of 95,301; ggH_0000: SR 32,335 (WPL; 26,796
  with WP90), refit failures 41 of 77,496.  Per-event error closure on ggH
  prompt SR events (core Gaussian, +-2 sigma iterative): pull widths
  1.10 / 1.21 / 1.14 (4mu / 4e / 2e2mu) unrefitted and 1.10 / 1.29 / 1.20
  refitted, i.e. the lambda-corrected errors underestimate the H -> 4l
  resolution by 10-30 % in MC (to be handled in the signal model, where the
  DCB width is parameterized against the per-event error from MC).  Control
  regions: data 2P2F / 3P1F / SS mean jet multiplicity 0.98 / 0.83 / 0.99
  (SR 0.46); ZZ MC 3P1F 0.455 against SR 0.432 (the per-candidate cleaning
  works); SS rows have finite D_bkg^kin in [0, 1] with the charge flip.
  Timing: data 0.28 ms per row, ZZ 2.8 ms, ggH 9 ms (MELA for every
  candidate); the production review measures GGZZ4Mu_0001 (278k rows) at
  about 64 min.
- 2026-09-24T21:41Z, stage 4a PRODUCTION REVIEW of h4l_select (independent agent; it read the source
  and binary from before the 21:05Z physics fixes): A1 stale MELA fields in SS
  rows (already fixed at 21:05Z: every MELA field is reset and SS rows get
  MELA with the flipped charge; made robust now by resetting the whole row),
  A2 incomplete MELA provenance (the dictionary LinkDef_out_rdict.pcm and the
  runtime data files were not verified: MCFM/JHUGen inputs, PDF grids,
  HiggsTotalWidth_YR3, inputs_4mu, 17 pAvgSmooth constant files, the two mJJ
  resolution files), A3 a MELA constant file that fails to load leaves K = 1
  silently (valid-looking but wrong D), A4 no configuration hashes in the
  report, A5 GGZZ4Mu_0001 (278,553 rows) takes 3,715 s (would be killed at the
  microcentury limit); B1 the Fortran code opens the shared MELA data files
  O_RDWR through the links MELA makes in the working directory, B2 SR and
  SRZ4l recompute MELA and the refit for the same candidate (factor ~1.9 on
  GGZZ4Mu), B3 ROOT errors during the MELA initialization were discarded, B4
  no loader pre-flight on workers; C: compiler pinned, counters initialized,
  truncation counters, line-shape asserts, plan-builder asserts, ROOT version
  and host, per-event file-key check.
- 2026-09-24T21:41Z, h4l_select FIXES (production review): the row starts from its defaults in
  every fill; a candidate chosen in several regions is computed once per event
  (row cache; only region, cr_type and ncand_region differ); the plan lists
  every MELA library, dictionary and data file (32 files) and the eight data
  files MELA opens for update, which the task copies into its empty working
  directory (it refuses a directory that already holds one of the names, so
  it never writes through a link into the installation) and verifies; after
  the initialization every task repeats a fixed six-candidate MELA self-test
  (4mu, 4e, 2e2mu without jets, a 1-jet, a VBF-like and a VH-like 2-jet event)
  against the reference values that make_select_plan.py computes with the
  program (--mela-selftest, fresh directory under production_v3/tmp/) and
  fails on a relative difference above 1e-6; the report (schema
  h4l_v3_select_report/2) carries fnv1a64 hashes of the embedded selection
  configuration, calibration payload, lambda payload, refit inputs and MELA
  block, the ROOT version, the host, the initialization and input-open ROOT
  errors (recorded, not fatal), the lepton/jet/candidate truncation counters
  and the decay-only MELA failures; every event's file key must be in the
  task's source keys; the Z1 line shape must have uniform bins; run_task.sh
  fails with 66 when ldd reports an unresolved library; the reconstruction
  Makefile pins /usr/bin/g++.  scan_select.py (new) binds each report to its
  plan task (schema, identity, the five hashes, inputs, source keys, MELA
  files and working copies, a passed self-test), recomputes the ROOT sha256,
  runs root_check, checks the tree counts against the counters (refit
  counters equal the filled rows, SR <= SRZ4l) and the per-sample coverage of
  the manifest files, and lists the outputs with logged ROOT errors.  Plan
  production_v3/h4l_select/v1/plan.json (74 tasks) written; local runs of
  three plan tasks started (production_v3/tmp/reco_dev/run4_*).
- 2026-09-24T22:03Z, h4l_select v2 SUBMITTED: plan production_v3/h4l_select/v2/plan.json (74 tasks;
  v1 is superseded unsubmitted: its self-test reference was computed before the
  canonical configuration hash, identical values), program frozen as
  program_inputs/h4l_select_v2/f561b55a8f44.  Regression: the fixed program
  reproduces the pre-fix outputs bit for bit on pseudo_data_0000, ggH_0000
  and ZZTo4L_0000 (all SR/SRZ4l/CR/ZL branches), 1.6-1.8x faster (ggH 380 s
  instead of 672 s, ZZ 410 s instead of 669 s).  The configuration hashes
  use a canonical text with %.17g floats (the shortest round-trip form is
  ambiguous at ties, e.g. 66.587677001953125, where nlohmann and Python
  differ); all five hashes verified against Python.  Smokes through the
  wrapper (PATH=/usr/bin:/bin): sel_pseudo_data_0000 49 s, 505 MB, 8.3 MB;
  sel_GluGluToHToZZ_M125_0000 5 min 30 s, 445 MB, 22 MB (both published).
  Dry run (72 jobs, Iwd, logs, proxy, TransferInput checked), then cluster
  1153978 (72 jobs, longlunch, 2000 MB memory, 2 GB disk); proxy 506,073 s.
- 2026-09-24T22:29Z, h4l_select v2 FAILURE and v3: sel_TTBar_0013 failed deterministically (three
  attempts, "lepton outside every lambda region"): one muon (pT 12.9 GeV,
  |eta| 0.45) of reco_TTBar_0013 has a NaN momentum error and a NaN SIP.  A
  scan of every field of all 74 event records (production_v3/tmp/reco_dev/
  scan_nonfinite.py) finds no other non-finite value.  Root-cause rule: a
  lepton without finite kinematics, a finite positive momentum error,
  isolation and SIP has no usable measurement and is not a loose lepton
  (counter lepton_unmeasured).  The program change requires a new plan:
  production_v3/h4l_select/v3/plan.json (74 tasks), frozen program
  program_inputs/h4l_select_v3; smokes sel_TTBar_0013 (25 s, 498 MB,
  lepton_unmeasured 1) and sel_pseudo_data_0000 (37 s, 506 MB); dry run;
  cluster 1153979 (72 jobs, longlunch).  The v2 run (cluster 1153978, 71 of
  74 outputs at the time) is superseded and never consumed.
- 2026-09-24T22:53Z, T&P v4 fits_i1 DONE (production_v3/tnp/v2/v4/fits_i1; histograms 21:22-21:24Z,
  fits 21:24-22:31Z, 8 jobs of 342-376 fits): status counts data/MC
  (ok/constrained/rejected/counting-only): mm id 135/48/5/0 and 131/41/1/15,
  mm sip 140/43/3/2 and 145/39/2/2, mm iso 146/38/1/3 and 150/35/0/3, mm full
  150/33/5/0 and 152/35/1/0, ee id 122/46/1/2 and 127/41/2/1, ee sip
  118/48/5/0 and 120/47/4/0, ee iso 104/43/4/20 and 105/42/3/21, ee full
  137/32/2/0 and 138/30/3/0.  Main reasons for "constrained": n_L (or n_R)
  at a bound, a background yield at zero, the fail exponential slope at its
  bound (fail_shape1), the fail resolution at its bound (ee iso); a few
  unconverged fits.  SF(full) at (45 GeV, 1.2): muon 0.99393 +- 0.00013,
  electron (WPL) 1.00026 +- 0.00036 (statistical and fit-model errors, first
  round, before the inspection).  Plots and 360 gallery pages:
  production_v3/tnp/v2/v4/fits_i1/plots/ (gallery/ 8 fits per page).  The
  complete visual inspection is split over eight independent inspectors
  (one per flavour and step, both roles, every page and fit, criteria: peak
  region described in pass and fail, fail signal not misattributed,
  physical background, no failed fit; never chi2 alone); their lists feed
  the fixes of fits_i2.
- 2026-09-24T22:58Z, CALIBRATION HALVES CLOSURES finished without converging (10 iterations;
  finalize_calibration.py: mean of iterations 7-10): the joint least squares
  is statistically inconsistent in halves mode, unlike the nominal run
  (nominal_r3 muon scale chi2/ndf 1.0-1.4 in every iteration): uniform
  halves muon scale chi2/ndf 27 at iteration 1 growing to 127-380, dominated
  by family-B categories with pulls of 10-47; payload oscillations up to
  0.012-0.032 in the scale (iteration rms).  Closure numbers at (45 GeV,
  1.2) with the inflated errors: uniform halves muon scale -0.01463 +-
  0.00237 (injected -0.015), smear 0.0058 (0.010; variance pull -1.85);
  electron scale -0.01253 +- 0.00309 (-0.012), smear 0.0131 (0.012);
  sloped halves muon scale -0.00782 +- 0.00265 (-0.010), smear 0.0106
  (0.008); electron scale -0.01367 +- 0.00202 (-0.014), smear 0.0109 (0.010).
  The halves closure cannot validate the nominal precision in this state;
  an independent investigation of the root cause (response, B-category
  errors, outlier handling, compositions) is running (scratch area
  production_v3/tmp/closure_debug/).  closure_report.py takes
  --closure-config (sha256 checked against the run record); the reports are
  in <run>/closure/.
- 2026-09-24T23:28Z, T&P fits_i1 INSPECTION (eight independent inspectors, every page and fit of
  the 360 gallery pages; reports summarized): (1) muon id: no visibly bad
  fit; (2) in all steps, errors collapse (0.0000-0.0002, below the binomial
  floor) or the fit is "rejected" whenever a background yield or slope sits
  at its bound (b_pass at pT 35-60 GeV; HESSE with parameters at bounds);
  (3) a low-mass hump in the FAIL spectra at 60-80 GeV whose position rises
  with the probe pT (~65 GeV at 20-22.5, ~71 at 25-27.5, ~77 at 30-32.5,
  merging into the Z at >= 35 GeV), identical in data and MC, strongest for
  electrons in the iso and full steps (in the 1.479-1.7 transition taller
  than the Z peak) and for muons at |eta| 1.2-1.6 (iso, full, weaker in
  sip); the shared-tail DSCB + exponential cannot follow it: the exponential
  is dragged through it or the fail signal drifts to its shift/width limits
  and merges hump and peak (eps biased by 1-5 %, same sign in data and MC);
  (4) electron fail backgrounds have a kinematic turn-on (combinatorial
  pairs peak near 2 sqrt(pT_tag pT_probe)) that a falling exponential cannot
  describe; (5) the |eta_SC| 1.479-1.7 column at pT 20-40 fails in every
  electron step; (6) cosmetic pull waves in high-statistics pass spectra;
  (7) bins with fewer than 50 fail entries are counting-only (MC mm id pT
  > 60; ee iso 20 bins).  Check in the DY MC (electrons, pT 25-27.5 GeV,
  |eta_SC| < 0.3, passing ID and SIP, failing isolation): 733 truth-matched
  prompt probes against 12 non-prompt, the prompt ones forming the hump
  (65-80 GeV) and the Z peak: the hump is Z -> ll gamma with the FSR photon
  inside the isolation cone (iso > 0.35 <=> m < m_Z/sqrt(1.35) = 78 GeV;
  no FSR photons for electrons in NanoAOD), i.e. genuine fail signal (the
  H -> 4l selection uses the same isolation).
- 2026-09-24T23:28Z, QUESTION (Chinese) and USER ANSWER: T&P signal model for the fail (and pass)
  spectra -> "MC模板⊗高斯（推荐）": the signal is the template of the
  truth-matched prompt probes of the DY MC per bin and category (smoothed),
  convolved with a Gaussian (shift and width per category); background
  CMSShape (erf x exp, EGM parameter ranges); alternatives: the analytic
  DSCB signal and a Bernstein background (fit-model systematic); only the MC
  shape is used, never an MC truth efficiency.
- 2026-09-24T23:53Z, h4l_select v3 COMPLETE and SCANNED: cluster 1153979, 72 attempts completed
  (no failure, no hold; GGZZ4Mu tasks 26-60 min on the workers), plus the two
  smokes; scan_select.py: 74 outputs valid (identity, five configuration
  hashes, MELA files, working copies and self-test, inputs, ROOT sha256,
  root_check, counters), coverage of every manifest complete, no ROOT error
  logged; production_v3/h4l_select/v3/scan.json.  Selected events (SR / SRZ4l):
  pseudo-data 838 / 1092; ZZTo4L 1,050,539 / 1,448,678; GGZZ4Mu 597,793 /
  619,197; GGZZ4E 27,641 / 28,265; GGZZ2E2Mu 2,792 / 2,900; ggH 32,335 /
  33,392; VBF 33,205 / 34,400; VH 30,587 / 31,539; DY 33 / 40; TTBar 111 /
  118.  analysis_v3/common/python/h4l_select_io.py (new): scan-based readers
  with the genWeight normalization.  Note: the calibration programs were
  rebuilt at 23:36Z because zpeak_model.h and binned_fit.h gained the
  template kernel, CMSShape, Bernstein-3 and optional MINOS (existing code
  paths unchanged; nominal_r3 stays valid with its frozen programs).
- 2026-09-24T23:53Z, T&P v5 (user model decision): tnp_histograms fills the signal templates
  tpass_/tfail_<f>_<step>_<bin | pt<row> | all> from the prompt-prompt pairs
  of the DY MC (both legs gen_flav 1, the per-(tag, step) choice closest to
  m_Z among them, 50-130 GeV at 0.25 GeV); fit_tnp: nominal = template
  (bin, else pT row, else inclusive below 500 effective entries; smoothed
  0.5 GeV) (x) Gaussian per category + CMSShape (alpha 30-85, beta
  0.005-0.2, gamma 0-1); alt_signal = stand-alone DSCB with the fail DSCB
  mixed with a Gaussian low-mass component; alt_background = template +
  Bernstein-3; efficiency error from MINOS (sandwich-scaled for weighted
  MC) with the binomial floor; Minuit status 1 accepted; background shape
  limits flagged only where the background exceeds 2 % of the category;
  narrow background maxima near m_Z rejected; categories below 20 entries
  counting only, below 100 with the inclusive background shape; curves of
  every model stored; plot_tnp.py adds the alternatives gallery (fits whose
  alternative deviates by more than max(3 sigma, 0.003)); run_tnp.py keeps
  alternatives excluded by the inspection (exclude_alternatives, with notes)
  out of the fit-model systematic.  Development check on electron and muon
  isolation bins with the FSR hump (production_v3/tmp/tnp_dev/fit5/): the
  template model follows the fail hump and the Z peak with flat pulls
  (e.g. ee iso pT 25-27.5 GeV, |eta_SC| 0-0.3: eps 0.9573 data, 0.9555 MC;
  v4 0.9819/0.9816); configuration tnp_ul16_v5.json; histograms v5
  (23:38-23:40Z); fits_i1 started 23:52Z (tmux tnp_v5).
- 2026-09-25T00:05Z, STAGE 4b SCAN i1 (optimize_selection.py, MC only; production_v3/optimization/v3/scan_i1/):
  Fisher information of a binned model per final state in m4l (refitted,
  105-140 GeV, 0.25 GeV) x D_bkg^kin (signal quintiles) x relative mass
  error (signal tertiles), mu = 1, m_H = 125 GeV, 20 fb-1; signal ggH + VBF
  + VH (0.3 GeV kernel), qqZZ + ggZZ (2 GeV kernel), Z+X from the fake-rate
  method on the DY and TTbar MC (loose + SIP denominator, OS formula);
  reference (AN values, WPL): sigma_mu 0.177, sigma_mH 0.294 GeV; yields in
  105-140 GeV: signal 14.5 / 7.8 / 21.3, ZZ 14.4 / 6.4 / 17.2, Z+X 3.6 /
  3.1 / 9.3 (4mu / 4e / 2e2mu).  J (reference 2.000): electron noIso WPL
  2.000, WP90 2.251, WP80 2.666, Iso WPL 2.009, WP90 2.205, WP80 2.597;
  isolation 0.15 / 0.20 / 0.25 / 0.30 / 0.35: 2.393 / 2.185 / 2.083 / 2.018
  / 2.000; SIP 2.5 / 3 / 3.5 / 4: 2.323 / 2.144 / 2.039 / 2.000; Z2 minimum
  4 / 6 / 8 / 10 / 12 / 14: 2.058 / 2.024 / 1.992 / 1.970 / 2.000 / 2.034;
  Z1 minimum 40 / 45 / 50: 2.000 / 2.002 / 2.006; leading pT 20 / 22 / 25:
  2.000 / 2.003 / 2.015; subleading 10 / 12 / 15: 2.000 / 2.002 / 2.015.
- 2026-09-25T00:05Z, QUESTIONS (Chinese) and USER ANSWERS: Z2 minimum -> "采用 10 GeV（推荐）";
  isolation / SIP (optimum at the edge of the testable range) -> "保持 0.35 / 4（推荐）".
  Selection configuration selection_ul16_v2.json (sel_v3): SR Z2 > 10 GeV
  (the smart cut with the same threshold), everything else as sel_an_v2;
  make_select_plan.py points to it.  Plan production_v3/h4l_select/v4 (74
  tasks; the program is the v3 one), smoke sel_pseudo_data_0000 (37 s; SR 80
  instead of 77, SRZ4l 105), dry run, cluster 1153983 (73 jobs, longlunch).
  h4l_select v3 stays the reference of the stage-4b scans (Z2 > 12 GeV SR).
- 2026-09-25T00:26Z, CALIBRATION CLOSURE INVESTIGATION (independent agent; scratch
  production_v3/tmp/closure_debug/): root cause of the halves-closure divergence in
  template_fit.cpp: weighted data (the closure data role, and the response passes of EVERY run,
  nominal_r3 included) used a scaled Poisson with a per-group scale n_b/V_b estimated from the
  observed data; a sparse group holding one low-weight ZZ/ggZZ pair (weights spanning five orders
  of magnitude) gets s = 1/w ~ 1e3-1e4 and dominates the likelihood (e.g. B_mm_10_19: fitted
  normalization 3.3x the data, ln k -0.0109 +- 0.00012 instead of -0.0028 +- 0.0014); 143/229
  family-B muon fits of the uniform closure broken; nominal main-pass fits (unweighted data)
  unaffected, but the nominal responses corrupted (baseline |ln k/sigma| median 5.6, max 558).
  Secondary: (D2) likelihood groups sized on kernel-smoothed contents (each event counted ~9
  times) and BB-lite variances ignoring the correlated shares of one event -> category errors
  too small by 1.3-1.4 (MC statistics effectively ignored); (D3) the response baseline fitted
  its own template events (D at its lower limit, voiding the category's responses).  Validation
  of the fix on the existing histograms: one corrected closure step gives scale chi2/ndf 1.18
  and smear 1.62 without outliers, u(45, 1.2) = -0.01491 +- 0.00011 (injected -0.0151), v =
  9.4e-5 +- 1.0e-5 (1.0e-4).  Impact on nominal_r3: scale at 30-60 GeV unchanged (< 1e-4 in u),
  soft/hard pT tails up to 8e-4, low-pT smear strongly changed, errors ~1.3x too small.
- 2026-09-25T00:26Z, CALIBRATION FIX (applied, confirmed bugs): template_fit.cpp: one
  effective-count scale per window for weighted data, likelihood groups >= min_group_width (1.2
  GeV) with MC effective entries of the event positions (no kernel), BB-lite variance with the
  shares of one event summed per group, Pearson chi2 in effective counts; run_calibration.py:
  common decorrelating smear r0 = 0.008 in every response variant, exclusion of categories whose
  normalization misses the data sum by more than 2 %, --config option.  Configuration
  calibration_ul16_v4.json (calib_v4); originals kept in production_v3/tmp/closure_debug/
  (template_fit_before_fix.cpp, run_calibration_before_fix.py).  Nominal calibration rerun
  started: run_calibration.py --version v4 --extract v2 --run nominal (tmux calib_v4).
  Downstream consequence: every product that applies the calibration payload is redone with the
  v4 payload (lambda, tag-and-probe histograms, h4l_select recalibrating from the stored raw pT
  and smearing deviates, signal model); the reported scale/smear come from v4.
- 2026-09-25T01:02Z, STAGE 5-7 CODE (development on h4l_select v3, v4 T&P SFs; not yet the final inputs):
  * analysis_v3/backgrounds/scripts/zx_estimate.py: Z+X OS method (AN Eq. 19; fake rates from the
    Z + 1 loose (+ SIP) lepton rows, |m_Z1 - m_Z| < 7 GeV, MET < 25 GeV, prompt ZZ subtracted with the
    ZZTo4L and ggZZ MC; no WZ sample: substitution) and SS method (AN Eq. 20; SS-space fake rates,
    electron conversion correction by the linear fake-rate vs mean-missing-hits relation of four Z + e
    samples evaluated at the SS control rows, (OS/SS)_MC from the DY and TTbar MC, SS rows with m4l > 100
    GeV extrapolated to m4l > 70 GeV with their own f3 f4 m4l distribution), bootstrap uncertainties,
    inverse-variance combination with the envelope uncertainty (AN 7.2.4), Landau shapes fitted to the
    signed OS distribution.  Development numbers (v3, 20 fb-1, m4l > 70 GeV): 4mu OS 4.3 +- 1.7 / SS
    8.5 (after the extrapolation fix), 4e OS 10.9 +- 1.6 / SS 11.7 +- 0.8, 2e2mu OS 26.5 +- 1.8 / SS
    21.8 +- 1.1; the 4mu OS estimate is low because the observed 3P1F (75 non-ZZ events) is below the
    2P2F extrapolation (111): the 3P1F term corrects the fake-rate composition; closure test on MC still
    to do.
  * analysis_v3/signal_model/scripts/yr4_xsbr.py: the YR4 spreadsheet (sha256 verified against the
    seed plan) read with the standard library; the mapped columns cover only 120-130 GeV (ggF 48.58 pb,
    VBF 3.782, WH 1.373, ZH 0.884, BR4l 1.24e-4 at 125 GeV); outside, a log-quadratic extrapolation of
    sigma_eff (documented).  production_v3/signal_model/yr4/v1/yr4.json.
  * analysis_v3/signal_model/scripts/signal_model.py: DCB per final state fitted to ggH + VBF (refit
    and not), the per-event-width DCB of the 3D fit (width = s x dm; s = 1.085 / 1.024 / 1.038 for
    4mu / 4e / 2e2mu, i.e. the lambda-corrected MC errors underestimate the resolution by 2-9 %), VH
    resonant DCB + non-resonant Landau (f_res ~ 0.95), rest-frame morphing (mean and width x k, the
    acceptance ratio from re-applying the kinematic selection to the scaled events; the shift morphing
    recorded for the systematic).  v3 development: DCB(refit) mean/width 124.90/1.25 (4mu), 124.55/1.75
    (4e), 124.78/1.44 GeV (2e2mu); signal yields in 105-140 GeV 13.99 / 7.85 / 20.95.
  * analysis_v3/common/python/{h4l_categories.py, h4l_sf.py, h4l_shapes.py}: the AN 6.1 categories
    (configuration categories_ul16_v1.json with the AN MELA-only working points), the per-lepton SF
    lookup (chain product; fallbacks recorded) and the event SF with per-flavour variations, the DCB with
    analytic normalization (unit-tested against numerical integration), Landau, Bernstein, weighted
    fits.  With the MELA-normalization constants D_WH and D_ZH are ~0, so the AN working points leave
    VH-hadronic empty: the category constants and working points are being optimized on MC by an
    independent agent (production_v3/tmp/category_opt/).
  * analysis_v3/inference/scripts/build_model.py and h4l_likelihood.py, fit_model.py: channels = final
    state x category in 105-140 GeV; yields with SFs and their flavour variations, eff_correction
    (development: 4mu 0.968, 4e 1.000, 2e2mu 0.984, inclusive 0.9815), Bernstein-3 ZZ shapes, Z+X
    Landau, D_bkg^kin templates P(D | m4l) (signal-quantile D bins) and relative-error templates; an
    extended unbinned likelihood (1D / 2D / 3D, refit or not), mu unbounded, m_H floating in 110-140
    GeV, nuisances (lepton scale and resolution per flavour, SF per flavour, Z+X per final state,
    morphing; set (a) adds luminosity, BR, qqZZ theory, ggZZ, signal QCD/PDF); MINOS intervals, the
    statistical interval with the nuisances fixed at their best fit, explicit fixed-nuisance impacts for
    the systematic part, the local significance with m_H floating; Asimov grid and unbinned toys.
    Development Asimov (3D, refit, v3 inputs): mu 0.995 +0.188/-0.171 (stat +0.187/-0.169, syst
    0.025), m_H 124.999 +- 0.24 GeV (stat 0.23, syst 0.073), expected Z 10.2.
- 2026-09-25T01:23Z, STAGE 7-9 CODE (development): make_pyhf.py (MODEL.json: the 2D model binned at the
  fitted m_H, final states x (m4l 2.5 GeV x 5 D bins), signal normfactor mu, normsys for the SFs and
  Z+X, histosys for the lepton scale/resolution, staterror; the pyhf MLE mu is checked against the
  reported mu: development check 0.971 vs 0.949, tolerance 0.096, consistent); run_toys.py
  (multiprocessing with the fork start method: the default forkserver binds a Unix socket that EOS
  FUSE refuses); toy generation rewritten after a first test showed an m_H pull width of 2.1: the
  per-event width was taken from the median error of groups of events instead of each event's own
  error, and the lowest error bin was sampled down to zero (now: the raw 0.5-99.5 % quantile edges,
  the signal m4l from the standard DCB scaled by each event's width, regenerated outside the window);
  the 3D likelihood no longer multiplies the per-event width by k = m_H/125 (the per-event error
  already scales with the mass; double counting).  48 development toys at mu = 1, m_H = 125 GeV:
  mu pull mean 0.10 +- 0.14, width 1.00; m_H pull mean 0.16 +- 0.12, width 0.85.  A development data
  fit was run on provisional inputs only to test the machinery (old calibration, v3 selection,
  provisional SFs and categories); nothing is tuned on it.
- 2026-09-25T01:31Z, H4L_SELECT v4 COMPLETE (cluster 1153983, selection sel_v3 with SR m_Z2 > 10 GeV, event
  records calibration v3/nominal_r3): scan_select.py --version v4 --manifests v2: 74 outputs valid,
  coverage complete, no logged ROOT errors (production_v3/h4l_select/v4/scan.json).  SR candidates:
  pseudo-data 869 (838 with m_Z2 > 12), ZZTo4L 1109422, ggZZ 4mu/4e/2e2mu 602285/27793/2807, ggH 33278,
  VBF 34261, VH 31406, DY 35, TTbar 112 (unweighted entries).  This is an intermediate selection: the
  final one (v5) recalibrates from the raw pT with the calibration v4 payload and the v4 lambda.
  T&P v5 fits_i1 galleries complete; the eighth inspector (ee_full) launched.
- 2026-09-25T12:33Z, CALIBRATION v4 FINAL, LAMBDA v4, Z+X CLOSURE, T&P v6, USER DECISIONS:
  * Calibration v4/nominal: 10 iterations, worst residual 180.6 -> 45.9 -> 9.07 -> 3.09 -> 1.68 -> 1.42 ->
    1.26 -> 1.08 -> 1.11 -> 1.45 -> 1.31 sigma (a plateau at the statistical level, not converged by the 1-sigma
    rule): finalize_calibration.py --last 4 (mean of the payloads applied in iterations 7-10; largest iteration
    rms: mm scale 1.3e-4, smear 1.4e-5; ee scale 4.5e-4, smear 2.0e-5, added to the statistical uncertainties).
    plot_calibration and diagnose_calibration --final (diag_final) done.
  * lambda v4 (run_lambda --version v4 --run nominal): the first attempt (label lambda_v2) stopped: electron
    |eta| 1.0-1.2 data fit status 1 (a background yield at zero made HESSE non-positive) and |eta| 1.57-2.0 MC
    fit EDM 1.6e4 (tail exponents running to their limit).  fit_lambda.cpp: iterative freezing (up to three
    rounds) of nuisance parameters at or next to a limit (1 % of the range for shapes, half an event for
    yields; lambda never fixed) and three alternative starting points of the shared DCB tails when the tails
    are free, the converged attempt with the lowest NLL kept (tail_start recorded).  Rerun lambda_v2_r2: all 12
    regions ok (muon data/MC 1.473/1.511, 1.267/1.295, 1.058/1.083; electron e.g. |eta| 0-0.8, d < 0.03:
    1.544/1.587; |eta| 1.57-2.0: 1.445/1.477).  lambda_v2b_r2 (method variant): the fits finished, the closure
    plot crashed on a closure category without an error (null): run_lambda.py draws such a point without an
    error bar; rerun.
  * Z+X development on h4l_select v4: the AN 7.2.3.2 per-process composition estimate implemented (DY, TTbar):
    4mu about 80 %, 4e 4-19 %, 2e2mu 4-44 %; the TTbar high-pT muon "fake rate" (0.3-0.8) is dominated by prompt
    W leptons in the Z1 combinatorics, so the estimate is ill-defined with the available samples.  MC closure
    (zx_closure.py, production_v3/backgrounds/v4/closure_v2; SS with the conversion correction as in the data):
    m4l > 70 GeV OS/SR = 0.54/0.71/0.71 (4mu/4e/2e2mu), SS/SR = 0.76/0.79/0.61, summed over final states 0.68
    for both (SR MC 64.3 +- 7.9); 105-140 GeV OS/SR 0.67/0.63/1.05; 3P1F from 2P2F observed/predicted data
    0.69/1.01/1.18 and MC 0.75/1.01/1.14 (the data behave like the MC).
  * QUESTION (Chinese) -> ANSWER "MC closure 合并（推荐）": the Z+X systematic is max(|1 - r|, sigma_r) of the
    combined closure (larger of the two methods, about 32 %), relative, for both methods and every final state,
    in quadrature with the bootstrap statistical uncertainty; central values not corrected; the AN
    composition estimate and the SS conversion-correction uncertainty kept as cross-checks.  Implemented in
    zx_estimate.py (mc_closure; the closure helpers moved there and shared with zx_closure.py); AGENTS.md
    updated.  Development run dev_zx3 (v4) started.
  * QUESTION (Chinese) -> ANSWER "采用 MC 优化方案（推荐）": categories_ul16_v2 (paper-convention D' = 1/(1 + c
    P_bkg/P_sig) from the stored raw MELA probabilities, c_2jet 0.1, c_1jet 0.0014, c_WH 0.16, c_ZH 0.4, working
    points 0.5, VH-MET MET > 90 GeV; the stage-4b optimization agent's proposal: J = sigma(mu_VBF)^2 +
    sigma(mu_VH)^2 relative 2.00 -> 1.70 +- 0.015, sigma(mu_VH) -15 %, VH-hadronic 0.81 signal events, purity
    33 %).  h4l_categories.py recomputes the discriminants when the configuration has a "discriminants" block
    (validated: 0 mismatches against the optimization's reference implementation on ggH, VBF, VH, ZZ and the
    data rows of v4 in 105-140 GeV); signal_model.py defaults to v2; AGENTS.md updated.
  * T&P v5 inspection: the eight inspector agents stopped at the usage limit; their partial findings (mm_id,
    mm_sip complete notes): coherent peak pulls of every pass fit above ~70k entries (the fixed 0.5 GeV template
    smoothing broadened the peak), pT-row templates too wide for central-|eta| fail peaks, template wiggles in
    sparse fail templates, CMSShape humps under the Z (alpha, beta at their upper bounds), spurious narrow
    low-mass Gaussians of the DSCB alternative.
  * T&P v6 (configuration tnp_ul16_v6.json, histograms with the calibration v4 payload): adaptive template
    smoothing (Abramson sample-point kernel, pilot bandwidth by Silverman's rule of the Z-peak entries within
    +-10 GeV of m_Z; a two-region variant was tried and rejected: the bandwidth step distorted the templates);
    template fallback bin -> |eta| column with neighbouring pT rows within a factor 2 in pT -> pT row -> all;
    EGM CMSShape ranges (alpha 50-80, beta 0.01-0.06, gamma 0.005-1); DSCB-alternative low-mass Gaussian mean
    65-88, width 3-12 GeV.  Development fits on the v5 histograms (mm_sip, mm_id, ee_iso problem bins): pass
    pulls without the coherent pattern, fail peaks followed.  fits_i1 done: SF(full) (45, 1.2) muon 0.99397 +-
    0.00015, electron 1.00038 +- 0.00038.  The plots crashed on one fit without an error (mm_id MC bin 128:
    MINOS flagged, HESSE NaN): plot_tnp.py compares such fits with the fixed 0.003 only; fit_tnp.cpp (next round)
    falls back MINOS -> HESSE -> a flagged MINOS interval that brackets the minimum (error_source recorded).
- 2026-09-25T12:44Z, H4L_SELECT v5 SUBMITTED: plan make_select_plan.py --reco v2 --version v5 --calibration v4/nominal
  --lambda-run v4/nominal/lambda_v2_r2 --refit v3 (74 tasks; calibration block same_as_event_records = false: every
  lepton recalibrated from its raw pT with the v4 payload); local smokes passed (sel_pseudo_data_0000 43 s, 506 MB;
  sel_GluGluToHToZZ_M125_0000 441 MB; smoke output vs v4: SR rows 79 vs 80, median pT/pT_raw muons 1.02001 vs
  1.02007, electrons 1.02021 vs 1.02006, the report records the v4 payload); dry run ok (Iwd, logs); 72 jobs
  submitted to cluster 1154052 (the two smoke tasks were published by their smokes).  Proxy 453,134 s left.
- 2026-09-25T15:00Z, T&P v6 INSPECTION -> v7 (user decision), LIKELIHOOD EXTENSIONS, Z+X SHAPE FIX:
  * T&P v6 fits_i1 (calibration v4 histograms, adaptive smoothing, EGM CMSShape ranges): SF(full) (45, 1.2)
    muon 0.99397 +- 0.00015, electron 1.00038 +- 0.00038.  Eight inspectors (notes files per step in
    production_v3/tmp/tnp_inspect_v6/fits_i1/, the inspectors now write notes after every page): no nominal
    fit failed outright; recurring: (1) data fail peaks narrower/taller than the MC fail template (the
    Gaussian can only broaden), (2) |eta| 1.4-1.6: the FSR-hump-to-peak ratio of the fail is smaller in data
    than in the MC template (a real data/MC fail-shape difference), (3) sparse single-bin fail templates at
    pT 15-35 carry MC statistical wiggles, wide merges at high pT distort, the inclusive template at low pT
    lacks the continuum, (4) the DSCB alternative puts the pass radiative continuum into the pass background
    and the fail continuum into its Gaussian or the background (most shifts cancel at the SF level; about 25
    visibly failing alternatives), (5) inflated or absurd weighted-MC errors (up to +-1.8) from the numerical
    sandwich; electron ECAL-transition bins (1.479-1.7, pT < 25) where the fail background shape decides
    the result (ee_full data bin 15: SF 0.923 vs ~1.00 from the alternatives).
  * Development tests (production_v3/tmp/tnp_dev/fit7-fit9): a template-to-data statistics ratio (4x) merged
    the dense pass templates and distorted them (rejected); Abramson alpha = 1 oversmoothed the pass tails
    (rejected); alpha 0.5 with a nearest-neighbour kernel floor of 100 effective entries keeps the pass fits
    and removes the fail wiggles.  The pass-like admixture phi of the fail template: MC fits give phi ~ 0
    (-0.06..0.05) with their own templates, data fits 0.17-0.9; the previously bad fail fits (hump too high,
    peak too wide) follow the data; efficiencies move by 0.1-1 % (low pT most).
  * QUESTION (Chinese) -> ANSWER "加入类 pass 成分（推荐）": nominal fail signal = (1 - phi) fail template +
    phi pass template of the bin, phi free in [-0.3, 0.9], convolved with the Gaussian; AGENTS.md updated.
  * v7 (tnp_ul16_v7.json, fit_tnp.cpp, zpeak_model.h, binned_fit.h, run_tnp.py): the admixture; kernel floor;
    min_effective_entries 1000; |eta|-column merges within a factor 2 then 4 in pT; effective-weight errors
    (s^2 = sum w^2 / sum w) for weighted MC instead of the sandwich (binned_fit WeightedErrors::Effective, T&P
    only); the error fallback MINOS -> HESSE -> flagged MINOS interval (error_source); the DSCB alternative with
    a pass low-mass Gaussian and the fail one down to 55 GeV; alternatives whose data or MC error exceeds 5x
    the nominal kept out of the systematic.  Run v7 (new histograms: the rebuilt tnp_histograms binary changed
    the frozen program of run v6) fits_i1 started.
  * Likelihood (h4l_likelihood.py, fit_model.py): POI schemes final_state, category, fv, mode (mu_ggH, mu_VBF,
    mu_VHhad, mu_VHlep from the generator-level production class of every signal yield), stxs0 (|y_H| < 2.5,
    forward parts at the SM), fid_fs and fid_<obs> (sigma_fid / sigma_fid^SM per generator final state or
    differential bin, non-fiducial signal a fixed fraction of the fiducial one per channel, VH non-resonant
    part at the SM); --fix-mh; --final-state (one final state's channels); the on-shell width (DCB x
    Breit-Wigner, h4l_shapes.dcb_bw_pdf, 1D/2D) with fit_width.py (Gamma_H scan, 95 % CL limit, no
    interference: substitution).  build_model.py: --channels categories|inclusive|pt4l|njets|ptj1 and the
    generator-level breakdowns of every signal yield (h4l_fiducial.py: the paper's differential bins, central
    jets pT > 30 GeV and |eta| < 2.5).  fiducial_xsec.py: SM fiducial cross sections from the skim GenTables
    (the delivered MC carries the production preselection: the tables hold 89-92 % of sum(genEventSumw), the
    missing events essentially all outside the fiducial volume; documented).  gof.py (saturated GoF with toys),
    z4l_mass.py (Z -> 4l m_Z per final state; development on v4: 91.23 +- 0.18 (stat) +- 0.07 (syst) GeV).
    Development Asimov fits on v4 (dev_model4, 3D refit): mu 0.996 +0.201/-0.181 (stat +0.182/-0.166, syst
    0.078), m_H 124.999 +- 0.235, Z 10.3; mode mu_ggH 0.997 +0.258/-0.241, mu_VBF 0.99 +1.64/-1.16, mu_VHhad
    0.98 +5.0/-3.6, mu_VHlep 0.99 +3.4/-2.3; fv mu_F 0.997 +0.249/-0.231, mu_V 0.99 +1.29/-0.94.
  * Z+X shape: the Landau fit of the signed OS-method distribution failed (4e/2e2mu chi2 194/78, 145/78; window
    fractions 0.43/0.56 against 0.25/0.31 in the rows themselves): replaced by the paper's procedure, a Landau
    plus exponential fitted to the average of the unit-normalized OS and SS distributions (70-400 GeV), and
    the window fraction taken from the predicted rows (weights of the yield combination); h4l_shapes.zx_pdf in
    the likelihood, the toys and the pyhf export.  dev_zx4 (v4): window yields 4mu 2.82, 4e 2.52, 2e2mu 6.98.
- 2026-09-25T15:03Z, H4L_SELECT v5 COMPLETE (cluster 1154052): scan_select v5: 74 outputs valid, coverage complete, no
  logged ROOT errors; SR entries pseudo-data 868 (v4 869), ZZTo4L 1109402, ggH 33268, VBF 34265, VH 31407.
  Downstream on v5 started (tmux): fiducial_xsec gen_v1, signal_model sm_v1 (categories v2), z4l_mass z4l_v1,
  zx_estimate zx_v1 and zx_closure closure_v1.  Calibration v4 halves closures (uniform, sloped; closure_ul16_v2)
  finished 10 iterations at the statistical plateau (worst residuals 1.2 / 1.7 sigma, scale chi2/ndf ~1.1-1.2):
  finalize and closure_report running.
- 2026-09-25T15:22Z, QUESTIONS (Chinese) on the RESULT.json calibration block: (1) how to combine the flavours
  -> ANSWER "全部分开填写，如果RESULTS等需求和我们等结果有矛盾的，记得问我而不是自己猜测着填写" (fill every
  per-flavour value separately; ask whenever a required format conflicts with our results); (2) what the
  required single-value fields calibration.scale_shift / calibration.smear hold -> ANSWER "填 null" (value and
  unc null, the structure kept; the per-flavour values as extra keys scale_shift_muon, scale_shift_electron,
  smear_muon, smear_electron).  AGENTS.md updated.  Systematics v2 (make_systematics.py, systematics_ul16_v2.json):
  lepton scale mu 3.7e-4, e 5.0e-4 (stat + iteration + the largest closure-profile rms below 100 GeV);
  resolution mu 2.6 %, e 3.4 % (smear-variance statistics and closure deviation propagated to the m4l width,
  plus the residual dilepton width difference after the calibration); flavour fractions from the MC m4l
  responses (2e2mu scale 0.47 mu / 0.53 e, resolution 0.35 / 0.65).  Calibration v4 halves closures at (45
  GeV, 1.2): uniform mu scale -0.014864 +- 0.000109 (injected -0.015), e -0.012123 +- 0.000212 (-0.012), smear
  mu 0.00987 (0.010), e 0.01154 (0.012); sloped mu -0.009814 +- 0.000107 (-0.010), e -0.014038 +- 0.000216
  (-0.014), smear mu 0.00730 (0.008), e 0.00718 (0.010, variance pull -2.2).
- 2026-09-25T15:25Z, USER: "记得及时更新要求，并且定期回顾要求、plan和log，避免搞错和忘记" (update requirements promptly;
  review the requirements, the plan and the log regularly).  AGENTS.md "Records and review" reaffirmed (also while
  long work runs); reviewed AGENTS.md, PLAN.md and this log: the decisions of today (categories v2, Z+X closure
  systematic, T&P pass-like admixture, per-flavour calibration values with the RESULT.json required fields null)
  are recorded in AGENTS.md, the PLAN decision table (updated: calibration block, T&P model and smoothing,
  weighted-MC errors, lepton scale/resolution nuisances) and here.  Pending from the review: eval_selfreport's
  single resSmear field conflicts with the per-flavour rule (ask before writing it); N-1 expected
  significances for the declared cuts; the stage validation plots of the v5 selection; the final reviews.
  Z+X v5 (zx_v1): 4mu OS 4.96 +- 2.10 +- 1.60, SS 9.31 +- 0.65 +- 3.00, combined 6.81 [2.33, 12.38]; 4e 10.91 /
  11.43 -> 11.16 [7.21, 15.21]; 2e2mu 26.92 / 21.32 -> 23.38 [14.36, 36.04] (m4l > 70 GeV); MC closure r = 0.679
  +- 0.099 (OS), 0.678 +- 0.086 (SS), relative systematic 0.322.
- 2026-09-25T15:38Z, QUESTION (Chinese) on eval_selfreport.json: the single measured.resSmear {value, unc} against
  the two per-flavour smears (mu 0.0109 +- 0.0008, e 0.0097 +- 0.0029) -> ANSWER "填 null + 额外键": resSmear value
  and unc null, extra keys resSmear_muon and resSmear_electron; muScaleShift/elScaleShift and
  muEffScale/elEffScale per flavour (the effScale fields = the T&P SF of the full single-lepton selection at 45
  GeV, |eta| 1.2).  AGENTS.md updated.  Development chain on v5 (dev_model5, provisional run-v6 SFs, systematics
  v2, 3D refit): data mu 0.885 +0.195/-0.174 (stat +0.178/-0.161, syst 0.072), m_H 124.99 +0.36/-0.39, Z 8.3;
  Asimov mu 0.996 +0.200/-0.180, m_H 124.999 +- 0.23, Z 10.2.  The Asimov mu of 0.996 (not 1.000) is checked
  next (the Asimov construction must return the injected point).
- 2026-09-25T15:58Z, FIX (3D likelihood normalization): the Asimov of dev_model5 returned mu = 0.996.  Cause: the 3D
  signal DCB width was s x dm (the absolute per-event error, dm = D_mass x m4l), which grows with m4l at fixed
  D_mass = dm/m4l (the third observable of the template P(D_mass)); the conditional m4l density at fixed D_mass then
  integrates to 0.993-0.998 over 105-140 GeV (Riemann check per D_mass bin, production_v3/tmp/diag/asimov_check2.py),
  so the signal pdf was not normalized (the Asimov carried c x nu_s signal events; the data fit was slightly
  biased as well).  Fix: the signal pdf is conditional on D_mass, width = s x D_mass x m_H (CMS L(m4l | m_H,
  D_mass)); signal_model.py fits s with D_mass x 125 GeV (sm_v2), h4l_likelihood.resonant_density uses
  s x (dm/m) x m_H, the toys generate m4l with s x e x m_H and store dm = e x m4l (dm/m4l = the drawn D_mass).
  Check on dev_model6 (sm_v2, v6 SFs): the Asimov grid sums equal the expected yields to 2e-5, the NLL is
  symmetric around mu = 1; Asimov fit mu 1.0000 +0.201/-0.180 (stat +0.183/-0.166, syst 0.076), m_H 125.0012 +-
  0.226 (the 1 MeV offset is the 0.25 GeV grid); data mu 0.887 +0.196/-0.174, m_H 124.965 +0.395/-0.392, Z 8.3.
  The per-event scale s of sm_v2: 4mu 1.083, 4e 0.886, 2e2mu 0.994.
- 2026-09-25T16:40Z, N-1 STAGE PREPARED (eval_selfreport N-1 expected significances): make_nm1_configs.py writes
  seven variants of selection_ul16_v2.json (analysis_v3/reconstruction/config/nm1/), each relaxing one requirement
  to the event-record floor of reco_ul16_v2.json (iso removed; SIP < 8; 30 < m_Z1 < 130; 3 < m_Z2 < 130 in every
  region, smart cut included; lepton pT 3/5 GeV; leading pT 15/7 GeV; m_ll(OS) > 2 GeV); lepton ID, trigger,
  dxy/dz (their floors are the cuts) and m4l > 70 have no variant.  make_select_plan.py: --selection-config,
  --samples (all outputs of each sample), versions with a suffix (v5_nm1_<cut>), the samples filter recorded
  (plans made before that change carry none; their tasks are exactly the nominal MC tasks).  Plans
  production_v3/h4l_select/v5_nm1_<cut> (63 MC tasks each; the v5 frozen program, event records, payloads, refit
  inputs and MELA block; only selection_config differs).  Smokes (sel_GGZZ4E_0000, all seven) passed: SR rows x1.064
  (iso), 1.095 (sip), 1.004 (z1), 1.023 (z2, SR = SRZ4l), 1.021 (pt), 1.001 (leadpt), 1.010 (osmass) of v5; peak
  memory 430 MB.  nm1_significance.py: MC counting in 118-130 GeV (m4l), S = ggH + VBF + VH, B = qqZZ + ggZZ + DY +
  ttbar MC, no SFs, Z = sqrt(2((S+B)ln(1+S/B)-S)), MC-stat recorded, provenance of every variant checked against
  the nominal (program, plan hashes, configuration differences limited to the relaxed keys).  h4l_select_io.read:
  optional lumi for MC-only scans (checked against the data luminosity when the scan has data).  Production review
  (subagent): no blocker; should-fix items done (records, provenance checks); noted: do not run make in
  analysis_v3/reconstruction before staging (headers changed; it would rebuild h4l_select and break the frozen
  binding); submit in waves; the cut named for "leadpt" must avoid the harness substring "pt".  Physics review
  pending; the relaxation-to-floor method is a method choice to put to the user.
- 2026-09-25T16:40Z, T&P EFFICIENCY CLOSURE (prompt section 8, thinned MC): tnp_histograms closure mode (MC jobs with
  "closure": parity of the original file key and keep_id per flavour: a passed identification of a probe is dropped
  with probability 1 - keep_id, per lepton from the new h4l::object_uniform of stream 4), run_tnp.py --closure/
  --closure-point (data role = even MC half with the thinning, MC role = odd half, each normalized to its own
  genEventSumw, templates from the odd DY half; closure bins of tnp_closure_ul16_v1.json; the bin-specific fit
  overrides and exclusions dropped; closure.json compares the recovered SF(full) chain and direct with keep_id);
  an existing run keeps the frozen tnp_histograms of its run.json (the rebuilt program must not break run v7).
  Points k1 (1, 1), k2 (mu 0.97, e 0.95), k3 (0.93, 0.90); k1 histograms: halves 21.98M / 21.89M combinations, the
  weighted id pass counts 36.97M / 36.98M, id efficiency 0.99497 / 0.99498.  Fits running (tmux tnp_closure).
- 2026-09-25T16:40Z, LIKELIHOOD: scheme fid_int (the paper's integrated fiducial fit: r_fid with the 4mu and 4e
  fractions floating, r_fid_2e2mu = (S r_fid - s_4mu r_4mu - s_4e r_4e)/s_2e2mu from the SM fiducial cross sections,
  fit_model.py --fiducial).  Development (dev_model6_incl, 1D m4l without the refit, m_H 125.09): Asimov r_fid
  1.000 +0.226/-0.198 (stat +0.194/-0.176); data r_fid 0.757 +0.201/-0.174 (r_4mu 1.17, r_4e 0.62).
  scan_likelihood.py (1D/2D profile scans, stat-only option).  The fiducial fits follow the paper: 1D m4l without the
  Z1 refit, no categories, no D_bkg^kin, m_H = 125.09 GeV (and profiled); substitution: the reconstruction-level
  candidate of the fiducial fits is the D_bkg^kin-chosen one of the other measurements (the paper uses Z1 closest to
  m_Z and the largest Z2 pT sum there), affecting only events with several candidates.
- 2026-09-25T16:21Z, N-1 WAVE 1 SUBMITTED: h4l_select_v5_nm1_iso cluster 1154122 (62 jobs + the smoke output), _sip
  cluster 1154123 (62 + 1); dry runs ok (Iwd and logs in the /eos/user run directories); proxy 440160 s.  Waves of two
  stages (production review: EOS load of the MELA hashing).  Physics review (subagent): configs, floors (h4l_reco.cpp),
  smart cut, region isolation, plans and normalization correct, the Condor stage can run; BLOCKING for the
  self-report numbers: the reducible part from the raw DY/ttbar SR rows (5 DY entries of weight 1.31 in the nominal
  window; sigma_Z(MC) ~ 0.3) cannot decide the small-effect cuts; suggested: the MC fake-rate method of stage 4b on
  each N-1 output's own ZL and CR trees (numerator/denominator following the variant; fake-rate pT bins below 5/7
  GeV for the pt variant), or the data-driven Z+X per variant; report sigma_MC of Z_with, Z_without, delta Z.
  Should-fix: provenance guards (done, extended to the manifests and normalizations next); the lepton ID can be
  varied (the records keep all loose electrons/muons; the tight ID is applied only in h4l_select; a program option
  would be needed); "without" is relaxed-to-floor (SIP only to 8).  Minor: the z2 variant's effective floor is 4 GeV
  (the OS-pair cut), the Z upper edges do not act in the window, m4l > 70 is determinable (no effect in 118-130),
  cross_clean_dr and event.trigger are inert keys, floors = thresholds while the payload differs (negligible).
  Questions to the user (Chinese): the reducible-background method of the N-1 significances, the meaning of
  "without" (floor relaxation), and a lepton-ID variant.
- 2026-09-25T16:35Z, QUESTIONS (Chinese) on the N-1 significances -> ANSWERS: (1) reducible background "数据驱动 Z+X"
  (the data-driven OS/SS Z+X of each N-1 selection; so every N-1 plan needs the pseudo-data tasks); (2) "放宽到下限并注明"
  (without = relaxed to the event-record floor, relaxed_to stated per cut, null with the reason where the floor is
  the cut); (3) "增加 ID 变体" (muon-ID and electron-ID N-1 variants through a program option; one frozen program for
  all N-1 selections plus a nominal reference checked against v5).  AGENTS.md updated.  Consequence: the MC-only
  plans v5_nm1_<cut> are superseded by full plans v5_nm1d_<cut> (data + MC, the new program); the running MC-only
  clusters 1154122/1154123 (iso, sip) are left to finish as a reproducibility cross-check of the MC part; the other
  five MC-only plans (z1, z2, pt, leadpt, osmass) are never submitted.
- 2026-09-25T16:40Z, USER: "为什么会这么久啊！... TnP拟合 ... 我之前都是一两个小时就全搞定了 ... 很多不都应该可以并行吗？赶紧的";
  "你就写AN就要6个小时？？？"; "对于筛选的优化，不需要那么追求极致 ... 优化的差不多即可了"; "一定要快一点，最后审核改成2个subagent吧".
  AGENTS.md: final reviews by two subagents; speed rule (parallelize every independent piece, Condor chunks, all
  local cores; the AN generated from the result files); no further selection fine-tuning (the N-1 table is
  reporting only).  Cause of the slow T&P: one local process per (flavour, step) fitted its 170-190 bins serially
  (six template-convolution fits per bin with multi-starts and MINOS).  Fix: fit_tnp runs as a Condor task
  (--task/--out-json aliases, the report records frozen_program, task_id and bins); run_tnp.py --condor plan writes
  the fits as chunks of 4 probe bins (each chunk refits "all" first: the per-bin fits start from the inclusive fit
  only, so the chunks reproduce the serial fits exactly), --condor merge validates and merges the chunk outputs
  into fits/fit_<f>_<step>.json and assembles the SFs.  v7 fits_i1c: 360 chunks (stage tnp_fits_v7_fits_i1c).
  Also: the 3D likelihood's constraints now have global observables (post-fit Asimov: centred at the fitted
  nuisances; toys: drawn from N(theta, 1)); scheme mh_fs (m_H per final state in one fit, mu common) and the
  post-fit Asimov option (fit_model.py --dataset asimov_postfit --postfit-from).
- 2026-09-25T16:43Z, SUBMITTED: T&P v7 fits_i1c on Condor (cluster 1154129, 359 chunks of 4 bins + the smoke chunk
  fit_mm_iso_c046, peak 356 MB; microcentury); N-1 selections with data and the new h4l_select (6c2a8b09...:
  leptons.muon.require_id / electron id "none"): v5_nm1d_ref (the nominal configuration, to be checked row by row
  against v5), _iso, _sip, _z1, _z2, _pt, _leadpt, _osmass, _muid, _elid = clusters 1154141-1154150 (73 jobs each +
  the data smoke sel_pseudo_data_0000 of every stage, all passed).  zx_estimate.py/zx_closure.py follow the scan's
  selection configuration (configure: the fake-rate denominator SIP threshold and the lowest pT bin from the lepton
  thresholds; unchanged for v5) and report the 118-130 GeV window fraction (table_window_fraction); build_model takes
  the fake-rate bins from zx.json.  zx_v2 (v5, the nominal re-run with the configured code) running.
  nm1_significance.py rewritten for the user decisions (data-driven Z+X per selection, reference check, ID variants).
  plot_tnp.py draws the galleries in parallel (one worker per flavour, step and role).  make_results.py and
  run_results.py written (every result job; RESULT.json, eval_selfreport.json, results_full.json, summary.md,
  key_numbers.csv, cuts.csv, systematics_coverage.csv, plots/paper_comparison); the paper values in
  analysis_v3/inference/config/paper_jhep11_2017_047.json.
- 2026-09-25T17:25Z, USER (after the usage-limit pause): fewer subagents (the 16 inspection agents stopped at the
  limit, their partial notes discarded); T&P iterations refit and inspect only the bins that are not fine; "我觉得现在的
  TnP拟合其实已经可以了，不用改进了": run v7 fits_i1c is the final tag-and-probe (SF(full) at 45 GeV, |eta| 1.2: muon
  0.99408 +- 0.00019, electron 1.00052 +- 0.00046; inclusive 0.99374 / 1.00004); simple image checks use a cheap
  model.  AGENTS.md updated.  The chunked Condor fits reproduce the serial fits exactly (mm_id, 463 entries, max
  |d eff| = 0); the local serial run fits_i1 was stopped (superseded by fits_i1c).  Final models (sm_v2, zx_v2,
  v7 fits_i1c SFs): model_cat_r0, model_incl_r0, model_pt4l_r0, model_njets_r0, model_ptj1_r0 (eff_correction
  inclusive 0.9959).
- 2026-09-25T17:40Z, RESULTS r1 (running; final models model_*_r0 = sm_v2, zx_v2, T&P v7 fits_i1c): nominal set (b)
  3D refit data: mu 0.888 +0.190/-0.169 (stat +0.178/-0.161, syst 0.057), m_H 124.964 +0.396/-0.390 (stat +0.394/
  -0.387, syst 0.050), Z 8.39; dominant impacts eff_mu +-0.041, eff_e +-0.035 on mu.  MODEL.json (make_pyhf, binned at
  the fitted m_H): mu_hat 0.921 against 0.888, tolerance 0.095: consistent.  Yields (dist_r1): m4l > 70 GeV expected
  848.2, observed 868; 118-130 GeV expected 62.2 (signal 41.5), observed 61.  T&P efficiency closure (k1/k2/k3, full SF
  at 45 GeV, 1.2): muon recovered 1.00003/0.97031/0.93038 +- 0.0003 for 1/0.97/0.93 (slope 0.995 +- 0.005), electron
  0.99989/0.94956/0.89888 +- 0.0004-0.0012 for 1/0.95/0.90 (slope 1.010 +- 0.007).  Toys (r1_toys_mu, r1_toys_mh,
  r1_paired) running; the N-1 selections scanned and their Z+X (zx_nm1, 100 bootstrap replicas) computed as the
  Condor clusters finish (tmux nm1_scan).
- 2026-09-25T17:55Z, QUESTION (Chinese) -> ANSWER "下限 0（同 paper）（推荐）": the per-category, per-mode, (mu_F, mu_V),
  STXS and fiducial POIs lie in [0, 20] as in the paper (the empty-category POIs had run to the -5 limit: without
  events the extended term mu s + b is unbounded below); the inclusive mu stays unbounded.  fit_model.Fitter: every
  extra strength/ratio POI limited to [0, 20]; the fv scan grid from 0.  The affected r1 fits (category, mode, fv,
  stxs0, fs_mu, fid_*) moved aside as <label>.superseded_poi_bounds_20260925 and rerun (run_results tag r1, the
  inclusive fits unchanged: their scheme has no extra POIs).  Observed so far (set b): per final state m_H 4mu
  124.60 +0.42/-0.40, 4e 129.17 +0.79/-1.01, 2e2mu 124.63 +0.76/-0.87 GeV (separate fits); Table 6: 1D refit
  124.56 +0.36/-0.33, 2D refit 124.56 +0.37/-0.34, 3D refit 124.96 +0.40/-0.39, 3D no refit 125.04 +0.41/-0.44;
  mu at 125.09 GeV 0.881 +0.187/-0.167; set (a) mu 0.876 +0.201/-0.175, m_H 124.96 +- 0.39.
- 2026-09-25T18:45Z, RESULTS r1 ASSEMBLED (production_v3/results/v5/r1/: RESULT.json, MODEL.json, eval_selfreport.json,
  results_full.json, summary.md, key_numbers.csv, cuts.csv, systematics_coverage.csv, plots/paper_comparison):
  set (b) mu 0.888 +0.190/-0.169 (stat +0.178/-0.161, syst 0.057), m_H 124.964 +0.396/-0.390 GeV, Z 8.39 (expected
  10.39), GoF p 0.407 +- 0.022 (499 toys, 735 bins); expected m_H +-0.224 (pre-fit), +-0.237 (post-fit); set (a) mu
  0.876 +0.201/-0.175; mu at 125.09 GeV 0.881 +0.187/-0.167.  Per final state m_H 4mu 124.60, 4e 129.17, 2e2mu 124.63
  GeV (compatibility p = 0.003: the 4e events near 129 GeV; Z -> 4e 90.51 +- 0.57 GeV shows no electron-scale
  shift); Z -> 4l combined 91.225 +- 0.182 +- 0.035 GeV (PDG p 0.84).  Width: the DCB (x) BW convolution by FFT
  with a per-final-state cache (identical to np.convolve to 1e-15, 3 ms instead of seconds; non-finite parameters
  rejected; the likelihood rejects any non-positive or undefined event density), best of three starts, grid to 10
  GeV (the first r1_width, grid to 5 GeV and a local minimum, moved aside as .superseded_grid_20260925): observed
  Gamma_H < 6.8 GeV (best 1.2), expected < 2.14 GeV at 95 % CL.  N-1 (nm1_v1, data-driven Z+X per selection, the
  reference identical to v5 row by row): Z with 7.38 +- 0.18; without iso 5.14, sip(8) 7.07, z1(30-130) 7.35, z2(3)
  6.38, pt(3/5) 6.44, leadpt(15/7) 7.38, osmass(2) 7.39, muon ID 7.17, electron ID 4.31.  Toys: mu slope 1.02 (0-3),
  m_H slope 0.993 (121-129); pulls ~0, widths 0.84-1.15 (mu >= 0.5); mu = 0: Z > 3 in 0.5 %, the mu interval not
  regular there (m_H undefined); paired injection dmu 1.009 +- 0.018 (1 -> 2), 1.063 +- 0.021 (0 -> 1: the mu = 0
  base with m_H floating on background, median mu -0.07).  Final reviews (two subagents: physics, deliverables)
  launched.
- 2026-09-25T19:05Z, FINAL REVIEW (deliverables, sonnet): no blocker; RESULT.json contract, numbers traced to their
  sources, eval_selfreport schema and cut-name mapping (no unintended harness matches), N-1 numbers, MODEL.json
  (independent pyhf refit mu 0.9226, within 0.095 of 0.888) all verified; should-fix: (1) the reconstruction step is
  absent from the SF chain without a note; (2) the README status table is stale.  Checked: the Z -> ll yield ratio
  data/MC of the calibration pairs is 0.982 (mm) and 0.969 (ee) against the T&P SF_full^2 0.988 and 1.001.
  QUESTION (Chinese) -> ANSWER "只用 T&P，产额比作交叉检验"; USER: "RECO efficiency of lepton can only be measured from MC",
  "TnP没办法测量，就用MC就行了": the reconstruction efficiency is the MC's, the yield ratio a documented cross-check.
  AGENTS.md updated; the results documentation (RESULT.json conventions, summary, AN) states it.
- 2026-09-25T19:40Z, FINAL REVIEW (physics): 2 blocking, 9 should-fix, 7 minor (production_v3/tmp/final_review/physics.md).
  Blocking: (B1) the 3D Asimov represented every D_mass quintile by its midpoint (the first bin at e1/2), so the 3D
  expected uncertainties were too optimistic (m_H 0.224 against the toy median 0.305); (B2) Table 6's 2D fit must be
  L(m4l, D_mass) (the code's 2D is L(m4l, D_bkg^kin)).  Should-fix: trigger efficiency; the Z+X -1 sigma yield
  (kappa_low used as central/(1 + kappa_low) instead of central (1 - kappa_low)); over-correlated efficiency
  nuisances (every bin's statistical error moved together); invalid width-scan points; fiducial Table 5 and the
  final-state compatibility; paper-style signal strengths from L(m4l, D_bkg^kin) at 125.09 GeV (paper Eq. 10.1);
  validation of the per-event-width 3D model; more toys; the m4l window.  Minor: mu limits, the clipped Z+X category
  fraction (+1.6 % 4mu; violates the signed-template rule), mu = 0 coverage, STXS = mode, njets bin 0, the comparison
  figure, the set (a) VH QCD scale.  QUESTIONS (Chinese) -> ANSWERS: trigger efficiency "用 MC（推荐）"; window "保持
  105–140（推荐）".  AGENTS.md updated.  Fixing every other item now (new models and a full rerun as tag r2).
- 2026-09-25T20:05Z, FIXES of the physics review, rerun as results tag r2 on models model_*_r1: (1) Z+X -1 sigma = central
  (1 - kappa_low) (likelihood and pyhf normsys lo); (2) efficiency nuisances = the coherent T&P fit-model shift in
  quadrature with the statistical part summed over independent bins (h4l_sf.efficiency_variation; 4e Untagged ggH
  +6.7/-6.5 % instead of +8.3/-7.9 %, the fit-model part dominating: its per-lepton mean 3-5 % below 20 GeV, 0.1 %
  above 30 GeV; the paper quotes 2.5-9 % on the yields); (3) the Z+X category fractions renormalized after zeroing a
  negative one (the final-state total preserved); (4) the 3D Asimov samples D_mass at E_SUB = 5 points inside every
  bin's raw quantile range (sum of weights = the expectation to 2e-6; the m_H curvature now gives 0.30 GeV, the toy
  value); (5) the dimension 2Dmass = L(m4l, D_mass) (Table 6's 2D; the old 2D = L(m4l, D_bkg^kin) kept as an extra
  row); (6) paper-style signal strengths (paper_mu, per final state, category, mode, fv, STXS) from L(m4l, D_bkg^kin)
  at m_H = 125.09 GeV (paper Eq. 10.1), the (mu_F, mu_V) contour likewise; the benchmark mu stays the 3D fit with m_H
  floating; (7) the width scan refits invalid points from other starts; (8) mu limits (-20, 50) (never reached);
  (9) systematics_ul16_v3.json: set (a) VH QCD scale 0.017 (sigma-weighted WH/ZH), the fits' default.
- 2026-09-25T21:40Z, RESULTS r2 (production_v3/results/v5/r2; models model_*_r1; 69 jobs, none failed): set (b) mu 0.903
  +0.189/-0.169 (stat +0.180/-0.163, syst 0.051), m_H 124.964 +0.396/-0.388 GeV (syst 0.051), Z 8.46 (expected 10.08);
  expected sigma(m_H) 0.303 pre-fit, 0.319 post-fit (now matching the toys); GoF p 0.422 +- 0.022 (491 toys); MODEL.json
  mu_hat 0.928 against 0.903 (tolerance 0.095): consistent; set (a) mu 0.890 +0.200/-0.175; paper-style mu (L(m4l,
  D_kin) at 125.09) 0.887 +0.186/-0.167.  Table 6 with the paper's 2D = L(m4l, D_mass): 1D refit 124.56, 2D refit
  124.94, 3D refit 124.96 GeV.  Width < 7.05 GeV observed, < 2.15 expected.  Fiducial: sigma_fid 2.50 fb (SM 3.19,
  m_H 125.09), final-state compatibility p 0.465; Table 5 ggH A_fid 0.474, eps 0.667, f_nonfid 0.045 (A_fid relative to
  the LO Pythia8 sample's sigma_eff).  Toys (500 per point, final model): mu 0.5/1/2/3 medians 0.499/0.994/2.004/2.976,
  pull widths 1.02/1.04/0.98/1.00, slope 1.007; m_H pulls -0.03 to -0.10, widths 0.94-1.04; mu = 0: Z > 3 in 0.2 %.
  D_mass validation (dmass_v1): the conditional model describes the signal MC peak per D_mass bin (rms within ~10 %,
  the lowest 4mu bin 1.65 MC / 2.26 model from the tails).  The first make_results attempt of r2 failed on a renamed
  scan label after writing RESULT.json/eval_selfreport.json: both moved aside as *.partial_failed_20260925.
- 2026-09-25T21:00Z, INCIDENT (EOS FUSE): the paper-style reruns failed with "No module named h4l_shapes": the FUSE view of
  analysis_v3/common/python/ showed the directory empty, then three entries (h4l_sf.py, h4l_shapes.py, h4l_style.py)
  with failing stat, while the EOS namespace (eos ls) listed every file intact.  Recovery: the three files copied out
  through XRootD (production_v3/tmp/fuse_rescue/, contents checked) and written back unchanged with xrdcp -f (sha256
  identical; EOS keeps its .sys.v# version entries); the FUSE view is correct again.  No content changed.  The 14
  paper-style jobs (L(m4l, D_kin) without the refit, set (a), 125.09 GeV) rerun.  MC-event closure (closure_mc_v1: the
  weighted MC events themselves as data, truth mu 1, m_H 125): 1D refit mu 1.0007 m_H 125.005; 1D no refit 0.9999 /
  125.002; 2Dmass refit 1.0047 / 125.013; 2Dmass no refit 1.0026 / 125.012; 3D refit (nominal) 1.0181 / 124.982; 3D no
  refit 1.0082 / 124.938: the nominal biases +1.8 % in mu (0.10 sigma) and -18 MeV in m_H (0.06 sigma), documented.
  Final-model toys: m_H 121-129 (400 per point) slope 0.998, pulls -0.06 to -0.01, widths 0.96-1.09; paired injection
  1 -> 2: dmu 0.999 +- 0.010 (0 -> 1: 1.116 +- 0.018, the m_H-floating mu = 0 base).
- 2026-09-25T23:40Z, USER (new benchmark version; verbatim essentials): the scoring outputs and the datasets follow the
  new `ref/Slopbench H4l (Copy).html` (EVAL_CONTRACT v0.4, DATASET.md of task h4l_ntuple); the only deviation: the
  scale, smear and efficiency-SF results per flavour inside `calibration` (no flavour-combined value, the earlier
  null fields and the outside keys removed), each value the data global average; the contract's "### Nuisance
  parameters (optional)" section rewritten accordingly (everything else unchanged) and shipped with the submission;
  the submission (run.sh etc.) runnable from the repository root or the location the HTML requires, local only (no
  Condor), at most ~20 min on a ~10 fb^-1 dataset; truth for the comparison: the injection record
  20260923T183922Z_917f9326efa44ea9.json (format from its generating program if unclear); push the final AN, main code,
  scoring results and the submission to github.com/yiyangzha/slopbench_h4lcontext branch `reference`, lean, nothing
  from ref/.  Answers to my questions: the datasets have exactly the format of our inputs (data and MC ROOT files as
  used here; luminosity and MC entries may differ, so the MC is renormalized); eval fit = binned m4l x D_mass (2D);
  data-weighted average with the Z->ll calibration-sample data leptons per bin; nested calibration format with the
  registered names unchanged.  All recorded in AGENTS.md ("Read access", "New contract", "Evaluation submission",
  "GitHub").
- 2026-09-25T23:45Z, USER (rule): no access outside this repository except what the user authorized; ask first.
  Record of what I accessed before this rule, while looking for the h4l_ntuple dataset schema (nothing of it used in
  any result): slopbench_code_main h4l_ntuplize.py (head), benchmark/branches.py, benchmark/build_mc.py, LATEST.md
  (head), a grep for "h4l_ntuple" in the local slopbench clones, their git branch lists, and directory listings of the
  yiyangz-contact group area top level, h4l_seeds_v3/ (top level) and /eos/home-y/yiyangz/codex/.  Recorded in
  AGENTS.md "Read access".
- 2026-09-25T23:30Z, FIX (fit_model.py, boundary minima): the paper-style category fit (set a, no refit, 125.09 GeV)
  stayed invalid (EDM 0.0019) with four strengths at their lower limit 0 (VH-hadronic, VH-leptonic, ttH, VH-MET).
  Diagnosis: Migrad converged (EDM 4e-5); HESSE at the limit of the bounded parameters' internal transformation raised
  the EDM above the goal and MINOS then refused every POI (the reported intervals were symmetric HESSE errors).  The
  earlier tolerance workaround (tol = 5 for EDM < 0.005) never triggered and is removed.  Root-cause treatment: after
  Migrad, every free parameter within 1e-4 of its range from a limit, where the NLL rises into the allowed region, is
  fixed at the limit and the rest re-minimized (the constrained minimum); its one-sided interval (0 on the limit side)
  is the profile-likelihood edge NLL + 0.5 found by bracketing and bisection (stat: nuisances fixed); MINOS runs on the
  free POIs.  Rerun with the new code: every fit with a POI at a limit (category data and Asimov, mode, stxs0,
  fid_njets, fid_ptj1 data; old outputs moved aside as *.superseded_boundary_20260925 or .superseded_edm2_20260925):
  all valid; category data: VBF-2jet 1.147 +0.752/-0.532, VBF-1jet 0.164 +0.275/-0.146, Untagged 1.080 +0.254/-0.222,
  VH-hadronic 0 +0.544, VH-leptonic 0 +1.784, ttH 0 +4.186, VH-MET 0 +6.141; the others reproduce their earlier MINOS
  intervals within 0.003 (e.g. mode VHhad 0 +1.816).
- 2026-09-26T00:30Z, USER (evaluation submission): the eval fit window stays 105-140 GeV with m_H in [110, 140] GeV (as
  the main analysis; asked because the new contract's physics envelope is 100-200 GeV); tag and probe may be simplified
  ("TnP的要求不要太高，差不多就行了，卷积等可以近似"); large ROOT files must be processed with ROOT, never Python
  ("大root文件处理绝对不能用python，必须用root"); work faster; ~20 min per ~10 fb^-1 eval run is fine; do not change
  validated methods, only equivalent faster implementations ("之前没问题的就基本不用改了，除了提升速度用等效的程序之外，都
  不要改了 ... 不要节外生枝"; the T&P simplification approved).  Consequences (recorded in AGENTS.md "Evaluation
  submission"): the uproot per-file reader replaced by the ROOT/C++ program my_analysis/src/h4l_reader.cpp (compiled at
  run time; 16 parallel processes over file chunks; a data shard ~9 s, an MC file ~0.4 s on EOS FUSE); the eval lepton
  calibration re-implemented as a port of the main analysis's run_calibration.py algorithm with the calibration_ul16_v4
  configuration (factorized model, families A/B, response pass with the decorrelating smear, joint least squares with
  outlier rejection, frozen selection after iteration 2, average of the last 4 applied payloads) with histogram
  convolutions instead of event-level kernels (an earlier simplified variant, whose smear did not converge because of
  pT-bin sculpting and a DY-only template, is kept as production_v3/tmp/calibration_eval_v1_superseded.py); the
  templates are DY + ttbar weighted to their cross sections, as the main analysis's MC role; T&P simplified (coarser
  bins, calibrated MC template without the Gaussian convolution, CMSShape nominal / exponential alternative, MC fits
  with negligible-background shapes fixed); toys and T&P fits run in parallel processes.  Timing on the 1 fb^-1
  development dataset: read 37 s, calibration 41-65 s, T&P 5 s, inference 23 s.
- 2026-09-26T00:45Z, FREEZE (before unblinding): production_v3/results/v5/final_v2/freeze_v2.json, sha256
  d42d3bb3069e057946a09073646d3c207a33fe34678491fc9ce9b0f2e2fe2584 (122 analysis_v3 source files, 15 binaries, the 9
  result files of final_v2: RESULT.json per EVAL_CONTRACT v0.4 with the per-flavour calibration, MODEL.json, ...).  The
  evaluation submission (my_analysis/) keeps the main analysis's methods; after this point it receives only crash fixes
  and speed changes, never a tuning to the truth.  Note: an accidental `pixi run --frozen` of the submission manifest
  created my_analysis/.pixi (environment directory, not deleted per the rules; excluded from git and from the
  submission copy).
- 2026-09-26T01:10Z, UNBLINDING (user-authorized record and generator configuration): truth m_H = 125 GeV, mu = 1 (every
  mode), 20 fb^-1, response profile medium (lepton scale k = 1 - 0.018 - 0.004 S(pT/100) - 0.008 S(|eta|/2.5), relative
  pT smear 0.01 plus 0.01 eta/phi smears), efficiency profile effA (thins muon loose/medium/tight/soft IDs, electron MVA
  WP90/WP80 and cut-based flags, single-lepton triggers; none of this selection's flags).  Truth summary in
  production_v3/results/v5/truth_v1/truth.json (make_truth.py): data-weighted truth scale_shift muon -0.01994, electron
  -0.01991; smear 0.01; sel_eff 1.  Comparison (main analysis, final_v2): mu 0.903 +- 0.179 (pull -0.54), m_H 124.964 +-
  0.392 (-0.09); scale_shift muon -0.01993 +- 0.00037 (+0.03), electron -0.01987 +- 0.00049 (+0.08); smear muon 0.0116 +-
  0.0005, electron 0.0113 +- 0.0021 (above the pT-only truth 0.01 as expected from the angular smears); sel_eff muon
  0.9941 +- 0.0018 (-3.3 sigma from 1), electron 1.0004 +- 0.0031.  AN with the truth section:
  deliverables/AN_h4l_ul16_pfnano_v2/main.pdf (make_an.py adapted to the final_v2 format).
- 2026-09-26T01:05Z, EVAL TEST (production_v3/tmp/eval_runs/test_10fb_4; ~10 fb^-1: 41 shards + all MC, DATASET.md
  layout by symlinks, analysis environment + lxplus ROOT): 1079 s (read 541, calibration 357, T&P 55, selection 94,
  inference 19); mu 0.878 +0.257/-0.224, m_H 125.558 +0.489/-0.480 GeV, Z 5.64 (expected 6.57), GoF p 0.040, coverage
  0.722; calibration muon -0.01996 +- 0.00037 / 0.0121 / sel_eff 0.9944, electron -0.01995 +- 0.00049 / 0.0131 / 1.0023.
  pyhf builds MODEL.json (after dropping the top-level _meta, which the pyhf schema rejects: now in MODEL_meta.json) and
  refits mu 0.8776 (reported 0.8778).
- 2026-09-26T01:15Z, GITHUB (user-authorized): branch `reference` of github.com/yiyangzha/slopbench_h4lcontext (commit
  00ce415, from a separate clone in production_v3/tmp/github_stage; this repository's git untouched): my_analysis/,
  analysis_v3/ (sources, no binaries), results/ (main_final, eval_test_10fb, truth), the AN (main.pdf, main.tex),
  AGENTS.md, PLAN.md, README.md, experiment_log.md, REFERENCE.md.  Verified: no file from ref/ (no ref/ path, no HTML,
  no reference PDFs).  my_analysis/pixi.lock force-added (the base .gitignore ignores pixi.lock).
