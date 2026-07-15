# Phase 1 physics re-review — calibrated CMS Open Data H→ZZ*→4ℓ measurement

**Verdict: ITERATE — not approved for Phase 1 advancement.**

**Commit:** Keep every binding decision, especially the real reconstructed-jet
VBF category and DY+jets nominal fake model. Do not substitute a flat-tree
proxy or advance to a data fit. Obtain the matching raw-data provenance first,
complete the MC-only preflights identified below, then submit a new
independent review.

## Scope and evidence

This is a physics-only re-review of the stated H→4ℓ goal, the original
physics review, the arbiter verdict, STRATEGY_v2, RAW_OBJECT_FEASIBILITY_v1,
COMMITMENTS, and the three named preflight JSON records. I did not use
methodology, conventions, source code, or unlisted artifacts. Phase 1
contains no figures, so there is no data/MC visual-normalisation check to
perform.

The revised package is materially more honest than the original one. In
particular, it explicitly says that its raw-object pilot proves an MC route
only, records A1-data as open in the ledger, and says that missing raw-data
provenance blocks rather than downscopes the VBF/FSR result. That distinction
is correct: raw-MC feasibility is not raw-data feasibility.

## Audit of every earlier A/B finding

| Earlier finding | Concrete v2 evidence | Re-review disposition |
|---|---|---|
| Physics A1 / arbiter A1: raw objects for VBF and FSR | The feasibility record finds no raw-data NanoAOD or manifest for either supplied fake-data product; their flat trees have no Jet, Photon, or FsrPhoton branches. STRATEGY_v2 says a one-to-one raw-data match is required before a VBF data fit and forbids a proxy. | **Still A.** The data prerequisite is externally blocked. In addition, the MC proof is not yet complete; see current A2. |
| Physics A2 / arbiter A2: predictive DY+jets fake model | v2 freezes L/P/F, 4P, 3P1F, 2P2F, and same-sign definitions; defines the transfer T by flavour/category/score/mass cell; retains DY as nominal; specifies prompt-ZZ and double-count subtraction, correlated rate and mass/score shape modes, TTBar separation, and closure/coverage tests. | **Strategy design resolved conditionally.** It must still be exercised with data once raw provenance exists. A small deterministic validation gap remains B1. |
| Arbiter A3: commitment and systematic ledger | COMMITMENTS is present, append-only, and tracks every D item, systematic source, validation, comparison, and flagship figure. It leaves unavailable items open rather than marking them downscoped. | **Resolved at Phase 1 strategy level.** |
| Physics B1 / arbiter B6: calibration convention, control model, and closure | D10 gives separate data and MC transformations to the common external Z reference; v2 defines the flavour/sample fit model, joint covariance including the common MZ term, rebuilds all candidate observables and templates, and predeclares injected e-only and μ-only closure. | **Resolved at Phase 1 strategy level.** No calibrated result is claimed before the specified fits and closures. |
| Physics B2 / arbiter A4: numerical 10 fb⁻¹ scale sanity | The preflight records 4,508,686 ee and 6,284,267 μμ broad-control candidates, statistical relative floors of 2.57e-5 and 2.07e-5, and channel mass floors of 3.22, 2.59, and 2.07 MeV. v2 predeclares joint-covariance propagation, non-overlapping ZFitModel alternatives, injected offsets, pull/coverage gates, and a profiled with/without-scale impact. | **Resolved at Phase 1 strategy level.** The alert thresholds are explicitly not assigned nuisance sizes, which avoids confusing a planning benchmark with a calibration result. |
| Physics B3 / arbiter B8: sparse simultaneous-fit resolving power | v2 fixes the base flavour×VBF categories, 1 GeV templates with 0.5/2 GeV tests, gamma constraints, morph-node closure, occupancy/Hessian/Fisher gates, and 500 event-disjoint toys at four (mH, μ) points. | **Resolved as a preflight contract, not as an achieved sensitivity result.** The current-flat inventory is only 9.24 expected signal events against about 101.77 listed backgrounds and is explicitly not a category table; the eventual result must honestly report limited μ resolving power if the ±20% injections cannot be separated. |
| Physics B4 / arbiter B9: comparison coverage and overlap | The comparison matrix gives HIG-16-041 and PDG targets, matched-mass-convention conditions for μ, an explicit non-comparable world-μ disposition, and a difference-not-pull rule until overlap covariance is known. | **Resolved at Phase 1 strategy level.** |
| Arbiter A5: HIG-16-041 method parity | The parity matrix names the published likelihood, constrained-mass, resolution, categorisation, and reducible-background differences. It binds an independent constrained-mass/reference-like pseudoexperiment comparison with bias and coverage criteria rather than asserting equivalence. | **Resolved as a required future cross-check.** |
| Arbiter B7: angular-NN decision and score transfer | v2 fixes the five rest-frame angles, excludes mass/identity/weights, gives deterministic train/validation/test splits and seeds, includes every fit background, and defines input, score, sculpting, occupancy, and expected-uncertainty gates. | **Resolved at Phase 1 strategy level.** Final score validation remains contingent on the raw-data route. |

## Remaining findings

### A1 — The requested data VBF/FSR programme is externally blocked

**Classification: A — must resolve; external dependency.**

RAW_OBJECT_FEASIBILITY_v1 finds that the supplied fake H→4ℓ and dilepton
products retain event identifiers but lack jets and photons, and that no
matching raw-data NanoAOD directory or recoverable manifest exists in the
supplied production/owner area. STRATEGY_v2 correctly observes that identifiers
cannot recreate these objects.

This blocks more than a jet plot. It blocks the requested reconstructed-jet
VBF category, electron and muon FSR recovery in data, data VBF migrations, the
VBF comparison to HIG-16-041, and the raw-object loose/not-tight fake-control
implementation. The package does not conceal this limitation, and it must not
be treated as a solvable flat-tree analysis task.

**Required external resolution:** the data producer must provide a raw-data
NanoAOD location/manifest that maps one-to-one to fake_data_10fb.root on
(run, luminosityBlock, event), preserves the required lepton, jet, b-tag,
FSR-photon, and generic-photon inputs, and permits the documented data pilot.
The re-review must inspect that persisted match before accepting any data VBF,
FSR, or fake-control result.

### A2 — MC VBF/FSR feasibility is promising but does not yet meet the full prior preflight

**Classification: A — must resolve; internally solvable with the supplied MC.**

The MC evidence is useful but narrower than the claim needed for an executable
fit:

- All 11 processes have a representative raw schema with the required branch
  groups, but the only persisted exact candidate match is 1,198/1,198 VBF
  candidates.
- That exact pilot maps a raw NanoAOD file to an h4l_mc_modified VBF flat
  product, while the supplied analysis MC is named h4l_mc_nominal. No
  corresponding exact-match evidence is shown for the nominal templates or
  the other ten fit processes.
- The 671 in-window two-cleaned-jet candidates are explicitly an object probe,
  not the frozen D_2jet VBF tag, its weighted 10 fb⁻¹ yield, purity,
  composition, migration, or profiled information.
- Muon-indexed FSR is preserved for 37 pilot candidates, but the required
  electron-photon association has not yet been validated even on MC.

The original A1 required stable identities and a quantitative VBF preflight
for all relevant MC processes. The current strategy and ledger correctly list
these as future gates, but their absence means the MC side is not fully
resolved either. This is not externally blocked: run the persisted
identity/object pilot for every actual nominal template and required
variation, validate electron FSR and angle/pairing preservation, implement the
frozen tag, and publish the weighted yield/efficiency/purity/composition,
JEC/JER/b-tag migration, and profiled-information table.

### B1 — The fake-transfer validation needs one fully deterministic, signal-like validation rule

**Classification: B — fix before a future PASS; internally solvable now.**

The new fake contract is a major improvement, but it permits coarser transfer
binning when a denominator has “negligible effective statistics” without a
numerical effective-count threshold or frozen merge priority. At this
luminosity the current-flat DY inventory is sparse in the fit window, so that
choice can materially change a mass/score template.

Also, the held-out same-sign region is a valuable orthogonal check but has a
different charge and potentially different conversion/fake composition from
the OS 4P target. It cannot by itself establish the predictive transfer into
the signal-like OS selection. Freeze a numerical N_eff/merge rule and a
blinded OS 4P sideband or other demonstrably composition-matched held-out
validation, separately from the regions used to derive the transfer. Execute
that check once A1 is resolved; do not call MC split closure an independent
data validation.

### C1 — Preserve the limited interpretation of the Z-control preflight

**Classification: C — apply during the next iteration.**

The broad-window means in the JSON are 85.70 GeV (ee) and 86.19 GeV (μμ), not
fitted Z-peak locations. v2 appropriately labels them only as broad-control
statistics and not as calibration corrections. Keep that discipline: before
calling the control selection clean or quoting a calibration, show the
flavour-separated peak fits, their goodness of fit and purity/contamination
assessment, and the before/after distributions required by the goal.

## Requested-scope assessment

| Requested scope | Physics assessment |
|---|---|
| Lepton calibration and correlated scale fit | The data/MC-to-external-Z convention, covariance, alternative model envelope, injections, and profile-impact definition are sound Phase 1 commitments. They are not yet calibration results, as they should not be. |
| DY+jets fake background | DY remains the required nominal model rather than being silently replaced by a data method. Its control-to-template and correlated nuisance plan is credible, subject to B1 and the raw-data blocker. |
| VBF category and FSR | The no-proxy rule is correct. Neither is feasible for the supplied fake data until external raw provenance arrives; MC has not yet supplied the complete tagged-yield/FSR preflight. |
| Four-lepton-rest-frame angular NN | The five-angle, no-mass design and deterministic quality gates are physically appropriate. Its data score and fake-control validation cannot be completed before A1. |
| Simultaneous mH, μ fit | The workspace/toy/morphing/template-statistics contract is appropriate and avoids claiming sensitivity from a formal likelihood alone. Actual resolving power is unproven and must be reported honestly. |
| Reference comparisons | The fixed 105–140 GeV window, HIG-16-041 parity cross-check, overlap rule, and PDG/non-comparable distinctions are now appropriate. |

## External versus solvable work

| Item | Status | Required action |
|---|---|---|
| Matching raw fake-data NanoAOD/manifest | **External blocker** | Obtain provenance from the data producer and run the documented one-to-one data pilot. |
| Per-process nominal-MC identity/preservation, tagged VBF yield, electron FSR validation | **Internally solvable** | Complete with the supplied raw MC before the next physics PASS request. |
| Deterministic fake-template merging and OS held-out validation definition | **Internally solvable** | Amend the contract now; execute on data after provenance arrives. |
| Z peak fits, NN validation, fit toys, and actual systematic impacts | **Planned downstream work, not evidence yet** | Perform only after the above prerequisites; do not relabel contracts as results. |

## Final decision

**ITERATE.** The revised strategy is honest about its central limitation and
repairs most of the original Phase 1 design deficiencies. It still cannot
advance because A1 prevents the requested data VBF/FSR/fake-control programme,
and A2 leaves readily solvable MC feasibility evidence incomplete. Resolve
A2 and B1 while requesting the external raw-data provenance; after the
one-to-one data pilot succeeds, re-review the complete VBF/FSR and fake-control
evidence before any data fit or claim of a complete H→4ℓ measurement.

