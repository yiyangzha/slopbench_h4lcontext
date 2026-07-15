# Phase 1 arbiter adjudication — calibrated CMS Open Data $H\to ZZ^*\to4\ell$ measurement

**Verdict: ITERATE — Phase 1 must not advance.**

The strategy has the correct intended endpoint and several binding decisions are
already sound: a detector-level simultaneous $(m_H,\mu)$ fit, the required
$105<m_{4\ell}<140$ GeV window, an independent flavour-separated $Z$-peak
calibration, an angular-NN approach with a cut-based cross-check, and DY+jets
as the nominal fake template. Those intentions are not yet enough to make the
Phase 1 programme executable or auditable. The Category A items below must be
resolved and all Category B items fixed before a future arbiter may PASS.

## Inputs, evidence limits, and direct checks

I read the prompt, committed strategy, all three Phase 1 reviews, the Phase 1
and root experiment/retrieval logs, methodology/03-phases.md,
methodology/06-review.md, both candidate conventions, and the arbiter
contract. Phase 1 has no figures, so the absent plot-validator review is
correct under methodology/06-review.md §6.2.

The unavailable corpus connector is explicitly recorded in
phase1_strategy/retrieval_log.md:3-7. It is not a reason to defer this
adjudication: the strategy contains primary-source fallback links. I did not
treat uncross-checked reference numbers as evidence of readiness. A direct
check of ntuplizer/h4l_ntuplizer.cpp:529-602 confirms that the current reader
has run/lumi/event identifiers and lepton collections, but declares no jet,
b-tag, photon/FSR, or MET object collection. A direct file inventory confirms
that phase1_strategy/COMMITMENTS.md is absent.

The method-parity finding is independently supported by the official
[CMS HIG-16-041 result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-16-041/):
the published mass study reports 1D, 2D, and 3D likelihood variants, including
the $\mathcal{L}(m'_{4\ell},\mathcal{D}'_{\rm mass},\mathcal{D}^{\rm kin}_{\rm bkg})$
configuration and an $m(Z_1)$-constrained mass. The strategy currently gives
only a generic binned likelihood, so it has not yet demonstrated parity or
defined the required reference-like cross-check.

## Structured adjudication

Findings with the same root cause are consolidated below; consolidation is not
a dismissal. No reviewer finding is dismissed.

| # | Finding | Source | Their Cat | Final Cat | Rationale |
|---|---|---|---|---|---|
| 1 | No demonstrated raw-object re-ntupling route for the binding VBF category and reference-like FSR reconstruction | Physics A1; Critical A1; Constructive A1 | A, A, A | **A** | The artifact itself acknowledges that the flat tree lacks needed objects (STRATEGY.md:31-38), while making VBF binding (:84-88) and FSR prerequisite reconstruction (:157-163). Direct source inspection corroborates the object deficit. This cannot be replaced by a $p_\mathrm{T}^{4\ell}$ proxy. |
| 2 | DY+jets is labelled as nominal fake template but lacks a predictive control-to-signal transfer, independent closure, shape coverage, and reducible-composition inventory | Physics A2; Critical A2; Constructive A2 | A, A, A | **A** | Naming loose/not-tight and same-sign controls (STRATEGY.md:249-256) does not define regions, transfer, pass criterion, or correlated rate/shape nuisance. The missing contract can bias both the fit and the NN trained against DY. |
| 3 | Mandatory commitment ledger and fully resolved systematic dispositions are missing | Critical A3; Constructive A3 | A, A | **A** | methodology/03-phases.md:771-785 requires COMMITMENTS.md at Phase 1 completion; the executor plan also promised it. Luminosity and generator/shower rows remain open inventories rather than explicit Will implement or justified Not applicable dispositions (STRATEGY.md:298-302). |
| 4 | The requested numerical 10 fb$^{-1}$ lepton-scale sanity test is not specified | Physics B2; Critical B2; Constructive A4 | B, B, A | **A** | I side with the higher constructive severity because the user explicitly requires a defensible numerical scale contribution for $m_H$ and $\mu$. The only quantitative check is the unrelated 0.38 GeV mass-statistical extrapolation (STRATEGY.md:52-57); it cannot bound or validate the scale contribution. |
| 5 | Statistical-method parity with HIG-16-041 is asserted rather than demonstrated | Critical A4 | A | **A** | This single-reviewer finding is valid under methodology/03-phases.md:207-223. The official CMS reference uses richer 1D/2D/3D likelihood variants and a constrained mass; the strategy has neither a material-feature comparison nor a reference-like cross-check. A generic claim to match likelihood treatment is not evidence. |
| 6 | $Z$-control selection, fit model, data/MC template convention, and injected-scale closure are insufficiently specified | Physics B1; Critical B1; Constructive B1 | B, B, B | **B** | The strategy correctly chooses an independent $Z$ anchor, but never fixes the data-versus-MC transformation convention or the quantitative control selection/fit/closure definition (STRATEGY.md:201-240). Without this, nominal and nuisance transformations could cancel or bias the mass fit. |
| 7 | The angular-NN has qualitative gates but no quantitative score-transfer/decision contract or complete fit-background treatment | Physics C2; Critical B3; Constructive B2 | C, B, B | **B** | The higher B assessment is warranted: the NN is binding primary inference information under [D4], yet it lacks acceptance statistics, score-level fake-mixture validation, class-mixture/weight disclosure, a treatment for every fit background, and a deterministic B-versus-A decision rule. |
| 8 | The simultaneous fit has no sparse-category, template-statistics, morphing, identifiability, or independent pseudoexperiment coverage preflight | Physics B3; Critical B4; Constructive B3 | B, B, B | **B** | Six flavour/category cells plus optional score bins at 10 fb$^{-1}$ cannot be accepted merely because a likelihood is written down. Expected per-bin yields, deterministic occupancy/binning rules, mass-morph closure, MC-statistics treatment, and unbiased joint $(m_H,\mu)$ recovery are all missing. |
| 9 | The comparison programme lacks a result-by-result comparability metric, overlap/covariance treatment, and audit-ready reference provenance | Physics B4; Critical B5, C2; Constructive B4 | B; B, C; B | **B** | The references and headline values are useful, but the analysis must say for each reported quantity/distribution whether the data overlap, which statistic is valid, and why a $\mu$ world comparison is or is not meaningful. Citation URLs alone do not supply the required table/HEPData-section provenance. |
| 10 | Common external-$M_Z$ reference uncertainty is not explicitly correlated across $ee$ and $\mu\mu$ calibrations | Physics C1; Critical C1 | C, C | **C** | The effect may be negligible, but it is shared by the two scale equations. Include it as a correlated component or document a numerical negligible-impact demonstration. |
| 11 | Luminosity and generator-coverage feasibility need an early decision deadline | Critical C1; Constructive C1 | C, C | **C** | [A3] and [L2] correctly reject invented uncertainties, but an early Phase 2 deadline will prevent an indefinitely conditional $\mu$ result or an untested model-uncertainty promise. |
| 12 | Re-ntupling should demonstrate that the currently feasible angular calculation is preserved after FSR/calibrated reconstruction | Constructive C2 | C | **C** | Matched-event before/after angle and pairing checks are a low-cost safeguard against accidentally breaking the valid candidate-level angular programme while solving VBF/FSR. |

## Independent binding-decision audit

| Binding item | Direct status | Adjudication |
|---|---|---|
| [D1] detector-level simultaneous template/profile likelihood | Present at STRATEGY.md:64-67 | Retain; method parity is A5. |
| [D2] $105<m_{4\ell}<140$ GeV for every fit | Present at :69-72 | Satisfied in strategy. |
| [D3] independent flavour-separated $Z$ calibration | Present in principle at :74-77,199-240 | Retain; operational definition is B6 and numerical sanity is A4. |
| [D4] angular NN primary with cut-based cross-check | Present at :79-82,165-197 | Retain; deterministic gates are B7. |
| [D5] reconstructed-jet VBF category | Binding at :84-88 | A1 blocks its feasibility; no proxy or silent omission is allowed. |
| [D6] DY+jets nominal fake model | Present at :90-92 | A2 blocks predictive use until transfer/closure is defined. |
| [D7] MC pseudo-data then deterministic 10% data | Present at :94-98 | Satisfied in strategy; later diagnostics remain binding. |
| [D8] supplied cross sections and Metadata.nEvents normalization | Present at :100-121 | Satisfied in strategy. |
| [D9] no particle-level/unfolded claim | Present at :105-107 | Satisfied in strategy. |

## Independent regression and motivated-reasoning check

This is Phase 1, so there is no earlier phase to regress to. The findings
require a Phase 1 **iteration**, not a regression ticket. I independently
checked the arbiter regression triggers:

| Check | Evidence in current Phase 1 record | Result |
|---|---|---|
| Validation failure without three documented remediation attempts | No closure, stress, or data/MC result has yet been claimed; only future tests are named | No current regression trigger; missing test contracts remain A2/B6-B8. |
| Single systematic >80% of total | No uncertainty breakdown or fitted result exists | Not evaluable yet; A4 requires a pre-fit scale-reach calculation and later impact audit. |
| GoF toy inconsistency | No toys or observed $\chi^2$ exist | Not triggered; B8 must predeclare toy-GoF and coverage evidence. |
| More than 50% bin exclusion | No bin-exclusion scheme exists | Not triggered; B8 must make occupancy/binning decisions deterministic. |
| Tautological comparison presented as validation | The $Z$-peak control is independent in intent and the nominal-model statement is labelled a limiting check, not a validation result | No current trigger, conditional on B6 preserving independent data/MC calibration and B9 separating closure from independent comparison. |
| Visually identical supposedly independent distributions or a trivial fit | Phase 1 contains no plots or fitted results | Not applicable at this gate. |

The narrative’s most consequential self-serving shortcut is to treat [L1]'s
luminosity-scaled mass-statistical benchmark as the requested lepton-scale
sanity check. It is not. A future promise to report the scale impact also does
not repair the missing predeclared numerical chain. Conversely, the strategy
does not currently claim an observed result, so no Phase 4 validation-target
pull/deviation rule is falsely passed at this stage.

## Dismissal audit

No finding is dismissed as out of scope or deferred merely because it needs
upstream work. In particular, A1 and A2 are not excused as Phase 2/3 work:
the small feasibility pilots, predeclared gates, and strategy commitments are
exactly the work required to make Phase 1 approvable. The dismissal cost,
physics-conclusion, and future-phase commitments of methodology/06-review.md
§6.5.1 therefore do not apply.

## Exact fix list for the iteration

### Category A — must resolve before Phase 1 may advance

1. **Raw-object VBF/FSR feasibility (A1).** Locate accessible, matching raw
   inputs for data and every fit MC process; make and record a small
   re-ntupling pilot preserving run/lumi/event, selected jets, JEC/JER,
   b-tag, FSR candidates, and required systematic inputs. Reconcile the
   contradictory timing in [A1] versus the handoff, predeclare object cleaning,
   FSR association/recovery, and VBF tag variables/thresholds, and state the
   viability evidence: expected 10 fb$^{-1}$ yields, purity/composition,
   migration uncertainty, and non-prior-dominated fit information. If the raw
   path is unavailable, formally revise the strategy; do not substitute a
   proxy.
2. **Predictive DY-fake contract (A2).** Keep DY+jets MC nominal, but define
   tight/loose and same-sign selections, orthogonality, transfer to every
   flavour/category/score template, independent closure region, expected
   composition/yields, pass statistic and coverage, and correlated rate/shape
   nuisances in both $m_{4\ell}$ and NN score. Inventory top, conversion/
   $Z\gamma^{(*)}$, and electroweak/nonprompt components with a template/control
   or quantified propagated negligible impact.
3. **Commitment and systematic ledger (A3).** Create
   phase1_strategy/COMMITMENTS.md using [x], [ ], and [D] statuses for every
   [A]/[L]/[D], systematic, validation, selection approach, reference
   comparison, and flagship figure. Amend the strategy/ledger so every
   applicable systematic has an explicit action or a documented justified
   downscope, including a decision tree for luminosity and a concrete inventory
   outcome for generator, shower, PDF, and scale coverage.
4. **Numerical scale-reach contract (A4).** Define a pre-fit, separate-flavour
   chain from expected/measured $Z\to ee$ and $Z\to\mu\mu$ fit covariance to
   $c_e,c_\mu$ variations, candidate/template rebuilding, and profiled
   $\Delta m_H,\Delta\mu$. Predeclare a plausible numerical range at 10
   fb$^{-1}$, independent pseudoexperiment/injected-offset coverage tests, and
   a non-overlapping covariance or single-envelope rule for fit model, binning,
   domain, and residual alternatives.
5. **Published-method parity (A5).** Add a source-specific HIG-16-041 parity
   table covering likelihood observables/dimensionality, constrained/refitted
   mass treatment, discriminants/categories, resolution treatment, signal
   shapes, and fake estimate. Either implement material features or commit a
   reference-like fit as an independent predeclared pseudoexperiment cross-check
   with bias and coverage for both $m_H$ and $\mu$.

### Category B — must fix before a re-review may PASS

1. **Calibration convention (B6).** Fix the $Z$ control trigger, object/vertex
   working points, fit domain, signal-plus-background/resolution model, and
   GoF/closure criterion. State the exact nominal and varied transformation for
   data and each MC template (including the role of its fitted MC peak), and
   validate injected electron/muon offsets in all three final states before
   allowing supported $\eta$/$p_T$ subdivisions.
2. **NN decision and score-transfer gates (B7).** Predeclare input and score
   agreement metrics, score–mass/sculpting limit, held-out criterion,
   architecture/seed-spread treatment, class weights/mixture, score-bin
   population rule, and exact formal-revision/fallback trigger. Validate the
   score in sideband and fake-control mixtures, including all fit backgrounds.
3. **Simultaneous-fit preflight (B8).** Publish expected signal/background
   yields per flavour/category/score bin, binning relative to mass resolution,
   minimum template-statistics treatment, $m_H$ grid/morphing closure, and
   independent pseudoexperiment evidence for unbiased joint recovery, finite
   pulls/coverage, toy GoF, no boundary saturation, and sensitivity to at least
   a predeclared ±20% $\mu$ injection.
4. **Comparison matrix and provenance (B9).** For $m_H$, $\mu$, and each key
   distribution, record the target/table, observable/selection equivalence,
   data-overlap and covariance assumption, valid comparison statistic, and
   reason if comparison is not meaningful. Keep HIG-16-041 and PDG-$m_H$
   distinct from calibration inputs, avoid an independence-implying pull if
   overlap covariance is unknown, and add source page/section or HEPData-table
   provenance for every numerical target.

### Category C — apply during the iteration

1. Include the common external-$M_Z$ uncertainty as a fully correlated
   electron/muon component or demonstrate its negligible impact numerically.
2. Add a Phase 2 feasibility deadline for a dataset-appropriate luminosity
   treatment and generator-variation availability, recording the consequences
   if $\mu$ remains conditional.
3. On matched events, compare pairing and angular variables before and after
   raw-object re-ntupling/FSR/calibration so that the VBF/FSR remedy does not
   silently alter the angular-NN basis.

## Required evidence for re-review

The next panel must receive a revised strategy, the new commitment ledger, a
raw-object feasibility/pilot record, a fake-region transfer-and-closure
specification, method-parity table, quantitative calibration/scale plan, NN
and fit pseudoexperiment acceptance criteria, and the comparison matrix with
source provenance. Per methodology/06-review.md:506-514, the next review must
be a written re-review that explicitly verifies every A/B item above.

## Final decision

**ITERATE.** Category A items are A1 raw-object VBF/FSR feasibility, A2
predictive DY-fake modelling, A3 the commitment/systematics ledger, A4 the
numerical lepton-scale sanity contract, and A5 published-method parity. No
regression or escalation is warranted on this first Phase 1 review cycle, but
the phase may not advance until a fresh independent re-review documents that
all A and B items are resolved.
