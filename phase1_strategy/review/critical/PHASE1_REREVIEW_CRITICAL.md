# Phase 1 critical re-review — calibrated CMS Open Data $H\to ZZ^*\to4\ell$

**Verdict: ITERATE / BLOCKED — Phase 1 must not advance.**

The revision repairs the original strategy-contract defects in the fake model,
calibration, NN, fit preflight, method parity, comparison, and commitment
tracking plans. It does not make the requested data-level VBF/FSR analysis
executable. The supplied fake-data products still contain no raw jets or
photons, and no matching raw-data NanoAOD location or manifest is available.
Because the user requires a reconstructed-jet VBF category and FSR recovery,
this is a remaining Category A blocker, not a permissible downscope.

## Scope and evidence base

I read CLAUDE.md, prompt.md, the Phase 1 instructions, the critical-reviewer
contract, methodology/02-inputs.md, Phase 1 in methodology/03-phases.md,
methodology/06-review.md §§6.1--6.5.1, methodology/12-downscoping.md,
methodology/appendix-sessions.md, and the candidate conventions. The selected
detector-level simultaneous $(m_H,\mu)$ parameter fit is neither a closed-form
double-tag extraction (conventions/extraction.md:8-21) nor a particle-level
corrected spectrum (conventions/unfolding.md:6-30). The search convention is
not governing: the endpoint is neither a limit nor a discovery significance
(conventions/search.md:7-16).

I also read the original physics/critical/constructive reviews, the arbiter
adjudication, STRATEGY.md, STRATEGY_v2.md, COMMITMENTS.md, the raw-object
feasibility record, all preflight JSON artifacts, Phase 1/root experiment logs,
fixer logs, and the pilot/preflight source code. The corpus connector remains
unavailable as recorded in phase1_strategy/retrieval_log.md:3-7; I did not wait
for it or use its absence to waive a requirement.

| Independent current-state check | Direct result | Consequence |
|---|---|---|
| Supplied production-root inventory | The root has nanoaod_modified, MC products, and fake_data, but no raw-data NanoAOD directory. | Raw MC availability does not establish a raw-data route. |
| Supplied H4l fake-data schema | fake_data_10fb.root has h4lTree;1 and Metadata;1; h4lTree has **0** Jet_* and **0** Photon_*/FsrPhoton_* branches. | It cannot form a reconstructed-jet VBF category or recover FSR. |
| Supplied dilepton fake-data schema | fake_data_dilep_10fb.root has dilepTree and Metadata; its candidate tree has **0** Jet_* and **0** Photon_*/FsrPhoton_* branches. | It supports a lepton-only $Z$ control but cannot remedy VBF/FSR provenance. |
| Persisted VBF MC pilot | Direct read of raw_object_pilot_vbf_v1.root gives 1,198 rows, 1,134 in 105--140 GeV, 671 in-window rows with $\ge2$ cleaned jets, and 37 with a muon-associated FSR photon. These equal raw_object_pilot_vbf_v1.json. | The MC re-ntupling claim is real, but it is explicitly not a data feasibility result. |
| Pilot implementation | src/raw_object_feasibility.py:22-33,161-189 fixes a VBF MC raw input and MC flat input; the record says data provenance is unresolved (RAW_OBJECT_FEASIBILITY_v1.md:14-19,121-127). | No one-to-one fake-data-to-raw-data match exists. |

The raw-data conclusion is also stated by the revision itself:
STRATEGY_v2.md:15-19 says Phase 1 remains blocked, and COMMITMENTS.md:18-23
leaves raw-data provenance and reference-like FSR open rather than falsely
marking them downscoped.

## Finding-by-finding resolution audit

“Resolved at Phase 1” means that the original finding required an executable,
binding strategy contract and that v2 now supplies it. It does **not** mean
later implementation or data validation has occurred. The source labels below
are the arbiter consolidations in PHASE1_ARBITER.md:43-56, covering every
original Category A/B finding from the physics, critical, and constructive
reviews.

| Arbiter finding and original sources | Direct v2 evidence | Re-review resolution and status |
|---|---|---|
| **A1 raw-object VBF/FSR feasibility** — Physics A1; Critical A1; Constructive A1 | V2 freezes raw-jet and reference-like VBF requirements (STRATEGY_v2.md:85-111); MC branch coverage and 1,198/1,198 matching are documented in RAW_OBJECT_FEASIBILITY_v1.md:35-96. The direct supplied-data schema check above still finds no required objects. | MC feasibility and the no-proxy rule are repaired, but the required raw-data source, match, raw-systematic inputs, VBF/FSR validation, and information test do not exist. **A — unresolved; blocks advancement.** |
| **A2 predictive DY fake model** — Physics A2; Critical A2; Constructive A2 | Frozen $L/P/F$, 4P/3P1F/2P2F/2P2LSS definitions: STRATEGY_v2.md:137-159; transfer formula, prompt subtraction, covariance-correlated rate/shape modes: :161-188; independent hash-split/held-out closure, $p\ge0.05$, coverage 0.62--0.74, and failure rule: :200-217; separate reducible inventory: :190-198. | The missing contract is resolved at strategy level. V2 correctly says it cannot validate a fake-data fit without the raw-data route (:157-159). **Resolved at Phase 1; data execution is gated by A1.** |
| **A3 commitment ledger and systematic dispositions** — Critical A3; Constructive A3 | COMMITMENTS.md:13-164 tracks constraints/decisions, calibration, fakes/NN, fit, comparisons, and figures. The systematic programme is binding at STRATEGY_v2.md:428-456. | The absent-ledger defect is resolved. Luminosity and generator/shower availability are honestly held as Phase 2 evidence deadlines and conditional disposition paths (:441-445; COMMITMENTS.md:24-30), not silently omitted. **Resolved at Phase 1; audit these open evidence deadlines before Phase 2 closes.** |
| **A4 numerical 10 fb$^{-1}$ lepton-scale sanity contract** — Physics B2; Critical B2; Constructive A4 | The JSON contains 4,508,686 $ee$ and 6,284,267 $\mu\mu$ candidates, with $2.57\times10^{-5}$ and $2.07\times10^{-5}$ scale-statistical floors and 3.22/2.59/2.07 MeV channel mass floors. V2 fences them as floors, not fit results (STRATEGY_v2.md:280-314). Joint covariance, a non-overlapping ZFitModel envelope, $\pm0.05\%/\pm0.10\%$ injected offsets, 500 toys, and profiled scale impact are fixed at :258-307. | The old 0.38 GeV mass-statistics proxy is no longer passed off as a scale test. **Resolved at Phase 1; no calibration or fit result has yet been validated.** |
| **A5 HIG-16-041 statistical-method parity** — Critical A4 | The parity matrix contrasts binned $m_{4\ell}$ with HIG-16-041 1D/2D/3D likelihoods, constrained mass, per-event resolution, categories, and data-derived fakes (STRATEGY_v2.md:406-418). An independent constrained-mass/resolution/kinematic-score pseudoexperiment test must agree within 0.2 combined standard deviations and pass coverage for both POIs (:420-426; COMMITMENTS.md:115-119). | Unsupported equivalence became an explicit non-equivalence plus a falsifiable cross-check. **Resolved at Phase 1; parity has not passed yet.** |
| **B6 calibration definition and data/MC convention** — Physics B1; Critical B1; Constructive B1 | Frozen trigger/ID/isolation/vertex/mass control selection, BW $\otimes$ double-Crystal-Ball plus background model, alternatives, and GoF: STRATEGY_v2.md:219-241. Separate absolute data/MC transformation to common $M_Z$, full rebuild, and joint covariance: :243-278. | The template-convention ambiguity, common-$M_Z$ correlation, subdivision gate, and injected closure requirements are explicit. **Resolved at Phase 1.** |
| **B7 angular-NN decision and score transfer** — Physics C2; Critical B3; Constructive B2 | Inputs, exclusions, deterministic split, and all fit backgrounds are fixed at STRATEGY_v2.md:316-337. Five quantitative quality gates and deterministic revision handling are at :339-367. | The score-level, mixture-aware decision contract exists; a quiet cut-based substitution is prohibited. **Resolved at Phase 1.** |
| **B8 simultaneous-fit sparse-category/coverage preflight** — Physics B3; Critical B4; Constructive B3 | Six flavour$\times$VBF/non-VBF cells, 1 GeV nominal/0.5--2 GeV alternative bins, gamma constraints, 120--130 GeV 0.1 GeV grid, and leave-node closure: STRATEGY_v2.md:369-390. Per-cell identifiability and 500 independent toys at $(124,1)$, $(125,0.8)$, $(125,1.2)$, $(126,1)$ with pull, coverage, boundary, GoF, and resolving-power gates: :392-404. | The old qualitative sparse-category promise is falsifiable and preserves VBF when score bins are sparse. **Resolved at Phase 1.** |
| **B9 comparisons/comparability/provenance** — Physics B4; Critical B5; Constructive B4 | STRATEGY_v2.md:458-472 provides targets, source provenance, equivalence and overlap treatment, and allowed comparison statistics for $m_H$, $\mu$, $m_{4\ell}$, VBF, NN, and $Z$ controls. It requires a descriptive difference rather than an independence-implying pull without overlap covariance, and explains why a world-$\mu$ comparison is not meaningful. | The comparison plan is now scoped and avoids treating a calibration input as independent validation. **Resolved at Phase 1.** |

## Systematic and parity disposition check

All 15 source groups at STRATEGY_v2.md:434-450 now have binding treatment,
including the formerly open luminosity, pileup/JEC/JER/b tagging, signal
acceptance/production/PDF/QCD, shower/UE, separate $q\bar q$/$gg$ ZZ, DY
transfer covariance, and top/conversion/electroweak rows. V2 requires every
variation to propagate through reconstruction, categories, templates, and the
simultaneous fit (:428-432), and prohibits an arbitrary flat fake error. That
meets the Phase 1 systematic-disposition requirement.

This is not evidence that exact-dataset luminosity, usable generator weights,
or raw-object VBF variations already exist. Those remain open ledger items and
must not become zero, borrowed, or silently absent systematics. Similarly, the
method-parity matrix is an adequate Phase 1 plan, not a passed parity test; the
raw-data gap blocks its required VBF component in a final data analysis.

## Remaining/new findings

| ID | Classification | Finding and required disposition |
|---|---|---|
| R1 | **A — must resolve** | **Matching raw-data NanoAOD/manifest is absent.** Before Phase 1/Phase 2 advancement, obtain a source one-to-one matching fake_data_10fb.root on $(run,lumi,event)$ and preserving selected jets, JEC/JER-relevant inputs, DeepJet inputs, FSR/generic photons, and lepton inputs. Re-run the data pilot, then execute the VBF/FSR, data fake-control, and category-information gates. Do not substitute $p_T^{4\ell}$, candidate multiplicity, or a flat-tree proxy. |
| R2 | **C — audit clarity** | STRATEGY_v2.md:21-23 says finding dispositions are in outputs/FIX_RESOLUTION_v1.md, but that file is absent. This re-review supplies the required audit, but the stale link should be corrected or replaced before the next revision. |

There are **no remaining Category B findings** in the Phase 1 strategy contract.
The Category A raw-data gap is sufficient to prevent a PASS and means the
otherwise repaired data-control portions of the fake, FSR, VBF, and final-fit
plans cannot yet be claimed as executed.

## Regression and completion audit

No Phase 1 regression ticket is appropriate: there is no earlier analysis phase
or result to regress. There are no Phase 1 figures, fit outputs, closure
results, or uncertainty breakdowns to validate visually or numerically, and
this review does not pretend unrun tests passed. The only valid next state is
to retain the full requested scope, obtain the external raw-data provenance,
and run the committed pilots and validation chain.

**Final critical decision: ITERATE / BLOCKED.** The revision is an honest,
substantially improved strategy, but it cannot advance toward the requested
VBF/FSR data analysis until raw-data NanoAOD provenance is supplied and
verified. A Phase 4b human-gate permission does not relax this prerequisite.
