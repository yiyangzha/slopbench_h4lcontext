# Phase 1 experiment log

## 2026-07-15 — Strategy session initialized

- Read the root user objective and Phase 1 instructions before analysis drafting.
- Created `plan.md` first, as required.  The planned primary technique is a simultaneous multi-category profile-likelihood template fit to calibrated $m_{4\ell}$, not an unfolded particle-level spectrum or a closed-form double-tag extraction; applicability of both candidate conventions will be stated and audited in the strategy.
- The retrieval corpus is not currently exposed in the available tool set.  The executor will log attempted corpus-access verification and use primary official CMS, HEPData, and PDG sources rather than treating remembered values as numeric inputs.

## 2026-07-15 — Constructive review session started

- Began an independent evidence-based audit of the committed Phase 1 strategy without modifying it.  The review scope explicitly covers raw-object VBF recovery, four-lepton-rest-frame angular NN feasibility, Z-control calibration and propagated residual scale uncertainty, DY+jets fake transfer and closure, simultaneous-fit coverage, reference/PDG comparison metrics, and 10 fb$^{-1}$ scale sanity.
- The corpus-connector limitation already recorded in `retrieval_log.md` will be treated as an uncertainty in source verification, not as evidence that the strategy requirements are met.

## 2026-07-15 — Constructive review completed

- Wrote `review/constructive/PHASE1_CONSTRUCTIVE_REVIEW.md` with an ITERATE verdict. The review finds the planned angular NN viable at candidate level, but requires a demonstrable raw-object VBF/FSR route, predictive DY-fake transfer and closure, `COMMITMENTS.md`, and a numerical 10 fb$^{-1}$ lepton-scale sanity test before Phase 1 can pass.

## 2026-07-15 — Arbiter adjudication completed

- Wrote `review/arbiter/PHASE1_ARBITER.md` after independently checking the strategy, all reviews, the applicable conventions and review methodology, the missing commitment ledger, and the current ntuplizer object availability.
- Verdict: **ITERATE**. The Phase 1 strategy must resolve raw-object VBF/FSR feasibility, predictive DY-fake modelling, the commitment/systematics ledger, the numerical lepton-scale sanity contract, and published-method parity; no Phase 1 regression trigger or escalation is applicable yet.

## 2026-07-15 — Critical re-review completed

- Wrote `review/critical/PHASE1_REREVIEW_CRITICAL.md` after reviewing the original panel/arbiter artifacts, `STRATEGY_v2.md`, `COMMITMENTS.md`, raw-object feasibility record, preflight JSON, source code, and experiment/session logs.
- A read-only direct check finds zero `Jet_*`, `Photon_*`, and `FsrPhoton_*` branches in both supplied fake-data candidate products, while the persisted raw-object pilot verifies an MC-only route. The revised strategy contracts resolve the original planning defects, but the missing matching raw-data NanoAOD/manifest remains Category A and blocks Phase 1 advancement without a proxy/downscope.

## 2026-07-15 — Pseudo-data raw-object regeneration started

- Under `REGRESSION_TICKET.md`, created two immutable response-repair configs
  with the exact original NanoAOD inputs, seed, branch map, scale, and smear.
  Their only configuration changes are a new versioned output root/suffix and
  fresh validation-summary directories; both retain `overwrite: false`.
- The documented CMS + LCG_108 `nanoaod_response --dry-run --threads 1`
  accepted both configs and planned the two absent versioned outputs. The
  ticketed corrupted files and all existing production outputs remain intact.
- The generic response staging helper rejects the immutable per-file task
  contract because it expects a directory input. Its exact error was captured;
  widening to a parent directory would stage unrelated files and was not done.
  A literal one-file JDL using the documented worker wrapper passed
  `condor_submit -dry-run` and submitted the DY test as cluster `952603`.
