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
