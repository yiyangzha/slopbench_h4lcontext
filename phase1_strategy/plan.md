# Phase 1 execution plan

## Scope and success criteria

This Phase 1 deliverable will define a reproducible strategy for a calibrated CMS Open Data $H\to4\ell$ measurement.  It will bind the later phases to: separate $Z\to ee$ and $Z\to\mu\mu$ scale calibrations, a VBF category, a four-lepton-rest-frame angular NN discriminator, DY+jets as the reducible/fake model, and a simultaneous multi-category $m_{4\ell}$ likelihood fit for $m_H$ and $\mu$.  It will also define the deterministic 10% data validation and the 10 fb$^{-1}$ full-data result sequence.

Success means that `outputs/STRATEGY.md` is self-contained and includes all Phase 1 mandatory deliverables, explicit [A]/[L]/[D] labels, a binding systematic enumeration, a 2--3-analysis reference table with extracted numerical targets, a concrete two-approach selection/MVA comparison plan, and six named flagship figures.  `COMMITMENTS.md`, phase and root experiment-log entries, a retrieval record, and an incremental executor log will provide the audit trail.

## Evidence to collect before drafting

1. Query the experiment retrieval corpus for CMS $H\to ZZ\to4\ell$ measurements, their systematic programs, and cross-experiment or world-average comparisons.  For every identified reference, retrieve the relevant paper sections.  If the corpus is unavailable, record the failed tool availability and use primary official CMS, HEPData, and PDG sources with direct citations instead.
2. Inspect `ntuplizer/h4l_ntuplizer.cpp` and `ntuplizer/dilep_ntuplizer.cpp` to establish which reconstructed leptons, jets, weights, and four-lepton-rest-frame angular quantities are actually available.  Inspect a small ROOT metadata/schema sample only after estimating inputs; do not assume branches.
3. Retrieve the official CMS HIG-16-041 / JHEP 11 (2017) 047 analysis, its mass-fit window, its published mass and signal-strength results, and the PDG Higgs and $Z$ reference values used only as cited validation targets.
4. Determine whether the chosen simultaneous template/profile-likelihood parameter extraction is governed by either candidate convention.  Record why `extraction.md` and `unfolding.md` apply or do not apply, then enumerate any applicable requirements without silently omitting a source.

## Artifact structure to produce

- `phase1_strategy/outputs/STRATEGY.md` — primary Phase 1 strategy.
- `phase1_strategy/COMMITMENTS.md` — binding [D], [A], [L], systematic, validation, comparison, and flagship-figure commitments.
- `phase1_strategy/experiment_log.md` and root `experiment_log.md` — append-only decisions and evidence.
- `phase1_strategy/retrieval_log.md` — corpus queries, responses, and primary-source fallback record.
- `phase1_strategy/logs/executor_<session>_<timestamp>.md` — incremental session notebook.

## Work sequence and checks

1. Build the source and ntuple evidence record; verify every numerical input is traceable to a source or explicitly deferred to later data inspection.
2. Draft the selection, calibration, MVA, categorization, fake-model, and fit plan; retain two qualitatively different selection approaches for Phase 3 with a shared quantitative comparison metric.
3. Write systematic sources as propagated nuisance/variation plans, including scale-fit statistical, model, binning, and residual components correlated across categories and channels as justified.
4. Define deterministic split rules, closure/stress tests, control/sideband regions, reference comparisons, expected/full-data staging, and planned visual evidence.
5. Self-audit against Phase 1 requirements and the candidate conventions; verify the six flagship figures and all commitments are discoverable in both strategy and commitments files.
6. Append logs, run markdown/reference/link consistency checks, inspect the git diff, and commit with a conventional Phase 1 commit.

## Explicit assumptions to test rather than silently adopt

- The supplied ROOT ntuples contain the reconstructed inputs required for the specified angular NN and VBF selection; branch-level confirmation is required before Phase 2/3 implementation.
- The effective cross sections in `prompt.md` are the supplied normalization inputs and the `Metadata.nEvents` sums are the only allowed generated-event denominators.
- A published CMS mass-fit window and numerical comparison targets can be recovered from official primary sources; if not, the inability will be logged rather than replaced with recalled numbers.
