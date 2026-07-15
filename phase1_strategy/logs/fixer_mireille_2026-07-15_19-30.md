# Phase 1 fixer session log — Mireille

### 19:30 — Session start and evidence audit

Read the root and Phase 1 instructions, fixer contract, prompt, committed
strategy, all three reviews, arbiter verdict, plans, phase/root logs, relevant
methodology, and candidate conventions. The task is an ITERATE fix: preserve
the original strategy in a new `STRATEGY_v2.md`, create the required
`COMMITMENTS.md`, and resolve every arbiter Category A/B item with auditable
evidence rather than replacing VBF or FSR with a proxy.

### 19:30 — Raw-object provenance milestone

The local EOS inventory exposes raw MC at
`.../cms_opendata_2017_full_production/nanoaod_modified/` for all eleven fit
processes. A 27 MB VBF NanoAOD pilot file has an `Events` tree with 8,000
events and event identifiers, `Jet_*` kinematics/JEC inputs/b tags, and
`FsrPhoton_*` plus `Muon_fsrPhotonIdx`; it is paired by UUID with an H4l flat
output. The same production root exposes only `fake_data` flat products, not
a raw-data NanoAOD directory, so a data-level matched re-ntupling route is
not yet established and will be reported as unresolved unless provenance can
be demonstrated.
