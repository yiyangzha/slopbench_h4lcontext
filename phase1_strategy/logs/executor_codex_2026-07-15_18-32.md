# Phase 1 executor session log — Codex — 2026-07-15 18:32 CEST

### 18:32 — Recovery and evidence milestone
Resumed from the committed Ada plan without changing any checkpoint file. Read the root and Phase 1 instructions, executor contract, applicable methodology, prompt, both candidate conventions, logs, retrieval record, and the predecessor log. The experiment-corpus MCP tools remain unavailable, so the documented fallback is being used: official CMS HIG-16-041/JHEP 11 (2017) 047, its HEPData record, and the PDG listing.

### 18:32 — Supplied-data capability check
Inspected the flat ntuplizer sources and opened the supplied `fake_data_10fb.root` with `uproot`. `h4lTree` has 652 entries and 111 scalar branches, containing four-lepton and lepton quantities but no jet, b-tag, or missing-transverse-momentum fields. The Phase 1 strategy will therefore bind re-ntupling that preserves those objects before VBF categorization, rather than falsely claiming that VBF can be recovered from the current flat tree.
