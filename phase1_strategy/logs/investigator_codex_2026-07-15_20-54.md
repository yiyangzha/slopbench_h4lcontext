# Investigator session log — Codex — 2026-07-15 20:54 CEST

### 20:54 — Session start and regression scope

Investigator assigned to trace the raw-object VBF/FSR provenance blocker without modifying strategy or analysis code. Read the root and Phase 1 instructions, the investigator contract, the Phase 1 re-reviews, `STRATEGY_v2.md`, `COMMITMENTS.md`, `RAW_OBJECT_FEASIBILITY_v1.md`, prompt, relevant logs, and the regression methodology.

### 20:54 — Fake-data provenance checkpoint

All eight JSON records under the supplied `fake_data` directory are dilepton mixer summaries and cite only modified MC inputs. The owner production records additionally prove that `fake_data_10fb.root` was created by `benchmark/h4l_mixer` from 400 `h4l_mc_modified` files, not collision data. The 652 fake H4l rows all occur in that documented source list; 18 `(run,lumi,event)` identifiers collide across source files, but each collision is resolved to exactly one source by an exact comparison of all 111 copied H4l branches.

### 20:54 — Raw-NanoAOD scan checkpoint

The resolved pseudo-data rows map by filename to 202 existing `nanoaod_modified` files. A complete identity/schema scan requires 210,195,551 raw `Events` rows; a 2,083,510-event pilot reads the three identity branches in 4.5 s, giving an approximately eight-minute serial estimate. Two inline parallel-launch attempts failed before raw coverage: Python 3.14 process workers cannot import an inline function, and the thread retry terminated without completing a file. No files outside this log were changed. A stable read-only batched scan remains in progress.

### 21:16 — Serial raw join completed

The single retained, serial Uproot scan read 207,806,094 `Events` rows across the 202 full-payload-resolved raw counterparts. It found 647 exact `(run, luminosityBlock, event)` matches for the 652 fake rows, with no duplicate raw matches. All 200 readable raw files contain the required identity, jet, FSR-photon, photon, and electron/muon four-vector branches; this was a branch-presence check only.

Two selected raw files open as empty ROOT files (no keys and no `Events` tree):

- `nanoaod_modified/DYJetsToLL/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear.root` — 3 target rows.
- `nanoaod_modified/GluGluToHToZZ/2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear.root` — 2 target rows.

The five unrecoverable fake identities are `(1,559,215650)`, `(1,390917,97729202)`, `(1,927325,231831246)`, `(1,672725,168181225)`, and `(1,1739,670997)`. The persistent result is `phase1_strategy/outputs/raw_pseudodata_provenance_v2.json`. The outcome proves a partial modified-MC pseudo-data raw-object route only; it does not resolve the collision-data prerequisite.

### 21:20 — Alternate-source and production-failure check

No alternate `nanoaod_modified` basename or replica exists for either UUID. The documented original NanoAOD inputs do contain the five exact identities and the complete object schema, but they are pre-response inputs rather than event-faithful replacements. Comparing the corresponding modified and nominal H4l rows shows 19 changed reconstructed-lepton/four-lepton branches for four rows and 37 for `(1,927325,231831246)`.

The modified-output Condor error logs record `SysError in <TFile::Flush> ... No such device` for both paths, explaining why the nominally large ROOT files have no keys. Therefore a 647-row subset may support only a clearly labelled engineering diagnostic with explicit owner approval; it is not a complete pseudo-data validation or physics result. The narrow blocker is production repair/regeneration of those two modified NanoAOD outputs, followed by event-level revalidation. Evidence is persisted in `phase1_strategy/outputs/raw_pseudodata_alternate_source_check_v1.json`.

### 21:25 — Reproducible repair handoff and ticket

The exact original NanoAOD inputs, immutable per-file task JSONs, C++ injector/source, response configurations, and downstream H4l task JSONs are locally readable. The injector requires the documented CMS + `LCG_108/x86_64-el9-gcc14-opt` environment: it fails in the default shell because GCC14/TBB libraries are absent there, but its `--help` succeeds after sourcing that setup. No regeneration was performed.

`REGRESSION_TICKET.md` now gives the exact two corrupt paths, their five rows, original inputs, Condor write-error logs, task/config/command provenance, non-destructive new-output rule, validation gates, cascade, and estimates. The ticket explicitly keeps the collision-data route blocked even if the synthetic pseudo-data repair succeeds.
