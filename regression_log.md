# Phase regression log

## 2026-07-15 — Phase 1 raw-object VBF/FSR provenance investigation opened

- **Trigger:** Phase 1 physics and critical re-reviews classify the absent raw-data NanoAOD/manifest for the requested VBF/FSR category as Category A.
- **Origin:** Phase 1 input/provenance contract for `fake_data_10fb.root`.
- **Initial evidence:** The owner production record shows that the supplied file is synthetic pseudo-data mixed from modified MC, not collision data. The flat fake-data tree has no Jet, Photon, or FsrPhoton branches.
- **Investigation status:** The documented modified-MC source chain reaches existing raw NanoAOD files, but it is not being treated as collision-data provenance. Exact raw-event and schema checks are being completed before the regression ticket is finalized.

## 2026-07-15 — Raw pseudo-data join result

- **Result:** A single serial scan of the 202 full-payload-resolved `nanoaod_modified` counterparts found 647/652 exact raw identities, no duplicates, and the required object branches in every one of 200 readable files.
- **Failure retained:** Two selected raw paths have no `Events` tree, leaving five named fake rows object-unavailable. They are documented in `phase1_strategy/outputs/raw_pseudodata_provenance_v2.json`; they have not been dropped or substituted.
- **Interpretation:** This is partial modified-MC pseudo-data recovery, not collision-data provenance. A complete VBF/FSR pseudo-data validation requires a documented replacement for the two unavailable raw counterparts or an explicit owner decision on the incomplete subset; the collision-data route remains unresolved.

## 2026-07-15 — Alternate source ruled out as a substitution

- **Search result:** There is no alternate modified-NanoAOD replica for either failed UUID. The documented original NanoAOD inputs contain all five rows and required object branches, but are pre-response inputs.
- **Compatibility result:** Each affected modified-H4l row differs from its nominal-H4l counterpart in reconstructed-lepton/four-lepton quantities; therefore the original raw objects cannot be silently mixed with the modified fake sample.
- **Root cause evidence:** The two `nanoaod_response` Condor error logs contain `TFile::Flush ... No such device`. Repair/re-run of these two production outputs is a narrow, reproducible blocker. A 647-row subset is engineering-only unless an owner explicitly redefines the pseudo-data validation scope.

## 2026-07-15 — Regression handoff prepared

- **Ticket:** `phase1_strategy/REGRESSION_TICKET.md` specifies a non-destructive two-file response regeneration path, exact source/config/command provenance, validation, impact cascade, and estimates.
- **Availability:** The original NanoAOD inputs, task configs, injector/source, and LCG_108 runtime are available. The injector needs the documented CMS+LCG environment; it is not runnable in the default shell because its GCC14/TBB runtime is not loaded.
- **Boundary:** Repairing the two modified-MC outputs can close only the pseudo-data object-coverage sub-blocker. It does not supply a collision-data source or authorize a data result.
