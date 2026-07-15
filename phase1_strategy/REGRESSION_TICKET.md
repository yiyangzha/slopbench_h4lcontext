# Regression ticket — Phase 1 raw-object provenance and corrupted modified NanoAOD

**Status:** Open production/provenance blocker. No downstream physics result has
been accepted, so there is no numerical result to retract.

**Origin phase:** Phase 1 strategy/provenance contract, with the immediate
technical fault in the upstream modified-MC production.

**Fix owners:** A production executor for the two corrupted NanoAOD response
outputs; then the Phase 1 strategy owner for scope/review updates. Do not
silently change the analysis endpoint from collision data to pseudo-data.

## 1. Origin and reviewer finding

The Phase 1 physics re-review says not to substitute a flat-tree proxy or
advance to a data fit until matching raw-data provenance exists
(`review/physics/PHASE1_REREVIEW_PHYSICS.md`, A1). The critical re-review
independently records that the supplied fake trees have no jets or photons and
that the production root has no raw-data NanoAOD directory
(`review/critical/PHASE1_REREVIEW_CRITICAL.md`, A1). `STRATEGY_v2.md` therefore
correctly freezes a no-proxy VBF/FSR requirement and requires a one-to-one raw
data match before a VBF data fit.

The investigation refines the provenance, without weakening that conclusion:

- `fake_data/fake_data_10fb.root` is **synthetic pseudo-data**, not collision
  data. The producer contract at
  `/eos/home-y/yiyangz/codex/slopbench_code/AGENTS.md` and
  `results/.../fake_data/mixer_10fb_summary.json` show that it is a mixture of
  400 `h4l_mc_modified` inputs. All eight supplied fake-data summaries likewise
  identify modified-MC mixer inputs.
- Its 652-row `h4lTree` has no `Jet_*`, `Photon_*`, or `FsrPhoton_*` branches.
  This is corroborated by `outputs/RAW_OBJECT_FEASIBILITY_v1.md` and the two
  reviews.
- Every fake row can be resolved to one documented modified-H4l source entry,
  but not with `(run,lumi,event)` alone: 18 identifiers have two or three
  candidates. Comparing all 111 copied output branches resolves each one
  uniquely. A persistent source-file/source-entry token is required for future
  production.
- The resolved sources point to 202 `nanoaod_modified` files. A single serial,
  read-only scan found 647/652 unique raw identities and the required object
  schema in all 200 readable files. The remaining five rows map to two ROOT
  files with no keys and no `Events` tree.

The full evidence is in
`outputs/raw_pseudodata_provenance_v2.json` and
`outputs/raw_pseudodata_alternate_source_check_v1.json`.

## 2. Exact fault, affected rows, and source provenance

| Sample | Corrupted modified raw path | Documented original NanoAOD input | Affected fake identities |
|---|---|---|---|
| DYJetsToLL | `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/nanoaod_modified/DYJetsToLL/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear.root` | `/eos/opendata/cms/mc/RunIISummer20UL17NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F.root` | `(1,390917,97729202)`, `(1,927325,231831246)`, `(1,672725,168181225)` |
| GluGluToHToZZ | `/eos/cms/store/group/phys_bphys/trigger/yiyangz-contact/h4l/cms_opendata_2017_full_production/nanoaod_modified/GluGluToHToZZ/2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear.root` | `/eos/opendata/cms/mc/RunIISummer20UL17NanoAODv9/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A.root` | `(1,559,215650)`, `(1,1739,670997)` |

Both original inputs are readable, contain their listed identities exactly once,
and have the identity, jet, FSR-photon, generic-photon, and lepton branch
groups. They are **not** valid replacements: the pseudo-data was made from the
scale/smear-modified chain. For every affected row the modified and nominal H4l
payloads differ in `mZ1`, `mZ2`, `m4l`, four-lepton kinematics, and all four
selected-lepton `pt/eta/phi` values (19 changed fields for four rows and 37 for
`(1,927325,231831246)`). Mixing original raw objects with the modified flat
candidate would be event-inconsistent.

The per-file response jobs record the likely failure directly:

- `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/nanoaod_response/runs/20260713T120732Z_2133781/logs/DYJetsToLL__280000__3183559F-C042-5944-84D4-63F2E7CF0D0F.949114.240.err`
- `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/nanoaod_response/runs/20260713T120732Z_2133781/logs/GluGluToHToZZ__2430000__78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A.949114.2.err`

Each contains `SysError in <TFile::Flush> ... No such device` for the named
`_scale_smear.root` output. The task-execution summaries nevertheless report
exit code zero, so successful process exit is not an adequate validation.

## 3. Regression classification

The operative trigger is the reviewers' Category-A A1 raw-object prerequisite:
the requested reconstructed-jet VBF and FSR programme cannot be formed from the
given flat fake product and must not be replaced by a proxy. The present
investigation additionally finds a concrete upstream output-integrity failure
in the only partial pseudo-data raw route.

None of the numerical §6.7 triggers (data/MC disagreement, failed closure,
GoF, parameter pull, bin exclusion, or systematic double counting) has yet
occurred: no data fit, closure, or accepted result exists. This ticket is
nevertheless required under the §6.7 procedure and upstream-improvement
cascade because a mandatory Phase 1 input/provenance condition is invalid and
would otherwise make all later raw-object consumers stale or misleading.

## 4. Impact trace and required cascade

| Phase / consumer | Affected work | Must re-run after a successful repair | Can remain unchanged now |
|---|---|---|---|
| Phase 1 strategy | A1-data, A2 FSR, D5 VBF, D7 staged-data wording, and feasibility evidence | Update provenance labels and append the repaired 652/652 evidence; re-review A1. Preserve the distinction between pseudo-data and collision data. | The no-proxy requirement, nominal-MC feasibility pilot, fit/NN design, and theory-reference plan. |
| Phase 2 inventory | Raw-input, source-kind, data-luminosity, and Z-control provenance inventories | Re-inventory the repaired pseudo-data mapping and explicitly record that the fake H4l and dilepton products are modified-MC mixtures. If the endpoint remains collision data, obtain a separate real-data product and raw manifest. | Nominal-MC cross sections and unrelated MC inventory entries. |
| Phase 3 reconstruction | Raw-object producer, FSR recovery, VBF variables/migrations, calibration rebuild, fake controls, NN data validation | Re-run all consumers of the repaired raw objects and demonstrate full source-token/identity coverage before any pseudo-data VBF/FSR result. | Pure strategy definitions and MC-only work proven independent of these two files. |
| Phase 4a | MC/pseudo-data validation | Re-run only if it consumes this fake-data product or regenerated H4l/mixer output. Record any 647-row engineering check as incomplete. | Independent MC-only toys that never read the pseudo-data product. |
| Phase 4b / 4c | Deterministic data gate and physics fit | Remain blocked. Modified-MC pseudo-data is not collision data; repair alone does not authorize a data fit or observed-data claim. | Nothing may be promoted to a data result. |
| Phase 5 / AN | VBF/FSR, calibration, NN, and final-fit figures/tables | Generate only after the preceding gates pass; label any interim 647-row study as engineering-only. | Existing planning diagrams and references that make no completed-result claim. |

## 5. Required scope decision

The fixer must not make this decision implicitly. The analysis owner must choose
one of the following endpoints before a physics result is advanced:

1. **Synthetic pseudo-data validation track.** Repair the two modified raw
   outputs and prove an exact 652/652 source-token/identity match. All results
   remain labelled modified-MC pseudo-data; no collision-data calibration,
   VBF result, or observed fit is claimed.
2. **Requested collision-data measurement.** The producer supplies a real
   collision-data H4l/dilepton product plus a raw-data NanoAOD manifest and
   deterministic source tokens. That is new input production, not a repair of
   the two MC files. It is required for the user-requested data VBF/FSR and
   calibration programme.

Until an owner explicitly approves a reduced engineering scope, a 647-row
subset is not scientifically sufficient: its missing 0.77% is process-specific
(DY and ggH) and can affect categories and composition. It may support a
clearly labelled implementation diagnostic only, never a complete validation
or physics result.

## 6. Non-destructive repair handoff for the pseudo-data track

### Available inputs and exact transformation provenance

All data/configuration/source inputs are locally readable and the documented
runtime is available:

| Item | Location / result |
|---|---|
| Original NanoAOD inputs | The two exact EOS paths in Section 2 are readable (1,283,487 DY and 827,970 ggH `Events`) and contain all five identities. |
| Immutable per-file response configs | `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/nanoaod_response/tasks/DYJetsToLL/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F/task.json` and `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/nanoaod_response/tasks/GluGluToHToZZ/2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A/task.json`. They fix seed `20260713`, tree `Events`, the branch map, saturating-downward scale, and pT/eta/phi smear. |
| Injector/source | `/eos/home-y/yiyangz/codex/slopbench_code/benchmark/nanoaod_response` and `/eos/home-y/yiyangz/codex/slopbench_code/benchmark/nanoaod_response.cpp`; default configuration is `/eos/home-y/yiyangz/codex/slopbench_code/benchmark/configs/lepton_response.json`. |
| Worker command | `/eos/home-y/yiyangz/codex/slopbench_code/benchmark/production/run_nanoaod_response_task.sh` sources CMS plus `LCG_108/x86_64-el9-gcc14-opt`, then runs `nanoaod_response --config TASK_CONFIG --threads 1`. |
| Runtime check | The injector's `--help` succeeds after sourcing `/cvmfs/cms.cern.ch/cmsset_default.sh` and `/cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc14-opt/setup.sh`. A bare shell lacks its GCC14/TBB runtime, so the documented setup is mandatory. |
| Downstream H4l configs | `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/h4l_ntuplizer/tasks/h4l_ntuplizer/modified/a6bc428b076e3a631098bb37d6f68bb7604a79d66332be1254544d2e66fe5e93/DYJetsToLL/280000/3183559F-C042-5944-84D4-63F2E7CF0D0F_scale_smear/task.json` and `/eos/home-y/yiyangz/codex/slopbench_code/results/cms_opendata_2017_full_production/condor/h4l_ntuplizer/tasks/h4l_ntuplizer/modified/a6bc428b076e3a631098bb37d6f68bb7604a79d66332be1254544d2e66fe5e93/GluGluToHToZZ/2430000/78FC46EA-1FB7-3B4B-8BDC-A2F62B73F44A_scale_smear/task.json`. |

The existing response task configs set `overwrite: false`. Do **not** delete,
overwrite, or rename the corrupted ROOT files. Create two new immutable repair
configs whose only intentional changes are a fresh output root/suffix and fresh
validation-summary directory; retain the original input, seed, branch map,
scale, smear, and one-file-per-task topology. First use `--dry-run`, then stage
a single-file test under the documented Condor workflow before running the two
repair tasks.

The command provenance for the executor is:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc14-opt/setup.sh
/eos/home-y/yiyangz/codex/slopbench_code/benchmark/nanoaod_response \
  --config /path/to/new_immutable_repair_task.json --threads 1
```

Use the existing wrapper/staging workflow where possible so the Condor record,
binary, input, output, and validation sidecars are preserved. The phase fixer
is authorized by this ticket to prepare/execute only after the analysis owner
authorizes external production; this investigation did not regenerate files.

### Validation required before accepting a repair

1. New files open with a non-zombie `Events` tree, with 1,283,487 DY and
   827,970 ggH entries; the original-to-repaired identity mapping is complete
   and unique.
2. Persist branch/schema checks for identity, `Jet_*`, raw/JEC/JER/DeepJet
   inputs, `FsrPhoton_*`, `Photon_*`, muons, and electrons. Verify preservation
   of all non-lepton event content and expected response only in the configured
   lepton fields.
3. Require an explicit write/close validation, not just exit status; inspect
   the response JSON/histogram sidecars and reject any `TFile` write error.
4. Re-run the two modified-H4l tasks to fresh paths. Compare every surviving
   source entry, especially the five fake rows, against the current
   `h4l_mc_modified` payload. If they are not event-level compatible, build a
   new versioned mixer manifest and a new pseudo-data file rather than altering
   `fake_data_10fb.root`.
5. Persist a fake-row-to-raw mapping with source file and source entry. Do not
   use `(run,lumi,event)` alone because 18 fake identities collide upstream.
6. Repeat the complete raw scan and require 652/652 exact, unique identities,
   zero unreadable selected raw files, and all required branch groups before
   VBF/FSR re-ntupling.
7. Re-run the committed FSR preservation, VBF yield/efficiency/purity,
   JEC/JER/b-tag migration, calibration rebuild, fake-control, angular, and
   closure gates that consume the repaired objects. Re-review Phase 1 A1/A2
   before proceeding.

## 7. Estimate

For the synthetic pseudo-data repair, estimate **3–6 active agent-hours** to
make immutable task configs, stage two jobs, validate raw/H4l/mapping outputs,
and update the evidence; EOS/Condor queue time is additional. If the repaired
H4l payload is not event-level compatible and a new mixer/product must be made,
allow **6–12 active agent-hours** plus re-running all downstream consumers.
Producing an actual collision-data route is external producer work and cannot
be credibly estimated from this workspace.

## 8. Acceptance and closure

Close this ticket only after the selected endpoint has passed the relevant
requirements above, Phase 1 has been re-reviewed, and every downstream consumer
has either been re-run or explicitly shown independent. A repaired modified-MC
route closes only the pseudo-data sub-blocker; it does not close the separate
collision-data A1 requirement.
