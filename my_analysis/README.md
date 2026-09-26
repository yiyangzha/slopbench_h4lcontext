# H -> ZZ* -> 4l evaluation submission (task h4l_ntuple, EVAL_CONTRACT v0.4)

    ./run.sh <dataset_dir> <output_dir>

writes `RESULT.json` and `MODEL.json` (the fitted binned likelihood as a pyhf workspace at the fitted m_H, buildable with
`pyhf.Workspace(spec).model()`), plus `MODEL_meta.json`, `summary.txt`, `results_full.json` and `logs/run.log`.
The environment is `pixi.toml` / `pixi.lock` (ROOT for the large input files, numpy/scipy/iminuit/uproot); `run.sh`
uses an already active environment with ROOT and those packages, else `pixi run --frozen`.  `EVAL_CONTRACT.md` is the
contract with the per-flavour calibration block used here.

Method: the main analysis of this repository (`analysis_v3/`; AN in `deliverables/`) with the same methods and
parameters in every step; the only simplification is the tag-and-probe nominal fit model; otherwise only equivalent
faster implementations and what the environment forbids (MELA):
1. `src/h4l_reader.cpp` (ROOT, compiled at run time, parallel over file chunks): the analysis trigger OR and a good
   vertex, AN-16-442 loose leptons with ghost cleaning and FSR photons, the Z -> ll calibration pairs (with the
   relative momentum errors and the charge x eta/phi slices of the legs), tag-and-probe pairs, loose-lepton event
   records (data; DY and ttbar control and Z + 1 lepton candidates; ZZ and signal candidates) and the generator
   bookkeeping (genWeight, genEventSumw).
2. MC normalization: genWeight-based, 1000 L sigma_eff genWeight / sum of genEventSumw of the files read (the
   count-based n_selected / n_preselection is recorded as a cross-check); half of the DY, ttbar and ZZ files are read,
   the ZZ control and Z + 1 lepton rows (prompt subtraction of the fake rates) come from a fifth of those.
3. Lepton calibration per flavour on Z -> ll (`h4l_eval/calibration.py`): the factorized scale/smear model, the two
   category families, the event-level template (frozen pair deviates carry the extra smearing, a fixed kernel of 0.3 %
   relative width, likelihood groups of at least 1.2 GeV; in histogram form with binned deviates), the response pass
   with the weighted MC as the data role, the joint least squares with outlier rejection and the iteration.
4. Per-event mass uncertainty (`h4l_eval/lam.py`): lambda measured on the calibrated Z -> ll pairs by the simultaneous
   BW (x) DCB fits in bins of the predicted uncertainty, data lambda for data and MC lambda for MC.
5. Tag-and-probe SF of the full single-lepton selection (`h4l_eval/tnp.py`): nominal simplified (calibrated MC template,
   CMSShape background), fit-model systematic from the stand-alone DSCB signal and the Bernstein background applied to
   data and MC alike.
6. Final selection (`h4l_eval/selection.py`): FSR-subtracted isolation on the calibrated leptons, AN tight leptons, ZZ
   candidates with the smart cut, the Z1 kinematic refit and the per-event mass uncertainty; the control regions
   (2P2F, 3P1F, same-sign) and the Z + 1 lepton rows as the main analysis's h4l_select (candidate choice: Z1 closest to
   m_Z, then the largest Z2 scalar pT sum; no MELA).
7. Z+X (`h4l_eval/zx.py`): OS and SS methods with prompt-subtracted fake rates (SS: conversion correction), the MC
   closure systematic on the DY and ttbar MC, bootstrap statistics, the combination and the Landau + exponential shape.
8. Signal model and systematics (`h4l_eval/templates.py`): ggH + VBF DCB per channel, non-resonant VH, rest-frame
   scaling with m_H, A x eff(m_H) with the morphing nuisance, YR4 sigma_eff; the scale and resolution magnitudes of
   the main analysis's make_systematics.
9. Binned likelihood m4l_refit (105-140 GeV, 0.5 GeV) x 5 D_mass bins x 3 final states (`h4l_eval/model.py`), mu and
   m_H floating (m_H in [110, 140] GeV), nuisances: lepton scale/resolution and morphing (histosys), efficiencies and
   Z+X (normsys); MINOS intervals, stat with the nuisances fixed, syst from the impacts, toys for the GoF p-value (5 GeV
   bins) and the coverage.

Calibration in RESULT.json: per flavour under the registered names (`scale_shift`, `smear`, `sel_eff`, each
`{"muon": {value, unc}, "electron": {value, unc}}`), the data-weighted averages over the Z -> ll calibration-sample
leptons.
