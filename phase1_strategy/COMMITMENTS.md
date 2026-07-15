# Phase 1 commitment ledger — strategy v2

This is the mandatory, append-only Phase 1 tracking ledger. `[x]` means the
item is resolved at the stated phase; `[ ]` means it remains binding and open;
`[D]` means formally downscoped with a documented reason. There are no `[D]`
entries in this version: unavailable raw-data provenance is an open blocker,
not permission to downscope the requested VBF/FSR programme.

Evidence and exact acceptance criteria are in
[`outputs/STRATEGY_v2.md`](outputs/STRATEGY_v2.md). Update every entry at each
phase boundary; no item may silently disappear.

## Binding constraints, limitations, and decisions

- [x] A1-MC: Raw MC NanoAOD coverage and a matched VBF jet/FSR pilot —
  resolved Phase 1 by `outputs/RAW_OBJECT_FEASIBILITY_v1.md` and
  `raw_object_pilot_vbf_v1.json` (11 fit processes; 1,198/1,198 pilot match).
- [ ] A1-data: Matching raw-data NanoAOD/manifest for `fake_data_10fb.root` —
  blocked externally; no raw-data NanoAOD exists under the supplied production
  root. Re-run the one-to-one event-match pilot when provenance is supplied.
- [ ] A2: Reference-like FSR candidate reconstruction — implement after the
  raw-data route exists; validate muon-indexed and electron-associated photons,
  pairing, and five-angle changes before fitting.
- [ ] A3: Dataset-specific luminosity uncertainty — Phase 2 must locate an
  applicable calibration/manifest or label $\mu$ conditional on supplied
  10 fb$^{-1}$ without importing a 35.9 fb$^{-1}$ uncertainty.
- [x] L1: 10 fb$^{-1}$ mass-statistics benchmark — recorded as 0.38 GeV only
  as a planning benchmark; it is not a scale uncertainty or result.
- [ ] L2: Generator/shower/PDF/scale coverage — Phase 2 inventory must locate
  usable weights/alternates and record an evidence-backed disposition.
- [x] D1: Detector-level simultaneous $(m_H,\mu)$ profile likelihood —
  retained in `STRATEGY_v2.md`; no unfolding or closed-form extraction claim.
- [x] D2: $105<m_{4\ell}<140$ GeV in every fit and variation — retained.
- [x] D3: Independent flavour-separated $Z$ calibration — exact data/MC
  transformation, covariance, and closure contract added.
- [ ] D4: Angle-only NN primary with cut-based cross-check — implement both;
  apply the deterministic gate before assigning the primary role.
- [ ] D5: Exclusive raw-jet VBF/non-VBF category — retain the reference-like
  $\mathcal D_{2\mathrm{jet}}^{\rm VBF}>0.5$ requirement; no proxy allowed.
- [x] D6: DY+jets MC is nominal fake template — retained with transfer,
  covariance, composition, and closure contract.
- [x] D7: MC pseudo-data then deterministic 10% data — retained.
- [x] D8: Supplied cross sections and Metadata `nEvents` normalization —
  retained; no data-derived normalization.
- [x] D9: No particle-level/unfolded claim — retained.
- [x] D10: Separate absolute data/MC calibration to common $M_Z$ — exact
  transformation and common-$M_Z$ correlation added.
- [ ] D11: Reference-like constrained-mass pseudoexperiment cross-check —
  implement and compare primary versus reference-like fit before data result.

## Reconstruction, calibration, and VBF validation

- [ ] Raw-object producer for each nominal and systematic MC input — persist
  event IDs, cleaned jets, raw/JEC/JER fields, DeepJet inputs, FSR/generic
  photons, and lepton inputs.
- [ ] Raw-data pilot — demonstrate one-to-one event match and the same object
  branches before using any VBF/FSR category in fake data.
- [ ] VBF viability — publish final-state/process yield, efficiency, purity,
  composition, JEC/JER/b-tag migration, and profiled Fisher-information table.
- [ ] Angular preservation — report matched-event pre/post-FSR/calibration
  pairing agreement and all five angle shifts by final state.
- [ ] $Z$ controls — use frozen trigger/ID/isolation/vertex/mass selection and
  fit data and MC separately in each flavour.
- [ ] $Z$ nominal model — extended unbinned BW $\otimes$ double-Crystal-Ball
  plus first-order background on 80--100 GeV; persist covariance and GoF.
- [ ] $Z$ alternatives — Voigtian/exponential, 0.25 GeV binned, and 78--104
  GeV domain form one non-double-counted `ZFitModel` envelope.
- [ ] $Z$ subdivision gate — require $N\ge10^4$ per proposed bin, global-vs-
  split $p<0.05$, two-standard-deviation correction difference, and injected
  closure before adding $\eta$/$p_T$ bins.
- [ ] Injected calibration closure — run $\pm0.05\%$ and $\pm0.10\%$ e-only
  and $\mu$-only injections in all final states on independent MC halves;
  require the predeclared pull/coverage criteria.
- [ ] Scale-impact report — show before/after Z peaks, joint covariance, and
  profiled $\Delta m_H^{\rm scale}$ and $\Delta\mu^{\rm scale}$; interpret
  `z_scale_sanity_preflight_v1.json` only as a floor/alert record.

## Fake-background and NN validation

- [ ] Tight/loose fake definitions — apply the frozen $L$, $P$, $F$, `4P`,
  `3P1F`, `2P2F`, and `2P2LSS` definitions by flavour/category/score bin.
- [ ] DY transfer — derive $T_{r,f,c,q,b}$ and the $f/(1-f)$ estimator with
  prompt-$ZZ$ subtraction and 2P2F double-count correction.
- [ ] DY closure — independent SHA-256 split and held-out same-sign closure
  for rate, $m_{4\ell}$, and NN score; report all cells and toy coverage.
- [ ] DY nuisance implementation — rate plus normalized correlated mass/score
  shape modes from the transfer covariance; no flat catch-all error.
- [ ] Reducible-composition inventory — control/model `TTBar`, conversion/
  $Z\gamma^{(*)}$, and electroweak multileptons separately or quantify their
  propagated effect; never leave “other nonprompt”.
- [ ] NN training — deterministic 60/20/20 split, five seeds, class-balanced
  signal versus explicit ZZ/ggZZ/DY/TT mixture, with score templates for every
  fit background.
- [ ] NN input and score gates — satisfy input $\chi^2/\mathrm{ndf}<5$ or
  correction, held-out/sideband/fake-control $p\ge0.05$, sculpting limit, and
  seed/architecture stability before primary use.
- [ ] NN score-bin gate — retain only bins meeting the frozen effective-count
  condition; merge a sparse score axis, never the VBF category.
- [ ] Selection decision — run both approaches with identical assumptions;
  NN is primary only if every quality and expected-uncertainty gate passes.

## Fit preflight and method parity

- [ ] Expected-yield table — publish every final-state/VBF/score cell’s signal,
  background hierarchy, effective MC count, nuisance impact, and bin map.
- [ ] Template statistics — use Barlow--Beeston-lite/gamma constraints for
  finite weighted MC templates.
- [ ] Mass binning and morphing — test 0.5/1/2 GeV bins and leave-one-node
  morph closure on the 120--130 GeV, 0.1 GeV mass grid.
- [ ] Independent fit toys — run 500 event-disjoint toys for $(124,1)$,
  $(125,0.8)$, $(125,1.2)$, and $(126,1)$ with pull, coverage, boundary, and
  toy-GoF acceptance tables.
- [ ] VBF information — demonstrate finite non-prior-dominated profiled VBF
  information; formal revision is required if it fails.
- [ ] Reference-like mass cross-check — implement $m'_{4\ell}$ constrained fit,
  per-event resolution proxy, and kinematic-score/MELA option; compare both
  POIs and coverage to the primary fit within the predeclared tolerance.
- [ ] HIG-16-041 parity table — retain source-by-source likelihood,
  constrained-mass, category, resolution, fake, and signal-shape comparison.

## Systematic sources

- [ ] $Z\to ee$ peak statistical covariance — joint flavour/category
  propagation.
- [ ] $Z\to\mu\mu$ peak statistical covariance — joint flavour/category
  propagation.
- [ ] $Z$ model/binning/domain/residual — one discrete/envelope group.
- [ ] Common external $M_Z$ — one fully correlated rank-one component.
- [ ] Lepton resolution — flavour-correlated line-shape nuisance.
- [ ] Lepton reconstruction/ID/isolation/trigger — data/MC scale factors.
- [ ] Luminosity — exact-dataset calibration or conditional-$\mu$ treatment.
- [ ] Pileup and JEC/JER — raw-object category-migration variations.
- [ ] b tagging/mistag — raw-object VBF/top-control migration variations.
- [ ] Signal acceptance/production mixture/PDF/QCD scale — evidence-backed
  in-file weights or alternate configuration.
- [ ] Shower/hadronization/underlying event — independent configuration or
  documented, tested limitation.
- [ ] MC statistics and mass morphing — gamma constraints and closure.
- [ ] $q\bar q$/$gg$ ZZ theory and sideband variations — separate processes.
- [ ] DY fakes — control statistics, composition, and transfer covariance.
- [ ] Top/conversion/electroweak reducibles — distinct template/control effect.
- [ ] NN modelling — calibration, replica, binning, and decorrelation shifts.
- [ ] Fit discretization/model — toy-backed predeclared alternatives.

## Comparisons and flagship figures

- [ ] $m_H$ comparison — HIG-16-041 Table 6 and PDG target, with overlap/
  covariance decision before any pull.
- [ ] $\mu$ comparison — HIG-16-041 Table 3 under matched mass convention;
  explicitly state why a world-$\mu$ comparison is not meaningful.
- [ ] $m_{4\ell}$ comparison — common selection/binning and covariance-aware
  statistic only if provenance permits it.
- [ ] VBF/jet comparison — HIG-16-041/HEPData only after raw-data route.
- [ ] NN comparison — state non-equivalence to published MELA discriminant;
  use score controls and reference-like fit instead of an invalid pull.
- [ ] Figure 1: before/after $ee$ and $\mu\mu$ peaks with residuals.
- [ ] Figure 2: calibrated $m_{4\ell}$ by final state and fit window.
- [ ] Figure 3: VBF jet/discriminant validation and VBF-category mass.
- [ ] Figure 4: NN score in fit, sideband, and fake-control regions.
- [ ] Figure 5: simultaneous $(m_H,\mu)$ profiles and scans.
- [ ] Figure 6: uncertainty breakdown, scale contribution, pulls, covariance.
- [ ] Diagram: calibration-to-fit chain.
- [ ] Diagram: rest-frame angle definition and NN inputs.
- [ ] Diagram: exclusive VBF/non-VBF/control-region topology.
