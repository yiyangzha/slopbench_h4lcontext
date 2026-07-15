# Phase 1 strategy: calibrated CMS Open Data $H\to ZZ^*\to4\ell$ measurement

## Summary

This strategy defines a detector-level measurement of the Higgs-boson mass
$m_H$ and signal-strength modifier $\mu$ in the $4e$, $4\mu$, and
$2e2\mu$ final states. The primary result is a simultaneous
multi-category profile-likelihood fit to calibrated $m_{4\ell}$ in the
published $105<m_{4\ell}<140$ GeV fit window. It includes separate
$Z\to ee$ and $Z\to\mu\mu$ peak calibrations, a VBF category, and a
four-lepton-rest-frame angular neural-network (NN) discriminant. DY+jets
simulation is the primary reducible/fake-background template.

The scale correction and its residual uncertainty are propagated by
recomputing four-lepton candidates and fit templates. No observed $m_H$,
$\mu$, or calibration number is claimed in Phase 1. Numbers below are
supplied normalizations or external validation targets and are cited in
[Sources](#sources).

## Retrieval status

The required experiment-corpus MCP operations were not exposed in this
session; the failed capability check is in ../retrieval_log.md. The
documented fallback is official CMS HIG-16-041 material, its HEPData
record, official CMS reference analyses, and PDG listings. This is a
source limitation, not a justification for using uncited remembered
constants.

## Binding constraints, limitations, and decisions

- **[A1] Raw-object prerequisite for VBF and FSR.** The delivered flat
  H4l tree has four-lepton and lepton scalar branches but no jets, b tags,
  missing transverse momentum, or FSR photons. This was confirmed from
  the supplied ntuplizer and from the data ROOT schema. Before Phase 3,
  re-ntupling from the corresponding NanoAOD must retain selected jets,
  jet energy/resolution inputs, b tags, FSR photons, and event identifiers.
  It is invalid to infer a VBF tag from $p_\mathrm{T}^{4\ell}$ or claim
  FSR correction from the current flat tree. [S2]

- **[A2] Reference-like candidate reconstruction.** The current ntuplizer
  admits $m_{Z_1}$ above 12 GeV; the reference baseline requires
  $m_{Z_1}>40$ GeV. Apply the stricter requirement now where possible and
  make the replacement ntuplizer form candidates after FSR recovery and
  calibrated lepton four-vectors. [R1, Sec. 4]

- **[A3] Luminosity uncertainty is unresolved.** The supplied metadata and
  task specify 10 fb$^{-1}$ but no uncertainty on it. Phase 4 must retrieve
  a dataset-appropriate luminosity calibration or report $\mu$ conditional
  on the supplied luminosity. It must not copy the 2.5% uncertainty from
  the distinct 35.9 fb$^{-1}$ CMS data set. [S1; R1, Sec. 9]

- **[L1] Statistical reach.** The supplied exposure is smaller than the
  35.9 fb$^{-1}$ HIG-16-041 sample. Since independent-count statistical
  uncertainties scale as $1/\sqrt L$, its 0.20 GeV mass statistical
  uncertainty scales to $0.20\sqrt{35.9/10}=0.38$ GeV for 10 fb$^{-1}$.
  This is a planning sanity benchmark, not a result or substituted
  uncertainty. [S1; R1, Sec. 10.3]

- **[L2] Independent-generator coverage is unverified.** Input directories
  are labelled nominal. Phase 2 must inventory alternate generator,
  shower, PDF, and scale information before a generator-model uncertainty
  is quoted. Missing alternatives do not permit omitting the test.

- **[D1] Technique.** The primary endpoint is a detector-level simultaneous
  template/profile-likelihood measurement of $(m_H,\mu)$ with correlated
  nuisances. It is neither a closed-form double-tag extraction nor an
  unfolded particle-level spectrum. [R1, Secs. 8--10]

- **[D2] Fit window.** The nominal and every systematic $(m_H,\mu)$ fit
  uses $105<m_{4\ell}<140$ GeV, as required by the task and used by
  HIG-16-041. Wider ranges are validation or sideband plots only. [R1,
  Sec. 10.4]

- **[D3] Independent scale calibration.** Separate $ee$ and $\mu\mu$
  $Z$ peaks, not a Higgs mass reference or fitted $m_H$, determine lepton
  scale corrections. The external target is
  $M_Z=91.1876\pm0.0021$ GeV. [R4]

- **[D4] Angular NN is primary.** The NN is actual inference information,
  not a plot. A transparent cut-based result is the independent cross-check.
  If the NN fails its data/MC input tests, the inputs must be calibrated or
  this strategy formally revised; silently substituting cuts is prohibited.

- **[D5] VBF category is binding.** After [A1], a mutually exclusive
  reconstructed-jet VBF category is fitted simultaneously with non-VBF
  categories. The HIG-16-041 corresponding category had about 49% VBF
  purity at 35.9 fb$^{-1}$; that is context, not an assumed purity here.
  [R1, Sec. 6]

- **[D6] Fake model.** DY+jets MC is the required primary fake-background
  template. Fake-enriched data control regions constrain and validate it;
  they do not replace it with an unrelated fake model. [S1]

- **[D7] Staged data.** Phase 4a uses MC-only pseudo-data. Phase 4b uses a
  reproducible tenth of the data selected by a documented hash of
  run, luminosity block, and event number, without post-inspection retuning.
  The user has authorized direct passage of the Phase 4b human gate once
  those validation checks are complete. [S1]

- **[D8] Non-circular inputs.** MC weights use supplied effective cross
  sections, supplied luminosity, and summed Metadata.nEvents. Luminosity
  and lepton scales must never be derived by forcing the Higgs fit to a
  published rate or mass. [S1]

- **[D9] No particle-level claim.** The result is calibrated detector-level
  parameter extraction. It cannot be called unfolded or fiducial without a
  new particle-level definition and validated response correction.

## Samples and normalization

The supplied data are 10 fb$^{-1}$ at 13 TeV. For sample $s$, use

$$
w_s=\frac{\sigma_sL}{\sum_i n_{\mathrm{Events}}^{(i)}}.
$$

The relation follows by requiring
$\sum_i n_{\mathrm{Events}}^{(i)}w_s=\sigma_sL$. The denominator is the
required sum of Metadata.nEvents, never selected flat-tree entries.
Metadata genEventSumw remains a diagnostic, not a silently substituted
denominator. [S1]

| Role | Input sample(s) | Effective cross section [pb] | Planned use |
|---|---|---:|---|
| Data | fake_data_10fb.root | -- | Four-lepton candidates |
| $Z$ control data | fake_data_dilep_10fb.root | -- | Separate $ee$ and $\mu\mu$ calibration fits |
| Signal | GluGluToHToZZ | $6.024\times10^{-3}$ | Inclusive signal template and NN signal training |
| Signal | VBF_HToZZ | $4.8794\times10^{-4}$ | VBF template and VBF efficiency |
| Signal | ZHToZZ, WPHToZZ, WMHToZZ | $9.8394\times10^{-5}$, $1.072352\times10^{-4}$, $6.706\times10^{-5}$ | Associated-production and migration component |
| Irreducible | ZZTo4L | $1.325$ | $q\bar q\to ZZ^{(*)}\to4\ell$ template |
| Irreducible | GGZZ2E2Mu, GGZZ4Mu, GGZZ4E | $3.185\times10^{-3}$, $1.575\times10^{-3}$, $1.619\times10^{-3}$ | $gg\to ZZ^{(*)}\to4\ell$ templates |
| Reducible/fake | DYJetsToLL | $5.396\times10^3$ | Required fake-lepton template |
| Reducible | TTBar | $5.270\times10^1$ | Top control and residual nonprompt component |

All values in this table are supplied inputs. [S1] An inspected VBF MC
file has Metadata.nEvents, nGenEvents, genEventSumw, effectiveXsecPb, and
selected-event counts but no jet branches in the flat H4l tree. [S2]

## Observables and shared reconstruction

The fitted observables are calibrated $m_{4\ell}$, a common signal
strength $\mu$, a four-lepton-rest-frame angular NN score
$D_\mathrm{NN}$, and VBF jet observables after [A1].

The NN uses only

$$
\{\cos\theta^*,\cos\theta_1,\cos\theta_2,\Phi,\Phi_1\},
$$

formed after assigning $Z_1$ as the pair closest to the nominal $Z$ mass
and boosting selected lepton momenta to the relevant rest frames. It has
no $m_{4\ell}$, energy-scale, event-ID, weight, or category-label input.
CMS documents this five-angle convention for $H\to ZZ\to4\ell$. [R2,
Fig. 8]

The raw-object replacement must implement a shared reference-like
preselection: trigger-matched, identified and isolated leptons; two
opposite-sign same-flavour pairs; pair masses in
$12<m_{\ell\ell}<120$ GeV; $m_{Z_1}>40$ GeV; every opposite-sign pair
above 4 GeV; at least two leptons above 10 GeV and one above 20 GeV; and
$m_{4\ell}>70$ GeV. FSR is recovered before candidate masses are formed.
These are HIG-16-041 values, not data-tuned thresholds. [R1, Sec. 4]

## Two qualitative selection approaches

### A. Transparent cut-based cross-check

Approach A applies the shared preselection and a fixed VBF-enriched
jet-topology category. It fits flavour times VBF/non-VBF calibrated
$m_{4\ell}$ templates without an NN-score axis. It gives a reproducible
cutflow and a direct cross-check of the primary result, at the expected
cost of less signal/background angular separation.

### B. Angular-NN primary selection and fit

Approach B uses the same preselection and VBF tag, then trains a compact NN
on the five rest-frame angles. Signal MC is separated from irreducible
$ZZ$ and DY fake templates using stratified process/final-state samples,
fixed recorded seeds, and an event-ID train/test split. The primary
likelihood uses calibrated $m_{4\ell}$ and $D_\mathrm{NN}$ templates in
non-VBF categories; the VBF category has its own mass template and gets
the score dimension only if expected population supports it.

Before B is allowed into the fit, Phase 3 must show:

1. held-out train/test score agreement and separation;
2. data/MC agreement for every NN input in sidebands and fake controls;
3. NN-score versus $m_{4\ell}$ correlation and mass-sculpting checks; and
4. stability of fitted $(m_H,\mu)$ under predeclared score binning and
   architecture/seed replicas.

B is primary if these gates pass; A is always fully run. The decision table
reports expected $\sigma(m_H)$, expected $\sigma(\mu)$, goodness-of-fit,
category occupancy, and input data/MC agreement under identical
calibration/background assumptions. A smaller expected uncertainty cannot
override poor model agreement.

## $Z$-peak calibration and residual scale propagation

The supplied dilepton tree contains separate $ee$ and $\mu\mu$ candidates
with ID/isolation, trigger bits, and $m_{\ell\ell}$. The clean control
selection requires trigger matching, an opposite-sign same-flavour pair,
tight flavour-appropriate ID/isolation, primary-vertex compatibility, and
a predeclared broad $Z$-pole fit domain. It is independent of the Higgs
mass fit, categories, and reference Higgs mass.

Fit data and MC separately in each flavour using a Breit--Wigner
line shape convolved with a resolution model, following the CMS
$Z$-based calibration strategy. [R1, Secs. 5 and 9] Let
$\widehat m_{Z,f}^{\mathrm{data}}$ denote the fitted peak and
$M_Z^\mathrm{ref}$ the PDG reference. Under a common flavour scale,
every lepton four-vector changes as $p^\mu\mapsto c_fp^\mu$, hence
$m_{\ell\ell}\mapsto c_fm_{\ell\ell}$. Enforcing the reference peak gives

$$
c_f=\frac{M_Z^\mathrm{ref}}{\widehat m_{Z,f}^{\mathrm{data}}},
\qquad p^\mu_{\ell,\mathrm{cal}}=c_fp^\mu_{\ell,\mathrm{reco}}.
$$

The limiting check is
$c_f\widehat m_{Z,f}^{\mathrm{data}}=M_Z^\mathrm{ref}$. MC peak fits
validate data/MC response and resolution; they never replace the external
target or tune $m_H$. Kinematic subdivisions are introduced only if
predeclared occupancy and nested-fit tests support them. [R1, Sec. 9; R4]

Apply corrections at lepton-four-vector level before recomputing
$Z_{1,2}$, $m_{4\ell}$, rest-frame angles, NN scores, VBF categories, and
templates. For $2e2\mu$, apply $c_e$ and $c_\mu$ to their corresponding
leptons. A post-hoc fractional mass shift is prohibited.

Residual scale nuisances come from the peak-fit covariance, justified
line-shape/tail alternatives, binned versus unbinned or predeclared binning
variants, fit-domain/residual tests, supported kinematic dependence, and
data/MC closure. Each is rerun through the full chain as named electron or
muon scale variations. Electron components are correlated across $4e$ and
$2e2\mu$; muon components across $4\mu$ and $2e2\mu$; both correlations
span VBF and NN categories. Resolution is a distinct signal-line-shape
nuisance. Phase 4 reports the profiled lepton-scale contribution to both
$m_H$ and $\mu$, before/after $Z$ peaks, and the comparison to [L1].

## Backgrounds, categories, and simultaneous fit

| Component | Classification | Nominal model | Validation/control |
|---|---|---|---|
| $H\to ZZ^*\to4\ell$ | Signal | Supplied ggH, VBF, VH templates with common $\mu$ and $m_H$ morphing | Split-MC closure and category/flavour consistency |
| $q\bar q\to ZZ$ | Irreducible | ZZTo4L templates | Sidebands and theory variations |
| $gg\to ZZ$ | Irreducible | Dedicated final-state GGZZ templates | Sidebands and theory variations |
| DY+jets fakes | Reducible/instrumental | Required DYJetsToLL template | Loose/not-tight and same-sign fake controls |
| $t\bar t$ and other nonprompt | Reducible | TTBar/control template | b-enriched and loose-lepton checks |

The DY template remains the primary fake model under [D6]. Control regions
constrain its rate and test transfer of its shape, but do not conceal a
shape disagreement. A 3P1F/2P2F-style and same-sign cross-check is retained
because the CMS reference derives reducible-background checks in such
regions. [R1, Sec. 7.2]

Fit categories are the three flavours crossed with VBF-enriched and
non-VBF selections; NN-score bins are added where templates are populated.
All share $m_H$, $\mu$, scale nuisances, and applicable process nuisances.

For category $c$ and bin $b$, the planned extended likelihood is

$$
\mathcal L(m_H,\mu,\boldsymbol\theta)=
\prod_{c,b}\operatorname{Pois}\left(
n_{cb}\mid\mu s_{cb}(m_H,\boldsymbol\theta)+
\sum_k b_{k,cb}(\boldsymbol\theta)\right)
\prod_j\pi_j(\theta_j).
$$

Here $n_{cb}$ is the observed count, $s_{cb}$ signal, $b_{k,cb}$
background, and $\pi_j$ a calibration, experimental, theory, or
control-region constraint. It is the extended binned analogue of the
published simultaneous mass/discriminant fit. If data equal the nominal
model and all nuisances are zero, the Poisson factors are maximized at the
nominal parameters, which is the required limiting check. [R1, Sec. 10.1]

Templates at multiple $m_H$ hypotheses are validated and morphed. Scanning
$m_H$ profiles $\mu$ and all nuisances; scanning $\mu$ profiles $m_H$ and
all nuisances. Outputs include joint and one-dimensional profile scans,
goodness-of-fit, nuisance pulls/constraints, and the change in both
uncertainties when lepton-scale nuisances are fixed. A reproducible pyhf
workspace, all variations, and calibration covariance are persisted.

## Binding systematic program

All relevant sources propagate through calibration, candidate rebuilding,
categories, templates, and fit. A copied flat error is not a systematic.

| Source | Variation or constraint | Correlation / status |
|---|---|---|
| $Z\to ee$ peak statistics | Data fit covariance | $4e$ and $2e2\mu$, all categories; will implement |
| $Z\to\mu\mu$ peak statistics | Data fit covariance | $4\mu$ and $2e2\mu$, all categories; will implement |
| Scale model/binning/residuals | Predeclared alternate fits and supported subdivisions | Separate correlated electron/muon variations; will implement |
| Lepton resolution | $Z$ width data/MC residual | Flavour-correlated line shape; will implement |
| Lepton ID, reconstruction, isolation, trigger | Tag-and-probe scale factors | Flavour/kinematic map; will implement |
| Luminosity | Dataset-specific documented calibration | Correlated MC rate; blocked by [A3], never guessed |
| Jet scale/resolution, pileup | Official raw-object variations | VBF/non-VBF migration; will implement after [A1] |
| b tagging/mistag | Official scale factors | VBF/top-control migration; will implement after [A1] |
| Signal acceptance, mixture, PDF, QCD scale | Official weights or documented alternate prediction | Process/category correlated; inventory required by [L2] |
| Shower, hadronization, underlying event | Independent configuration or defensible variation | Category and NN shape; inventory required by [L2] |
| MC statistics and mass morphing | Per-bin modifiers and morphing closure | Shape/statistical; will implement |
| $q\bar q$/$gg$ ZZ | Theory and sideband variation | Process rate/shape; will implement |
| DY fakes | Control-statistics, composition, and transfer tests | Final-state correlated template variation; will implement |
| Other reducible backgrounds | Cross-section/control variation | Separate from DY unless transfer proves correlation; will implement |
| NN modelling | Input calibration, replicas, binning, decorrelation | Score shape, no data-tuned threshold; will implement |
| Fit discretization/model | Predeclared binning/template variants | Closure-backed envelope/discrete study; will implement |

This matches the reference program's lepton scale/resolution, efficiency,
luminosity, reducible background, theory scale/PDF, category migration, jet,
and b-tag treatment. [R1, Sec. 9]

## Candidate-convention audit

### Extraction convention

The extraction convention applies to closed-form tagged/untagged
hemisphere counting, not [D1]. It is not the governing technique, but every
required source is audited:

| Required source | Applicability | Binding disposition |
|---|---|---|
| Tag/selection efficiency | Applicable analogue | Lepton and VBF-tag efficiency variations |
| Efficiency correlation | No double-hemisphere formula | Not applicable; shared object efficiencies are correlated nuisances |
| MC efficiency model | Applicable analogue | Generator, shower, PDF, scale acceptance/migration tests |
| Non-signal contamination | Applicable | Explicit ZZ, DY, top templates |
| Background composition | Applicable | DY/top/control composition variations |
| Hadronization model | Applicable | Shower/hadronization variations in category and NN templates |
| Physics parameters | Applicable | Theory acceptance and production-mixture inputs |
| Flavour composition | No inclusive closed-form input | Not applicable; explicit final-state templates |
| Production fractions | Applicable | Correlated mode-mixture and category-migration variations |

Its pseudo-data, independent closure, operating-point stability, and
ten-percent diagnostic principles are adopted as analogous fit validations.

### Unfolding convention

The unfolding convention is not applicable to the [D9] detector-level fit.
Its requirements are nevertheless enumerated:

| Requirement/gate | Status |
|---|---|
| Particle-level definition | Not applicable; no unfolded/fiducial claim |
| Response correction and alternative correction | Not applicable; scale calibration is not unfolding |
| Unfolded covariance and PSD test | Not applicable; persist workspace, nuisance correlations, and scans instead |
| Closure, stress, prior/model dependence | Applicable analogue; split-MC and injected-signal/scale tests are mandatory |
| Data/MC input validation | Applicable analogue; mandatory for calibration, VBF, and every NN input |

Any later fiducial spectrum reopens this decision and must meet every
unfolding requirement, including particle definition, covariance,
independent correction cross-check, closure, and stress tests.

## Reference analyses and comparison targets

| Reference | Method parity and systematic program | Numerical target and comparability |
|---|---|---|
| CMS HIG-16-041, JHEP 11 (2017) 047 | Direct 13 TeV mass/rate fit with categories, discriminant, $Z$ calibration, and simultaneous likelihood. Match lepton scale/resolution, ID/efficiency, luminosity, fakes, theory, jet/b-tag, migration treatment. | $m_H=125.26\pm0.20$ (stat) $\pm0.08$ (syst) GeV; $\mu=1.05^{+0.19}_{-0.17}$ at 35.9 fb$^{-1}$. Directly comparable observable/energy, not equal statistical power or production. [R1] |
| CMS HIG-13-002, Phys. Rev. D 89 (2014) 092007 | Earlier $4\ell$ mass fit with per-event mass, angular observables, $Z$ calibration, and dijet information. Match calibration, lepton efficiency, reducible-background, resolution, and category tests. | $m_H=125.6\pm0.4$ (stat) $\pm0.2$ (syst) GeV; rate $0.93^{+0.26}_{-0.23}$ (stat) $^{+0.13}_{-0.09}$ (syst) times SM. Different energy/luminosity mixture makes rate contextual. [R2] |
| CMS HIG-21-019, Phys. Rev. D 111 (2025) 092014 | Modern 13 TeV mass/width fit with per-event resolution, final-state categories, and profile scans. Match or improve scale/resolution, beam-spot/line-shape, and likelihood treatment. | $m_H=125.04\pm0.11$ (stat) $\pm0.05$ (syst) GeV at 138 fb$^{-1}$. Same channel/energy but much larger data and advanced reconstruction; precision context, never calibration input. [R3] |

HIG-16-041 HEPData supplies the published fiducial/differential values and
correlations needed for later jet-distribution comparisons. [R1H] The PDG
2025 Higgs listing gives $m_H=125.20\pm0.11$ GeV. The final result reports
the difference and covariance-aware compatibility, while acknowledging a
world average is a validation target, not independent calibration data.
[R5]

## Flagship figures and diagrams

The six binding Phase 5 flagship figures are:

1. Before/after $m_{ee}$ and $m_{\mu\mu}$ peaks with fit residuals.
2. Calibrated $m_{4\ell}$ in all flavour channels and the fit window.
3. VBF-score/jet validation and VBF-category $m_{4\ell}$.
4. Rest-frame angular-NN score in fit, sideband, and fake-control regions.
5. Simultaneous $(m_H,\mu)$ profile contour and one-dimensional scans.
6. $m_H$/$\mu$ uncertainty breakdown with lepton-scale contribution,
   pulls, and scale correlation matrix.

Phase 5 also supplies diagrams of the calibration-to-fit chain, the five
rest-frame angles/NN inputs, and exclusive VBF/non-VBF/control-region
topology. Every flagship plot is generated from persisted artifacts in PDF
and PNG and receives a self-contained analysis-note caption.

## Handoff and success criteria

Phase 2 inventories trees, metadata, all MC variations, and raw-NanoAOD
locations; verifies event-ID and Metadata.nEvents preservation after [A1];
and quantifies data/MC agreement for every selection and NN input. Phase 3
implements both approaches, raw-object VBF/FSR, and the independent
$Z$ calibration. Phase 4 writes the workspace, all systematic templates,
calibration covariance, fit outputs, and comparison tables. No later phase
may replace a [D] commitment without a documented strategy revision and
review.

No Phase 1 analysis script or figure is created: this is the written
strategy and evidence handoff. The schema inspection used Pixi and uproot.
Future scripts must be Pixi tasks, and the final all task must reproduce
re-ntupling, calibration, selection, systematic variations, fit, and plots.

## Sources

- **[S1] Supplied task specification:** ../../prompt.md. Authority for data
  paths, 10 fb$^{-1}$, cross sections, and the Metadata.nEvents
  normalisation requirement.
- **[S2] Supplied ntuplizer and inspected metadata:**
  ../../ntuplizer/h4l_ntuplizer.cpp and the supplied H4l ROOT schema,
  inspected on 2026-07-15. Authority for current branch availability.
- **[R1] CMS HIG-16-041:** [official CMS result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-16-041/) and [open-access paper](https://cds.cern.ch/record/2272260/files/arXiv%3A1706.09936.pdf), JHEP 11 (2017) 047, [doi](https://doi.org/10.1007/JHEP11(2017)047).
- **[R1H] HIG-16-041 HEPData:** [record 80189](https://www.hepdata.net/record/ins1608162).
- **[R2] CMS HIG-13-002:** [official CMS result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-13-002/index.html), Phys. Rev. D 89 (2014) 092007, [doi](https://doi.org/10.1103/PhysRevD.89.092007).
- **[R3] CMS HIG-21-019:** [official CMS result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-21-019/), Phys. Rev. D 111 (2025) 092014, [doi](https://doi.org/10.1103/PhysRevD.111.092014).
- **[R4] PDG $Z$ listing:** [2025 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf).
- **[R5] PDG Higgs listing:** [2025 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-higgs-boson.pdf).
