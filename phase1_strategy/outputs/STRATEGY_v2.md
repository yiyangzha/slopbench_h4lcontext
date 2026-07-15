# Phase 1 strategy v2: calibrated CMS Open Data $H\to ZZ^*\to4\ell$ measurement

## Status and scope

This version supersedes the Phase 1 plan in `STRATEGY.md` for downstream
work. It retains the requested detector-level simultaneous measurement of
$(m_H,\mu)$ in $4e$, $4\mu$, and $2e2\mu$, the CMS HIG-16-041 mass-fit window
$105<m_{4\ell}<140$ GeV, the four-lepton-rest-frame angular NN, a real
reconstructed-jet VBF category, flavour-separated $Z$-peak calibration, and
DY+jets MC as the nominal reducible-background template. It does **not** claim
an observed calibration, $m_H$, or $\mu$ result.

This repair adds the contracts requested by the first review. It does not
overwrite the original strategy or the committed pilot evidence. In particular,
the raw-object pilot proves an MC re-ntupling route, while a matching raw-data
NanoAOD route for the supplied fake-data product remains unavailable. Therefore
this version does **not** declare the final data VBF/FSR programme feasible or
replace it with a proxy. Phase 1 remains blocked from advancing until that
external provenance is supplied and independently rechecked.

The binding status of every item is tracked in
[`../COMMITMENTS.md`](../COMMITMENTS.md). The review-finding disposition is
recorded in [`FIX_RESOLUTION_v1.md`](FIX_RESOLUTION_v1.md).

## Retained binding decisions and new operational decisions

- **[D1] Detector-level simultaneous fit.** The primary endpoint remains a
  simultaneous profile-likelihood extraction of $(m_H,\mu)$ with shared,
  correlated nuisances. It is not an unfolded or fiducial measurement and not
  a closed-form counting extraction.
- **[D2] Reference mass window.** Every nominal and systematic fit uses
  $105<m_{4\ell}<140$ GeV. Wider ranges are only sideband or validation plots.
  The same window is used in CMS HIG-16-041 for its on-shell mass/width work
  [R1, Sec. 10.4].
- **[D3] Independent $Z$ calibration.** The $ee$ and $\mu\mu$ peaks are fit
  independently and anchored to the external PDG $M_Z$, never to a Higgs mass,
  signal strength, or published Higgs result [R4].
- **[D4] Angular NN plus cut-based cross-check.** The angle-only NN is planned
  primary inference information, while a transparent mass-only cut-based fit
  is a fully executed independent cross-check. The deterministic decision
  gates are in [Angular-NN decision contract](#angular-nn-decision-contract).
- **[D5] Reconstructed-jet VBF category.** An exclusive VBF/non-VBF split is
  fitted simultaneously after raw-object re-ntupling. $p_{\mathrm T}^{4\ell}$,
  candidate multiplicity, or any flat-tree proxy is forbidden as a replacement.
- **[D6] DY+jets nominal fake model.** `DYJetsToLL` remains the nominal
  reducible/fake shape. Data control regions calibrate its transfer and
  constrain its nuisance parameters; they do not replace it with a separate
  nominal fake model.
- **[D7] Staged data.** Phase 4a uses MC-only pseudo-data. Phase 4b uses the
  documented deterministic tenth of data. The user has authorized direct
  passage through the Phase 4b human gate only after its review and validation
  conditions are genuinely met.
- **[D8] Non-circular normalization.** MC weights are
  $w_s=\sigma_sL/\sum_i n_{\mathrm{Events},i}$, using the supplied effective
  cross sections, $L=10\,000\ \mathrm{pb}^{-1}$, and Metadata `nEvents` sums.
  Neither luminosity nor lepton scale is inferred by forcing a Higgs result to
  agree with a reference.
- **[D9] No particle-level claim.** No particle-level definition or response
  correction exists, so neither “unfolded” nor “fiducial” may describe this
  result. The candidate conventions are analogous validation sources only.
- **[D10] Absolute data-and-MC calibration convention.** Data and every MC
  template are separately transformed to the common external $M_Z$ reference;
  the full joint calibration covariance, including the common $M_Z$ component,
  is propagated once through both data and templates.
- **[D11] Reference-like mass-fit cross-check.** The binned primary fit is
  tested against a predeclared reference-like constrained-mass fit on
  statistically independent pseudo-data; it is not asserted to be equivalent
  without that test.

## Evidence already obtained and what it does not prove

### Raw objects and VBF/FSR

The committed record [`RAW_OBJECT_FEASIBILITY_v1.md`](RAW_OBJECT_FEASIBILITY_v1.md)
and its machine-readable pilot establish the following MC facts:

| Evidence | Result | Meaning |
|---|---:|---|
| Fit-MC raw coverage | 11/11 supplied processes have an `Events` tree and the required identity, jet, b-tag, FSR-photon, and generic-photon branch groups | MC re-ntupling is feasible for every currently named fit process. |
| VBF pilot event match | 1,198 / 1,198 flat candidates exactly match raw `(run,lumi,event)` | The pilot is a real MC re-ntupling, not a schema-only assertion. |
| In-window pilot candidates | 1,134 in $105<m_{4\ell}<140$ GeV; 671 have at least two cleaned raw jets | This is an object-preservation check, **not** VBF purity, efficiency, or a final VBF yield. |
| FSR inputs | 37 matched candidates have `FsrPhoton`, all with a muon-indexed association | Muon FSR inputs are preserved; electron-photon association still needs its dedicated validation. |
| Supplied fake data | Flat `h4lTree`/`dilepTree` products retain event identifiers but have no `Jet_*`, `Photon_*`, or `FsrPhoton_*` branches; the production root exposes no raw-data NanoAOD directory | Event identifiers cannot reconstruct missing objects. The data VBF/FSR route remains unresolved. |

The jet object baseline is fixed before signal-region fitting: cleaned jets
must have $p_{\mathrm T}>30$ GeV, $|\eta|<4.7$, `Jet_jetId>=2`, and
$\Delta R(\mathrm{jet},\ell)>0.4$. The final exclusive VBF tag is frozen to
the HIG-16-041 two-jet form: exactly four selected leptons; either two or
three selected jets with at most one medium-DeepJet b tag, or at least four
selected jets with no medium-DeepJet b tag; and
$\mathcal D_{2\mathrm{jet}}^{\rm VBF}>0.5$ [R1, Sec. 6]. The 2017
medium-DeepJet working point is taken from the applicable official calibration
payload during the Phase 2 raw-input inventory and recorded with that payload
version. The pilot baseline is not that final tag.

Before any VBF data fit, the raw-data manifest must provide an input that
matches `fake_data_10fb.root` one-to-one on `(run,lumi,event)`. The same
producer must preserve selected jets, raw/JEC/JER-relevant fields, DeepJet
inputs, `FsrPhoton`, generic photons, and all lepton inputs. On matched MC it
must then report, per process and final state, VBF efficiency, purity,
background composition, JEC/JER and b-tag migration, and profiled VBF Fisher
information. A category with zero or prior-dominated profiled information
requires a formal strategy revision; it may not be silently omitted or
relabelled as a $p_{\mathrm T}^{4\ell}$ category.

FSR recovery is applied before candidate pairing and calibrated masses:
muon photons use `FsrPhoton_muonIdx`; electron recovery is tested with retained
photons and a documented geometric/electron association. The matched-MC
validation must compare pre/post-recovery pairing and all five NN angles by
final state. The same comparison is repeated on data only after raw-data
provenance is proven.

### Current-flat yield inventory

[`mc_preselection_inventory_v1.json`](mc_preselection_inventory_v1.json)
uses the supplied normalization and gives a useful, deliberately limited
10 fb$^{-1}$ hierarchy in the fit window. It has no raw jets, FSR recovery,
VBF split, or NN score; it must never be presented as the final category table.

| Final state | Signal total | $q\bar q\to ZZ$ | $gg\to ZZ$ | DY+jets | $t\bar t$ |
|---|---:|---:|---:|---:|---:|
| $4\mu$ | 3.07 | 12.01 | 0.36 | 1.57 | 0.19 |
| $4e$ | 1.89 | 5.71 | 0.21 | 29.38 | 0.90 |
| $2e2\mu$ | 4.27 | 16.37 | 0.47 | 33.05 | 1.54 |
| All states | 9.24 | 34.09 | 1.04 | 64.00 | 2.63 |

The large DY estimate and the nonzero top component make a predictive
transfer-and-closure contract essential. The entries are current-flat
preselection counts only; all future category yields are recomputed after the
raw-object and calibration chain.

## Exact DY+jets fake transfer and closure contract

The control contract preserves DY+jets MC as nominal while preventing its rare
four-lepton tail from being used as an untested shape.

### Frozen object and region definitions

The raw-object producer first applies the reference four-lepton preselection:
trigger match, two OSSF pairs, $40<m_{Z_1}<120$ GeV,
$12<m_{Z_2}<120$ GeV, every OS pair above 4 GeV, $\Delta R(\ell_i,\ell_j)>0.02$,
leading/subleading $Z_1$ lepton $p_{\mathrm T}\ge20/10$ GeV, and
$m_{4\ell}>70$ GeV [R1, Sec. 4]. It defines the following lepton states from
the retained NanoAOD branches before looking at data in the mass-fit window:

| State | Muon definition | Electron definition |
|---|---|---|
| Loose ($L$) | $p_T\ge5$ GeV, $|\eta|\le2.4$, $|d_{xy}|\le0.5$ cm, $|d_z|\le1$ cm, $\mathrm{SIP}_{3D}\le4$, `pfRelIso03<=0.35`, and PF-or-global | $p_T\ge7$ GeV, $|\eta|\le2.5$, the same impact-parameter/SIP/isolation limits |
| Tight ($P$) | Loose plus `mediumId` and `pfRelIso03<0.15` | Loose plus MVA-WP90 and `pfRelIso03<0.15` |

`4P` is the signal selection. `2P2F` contains a tight OSSF $Z_1$ and two
additional loose-but-not-tight leptons (`F`); `3P1F` contains a tight $Z_1$,
one additional tight lepton, and one additional loose-but-not-tight lepton.
They are mutually exclusive with `4P`. `2P2LSS` contains a tight $Z_1$ and
two additional loose same-sign, same-flavour leptons. It is disjoint from the
OS signal selection and is an independent validation region. Each definition
is retained separately by final state, VBF/non-VBF category, and planned score
bin. If a raw-data route is not available, these definitions are exercised on
MC but cannot validate a fake-data fit.

### Transfer into the fitted DY template

The nominal DY $m_{4\ell}$ and $D_{\rm NN}$ shapes come from weighted
`DYJetsToLL` MC after the same category and score assignment as `4P`. For a
flavour/category/score cell $(f,c,q)$, define a MC transfer

$$
T_{r,f,c,q,b}=\frac{N^{\rm DY,4P}_{f,c,q,b}}
{N^{\rm DY,r}_{f,c,q,b}},\qquad r\in\{\mathrm{3P1F},\mathrm{2P2F},\mathrm{2P2LSS}\},
$$

with a coarser, predeclared binning only when the denominator is empty or has
negligible effective statistics. The bin map and every merge are frozen from
MC before any fit-window data is used. The control-derived misidentification
probability is $f_\ell(p_T,|\eta|,\mathrm{flavour})$, measured in a
Z-plus-loose-lepton data control sample after prompt contamination subtraction.
The transfer estimator uses the CMS 3P1F/2P2F weights $f/(1-f)$, includes the
prompt-$ZZ$ subtraction in 3P1F, and subtracts the 2P2F double-counting term
as specified in HIG-16-041 [R1, Sec. 7.2.1]. It calibrates the DY control-to-
signal normalization; it does not replace the DY MC nominal template.

The workspace has one DY rate nuisance per correlated flavour/control source
and normalized shape modes derived from the full covariance of
$T_{r,f,c,q,b}$ and the measured fake probabilities. Rate modes preserve the
template shape; shape modes preserve its total rate. The same nuisance vector
acts jointly on the $m_{4\ell}$ and $D_{\rm NN}$ dimensions, so a score-shape
failure cannot be hidden by a mass-only normalization constraint. An arbitrary
flat fake uncertainty is prohibited.

### Background-composition inventory

`TTBar` is a distinct template with a b-enriched control region (at least one
selected b-tagged jet) and is not merged with DY absent a demonstrated common
transfer. Conversion/$Z\gamma^{(*)}$ and electroweak multilepton components
are separate required inventory rows: Phase 2 must locate a sample/control
constraint or quantify their propagated effect in every affected cell. “Other
nonprompt” is not a valid disposition. An unavailable component cannot be set
to zero; it triggers an evidence-backed strategy update before the fit.

### Independent closure, coverage, and failure rule

The DY MC is split by `SHA-256(run:lumi:event:DYFAKEv1) mod 2` before any
transfer is derived. One half derives $f_\ell$-equivalent transfer corrections
and covariance; the other is the closure target, then the roles are exchanged.
The 2P2LSS region is held out as an independent shape validation. For every
populated $(f,c,q)$ cell, validate both the total rate and the binned
$m_{4\ell}$ and $D_{\rm NN}$ shapes. Use a covariance-aware $\chi^2$ when
asymptotic bins are adequate and a saturated-likelihood toy $p$ value when
they are not.

The predeclared pass conditions are a rate pull below 2, shape $p\ge0.05$ in
each tested observable, and nominal 68% interval coverage between 0.62 and
0.74 in 500 independent Poisson/control-statistics toys. Report the full
cell-by-cell table, including cells that fail. A failure requires a traced
composition or transfer-model remedy and a repeated independent closure. It
cannot be absorbed into an unconstrained normalization, a broad catch-all
nuisance, or a silent replacement of DY MC.

## $Z$-calibration convention, closure, and scale reach

### Frozen control selection and peak model

The numerical preflight in
[`z_scale_sanity_preflight_v1.json`](z_scale_sanity_preflight_v1.json) fixes a
reproducible starting control selection: supplied OSSF candidate,
$70<m_{\ell\ell}<110$ GeV, $n_{\rm PV}>0$,
$|d_{z,1}-d_{z,2}|<0.1$ cm, and both `relIso03<0.15`. The $\mu\mu$ control
requires a recorded single/double/triple-muon trigger bit, two medium-ID
PF-or-global muons; the $ee$ control requires a recorded electron trigger bit
and two MVA-WP90 electrons. These exact control cuts are independent of
$m_H$, VBF category, NN score, and the Higgs reference mass.

For each flavour and for data and MC separately, the nominal calibration fit
is an extended unbinned likelihood on $80<m_{\ell\ell}<100$ GeV with a
relativistic Breit--Wigner convolved with a double-sided Crystal-Ball response
and a first-order positive background component. The alternative-model set is
predeclared as (i) Voigtian plus exponential background, (ii) a binned fit
with 0.25 GeV bins, and (iii) the $78<m_{\ell\ell}<104$ GeV domain. Each fit
must converge, have a finite covariance, and pass a binned goodness-of-fit
$p\ge0.05$; failed alternatives are diagnostic failures rather than sources
to be summed automatically.

Let $\widehat m^{X}_{Z,f}$ be the fitted location for flavour $f$ and sample
$X\in\{\mathrm{data},\mathrm{MC}\}$. The exact nominal transformation is

$$
c_f^X=\frac{M_Z^{\rm ref}}{\widehat m^{X}_{Z,f}},\qquad
p^{\mu,X}_{\ell,\rm cal}=c_f^X p^{\mu,X}_{\ell,\rm reco}.
$$

Thus data and every MC template are each placed on the same external scale;
the MC peak is retained as a measured response input, not discarded and not
used to replace $M_Z^{\rm ref}$. After this transformation, rebuild FSR,
pairing, $m_{4\ell}$, rest-frame angles, NN score, VBF variables, and all fit
templates. For $2e2\mu$, use $c_e^X$ for electrons and $c_\mu^X$ for muons.
A post-hoc fractional shift of only $m_{4\ell}$ is forbidden.

The calibration parameter vector is

$$
\mathbf q=(c_e^{\rm data},c_\mu^{\rm data},c_e^{\rm MC},c_\mu^{\rm MC},M_Z^{\rm ref}).
$$

All nominal and varied rebuilds draw/scan this **joint** covariance. The PDG
$M_Z$ component is a single rank-one, fully correlated electron/muon component;
it is propagated once rather than accidentally counted once per flavour or
cancelling without being documented. Electron components correlate $4e$ with
$2e2\mu$, muon components correlate $4\mu$ with $2e2\mu$, and each spans VBF
and NN cells.

Peak-fit statistical covariance is one nuisance group. The correlated
model/binning/domain/residual alternatives above form one discrete
`ZFitModel` envelope after removing the common statistical covariance; they
are not independent errors summed in quadrature. Supported $\eta$ or $p_T$
subdivisions are permitted only if every proposed bin has at least $10^4$
selected candidates, the split versus global likelihood-ratio test has
$p<0.05$, the correction difference is significant at two standard deviations
using the joint covariance, and independent injected-scale closure succeeds.

### Injected-scale closure and numerical sanity interpretation

Use independent event-hash MC halves for calibration derivation and closure.
Inject $\pm0.05\%$ and $\pm0.10\%$ electron-only and muon-only scale offsets,
then run the full calibration, reconstruction, template, and simultaneous-fit
chain in $4e$, $4\mu$, and $2e2\mu$. A valid configuration recovers each
injection with absolute mean pull below 0.2, pull width in [0.8, 1.2], and
68% coverage in [0.62, 0.74] over 500 independent toys. This test is required
before a kinematic subdivision or a calibrated data fit is accepted.

The committed preflight is intentionally modest. It selected 4,508,686 $ee$
and 6,284,267 $\mu\mu$ broad-window candidates. Its RMS/$\sqrt N$ proxy gives
relative statistical floors of $2.57\times10^{-5}$ ($e$) and
$2.07\times10^{-5}$ ($\mu$), corresponding to channel mass floors of
3.22 MeV ($4e$), 2.59 MeV ($4\mu$), and 2.07 MeV ($2e2\mu$). These are not
fitted peak positions, calibration corrections, or total scale uncertainties;
the broad-window means in that JSON must not be interpreted as either.

The same JSON records 0.152 GeV for $m_H$ and 0.360 for absolute $\mu$ as
**reference-scaled investigation thresholds**, derived from HIG-16-041 total
uncertainties at 35.9 fb$^{-1}$, not assigned nuisances and not a scale-error
prediction. The final scale-only impact is instead computed from the
profiled fit with and without the joint scale nuisance group,

$$
\Delta x_{\rm scale}=\sqrt{\sigma_x^2(\mathrm{all})-
\sigma_x^2(\mathrm{scale\ fixed})},\quad x\in\{m_H,\mu\}.
$$

An impact below the relevant control-statistics floor, an $m_H$ impact above
the 0.152 GeV alert threshold, or a $|\Delta\mu_{\rm scale}|$ above 0.360
requires a documented calibration/model investigation. It is not a licence to
assign either threshold as an uncertainty. The Phase 4 result must show the
before/after $Z$ peaks, calibration covariance, fitted $(m_H,\mu)$, and this
scale-only contribution.

## Angular-NN decision contract

The NN inputs are exactly

$$
\{\cos\theta^*,\cos\theta_1,\cos\theta_2,\Phi,\Phi_1\},
$$

calculated from calibrated, FSR-recovered four-lepton candidates in their
defined rest frames. The NN never receives $m_{4\ell}$, lepton momentum scale,
event identifiers, event weights, category labels, or a proxy for them.
Candidate pairing/angle conventions are first tested on matched MC before and
after raw-object reconstruction as required above.

Training uses a deterministic SHA-256 event split (60% train, 20% validation,
20% held-out test) and five recorded random seeds. The signal class combines
all supplied Higgs modes with their supplied 10 fb$^{-1}$ weights. The
background class explicitly includes $q\bar q\to ZZ$, $gg\to ZZ$, DY+jets,
and $t\bar t$ with their expected mixture, then is class-balanced only at the
class-total level. Every fit background is nevertheless evaluated through the
trained score and receives its own template/statistical treatment; no fit
background is implicitly absent because it was subdominant in training.

The NN may be primary only if all of the following are documented before the
data fit:

1. Every input has sideband and fake-control data/MC $\chi^2/\mathrm{ndf}<5$
   or a validated correction; score-level sideband and fake-control tests have
   covariance-aware or toy $p\ge0.05$.
2. Held-out versus training score distributions pass a two-sample $p\ge0.05$;
   seed/architecture replicas change each expected uncertainty by less than
   10% and do not change the selected score-bin map.
3. In background MC and controls, $|\rho_{\rm Spearman}(D_{\rm NN},m_{4\ell})|<0.10$
   and score-sliced mass-shape tests pass $p\ge0.05$. A failed sculpting test
   requires decorrelation/retraining, not a score threshold tuned on the
   signal window.
4. The score map is fixed from the training/validation MC: start with two
   equal-signal-quantile non-VBF bins and retain it only if each bin has
   $N_{\rm eff}\ge3$ signal and $N_{\rm eff}\ge10$ total background before
   mass binning. The VBF mass category remains present; it receives a score
   axis only if the same condition holds after its own category selection.
   Otherwise that score axis is merged, not the VBF category.
5. With identical calibration, fake, and fit assumptions, the expected
   NN-assisted fit has no worse expected $\sigma(m_H)$ **and**
   $\sigma(\mu)$ than the cut-based fit, passes its toy GoF/coverage gates,
   and has no model-quality failure.

The decision is deterministic: all five gates pass $\Rightarrow$ Approach B
is primary and Approach A is reported as cross-check; any gate fails
$\Rightarrow$ no silent cut-based substitution is allowed. The failure,
diagnostic, attempted remediation, and a formally reviewed strategy decision
are required before the analysis can change the NN's binding role.

## Simultaneous-fit preflight and toy plan

The base workspace contains the three final states crossed with exclusive
VBF/non-VBF categories. It shares $(m_H,\mu)$, the joint lepton-scale
covariance, fake transfer nuisances, process uncertainties, and MC statistical
constraints. The VBF split remains even if it is sparse. NN score bins are
added only through the frozen gate above.

The primary mass templates use 1 GeV bins in the fixed 105--140 GeV window;
the predeclared 0.5 and 2 GeV alternatives, and a rule prohibiting bins finer
than half the calibrated final-state mass resolution, test discretization.
Every finite-MC process/bin receives a Barlow--Beeston-lite (gamma) constraint
based on its effective weighted count. Empty components are explicitly zero
only after the process/category inventory has established that they have no
support; they are never silently folded into another background.

Mass templates are generated on a 120--130 GeV grid with 0.1 GeV spacing and
morphed only after leave-one-grid-node closure. The interpolation residual in
each populated bin must be within its MC statistical uncertainty and the
injected-mass bias must be below $0.1\,\sigma_{m_H}^{\rm stat}$. The template
grid, bin map, effective counts, expected per-process yields, and every merge
are persisted before data fitting.

The preflight has three independent layers:

| Layer | Required test | Acceptance rule |
|---|---|---|
| Occupancy and identifiability | Asimov table for every final-state/category/score cell, including expected signal/background yields, dominant process, $N_{\rm eff}$, and nuisance impact | All retained cells have a defined template-statistics treatment; profiled Hessian is finite; the VBF category contributes finite non-prior-dominated information. |
| Independent recovery | 500 toys generated from event-hash-disjoint MC at $(m_H,\mu)=(124,1),(125,0.8),(125,1.2),(126,1)$ | For both POIs, mean pull $|\bar p|<0.2$, width [0.8,1.2], 68% coverage [0.62,0.74], and no more than 5% of $\mu\ge0$ toys at the boundary. |
| GoF and resolving power | Saturated/deviance GoF from the toy ensemble and fixed-versus-profiled nuisance scans | Observed/validation deviance lies in the central 95% toy interval; the expected $\mu=0.8$ and $\mu=1.2$ injections are distinguishable from $\mu=1$ according to their reported profile intervals, or the limited resolving power is explicitly stated. |

The fit reports joint and one-dimensional profile scans, nuisance pulls and
constraints, template-statistics impacts, category Fisher information, and
the with/without-scale uncertainty comparison. A fit with zero $\chi^2$,
identical toy outputs, wholesale score-bin exclusion, or a profile determined
by priors rather than events triggers an investigation rather than acceptance.

## HIG-16-041 method-parity and cross-check matrix

The primary fit is deliberately not called identical to the CMS measurement.
The matrix below makes the gap testable and commits a reference-like
cross-check instead of asserting parity.

| Published HIG-16-041 feature | Source provenance | v2 primary method | Binding parity action |
|---|---|---|---|
| 1D, 2D, and 3D mass likelihoods using $m'_{4\ell}$, per-event mass information, and $\mathcal D^{\rm kin}_{\rm bkg}$ | [R1, Sec. 10.3, Table 6] | Binned calibrated $m_{4\ell}$ with optional angular-NN score | Build an independent constrained-mass 1D/2D/3D-like pseudoexperiment fit: $m'_{4\ell}$; then $m'_{4\ell}$ plus per-event calibrated mass-resolution proxy; then add the declared kinematic score. If a cited maintained MELA implementation is available, repeat with $\mathcal D^{\rm kin}_{\rm bkg}$; otherwise label the score replacement as a non-equivalent cross-check. |
| $m(Z_1)$ constrained refit | [R1, Sec. 10.3] | Direct calibrated four-vector mass | Implement and validate the constrained refit after FSR/calibration; compare primary and constrained-fit bias/coverage on independent toys. |
| Simultaneous categories and shared nuisance profiling | [R1, Secs. 6, 10] | Three flavours × VBF/non-VBF, shared $(m_H,\mu,\theta)$ | Retain the required VBF category. The smaller 10 fb$^{-1}$ programme does not claim the reference's category granularity; any coarsening is documented by the occupancy preflight. |
| Per-event signal resolution and signal shapes | [R1, Sec. 8, Table 6] | Calibrated mass templates with a separate resolution nuisance | Persist the resolution proxy/response inputs and test mass-morph closure. Missing resolution inputs block the reference-like 2D/3D cross-check rather than being silently assumed. |
| Data-derived reducible Z+X estimate | [R1, Sec. 7.2] | User-required DY MC nominal shape | Keep DY nominal, but enforce the explicit 3P1F/2P2F/SS transfer, composition, and independent closure contract above. |

For each independent pseudoexperiment configuration, report
$(m_H^{\rm primary}-m_H^{\rm reference-like})/\sigma_{\rm combined}$ and the
corresponding $\mu$ difference. Both absolute differences must be below 0.2
combined standard deviations, with valid coverage in both fits, before the
binned primary method can be treated as adequate for this 10 fb$^{-1}$
measurement. A failure is a method-parity finding, not an invitation to drop
the check.

## Systematic programme and Phase 2 feasibility deadlines

Every row below is a binding disposition, not an open inventory placeholder.
All variations are propagated through reconstruction, categories, templates,
and the simultaneous fit.

| Source | Disposition and correlation |
|---|---|
| $Z\to ee$/$Z\to\mu\mu$ statistical calibration covariance | Will implement as the joint $\mathbf q$ covariance, correlated by flavour across all categories. |
| $Z$ fit model, binning, domain, and residual dependence | Will implement as one non-overlapping `ZFitModel` envelope/discrete profile after common covariance removal. |
| External $M_Z$ | Will implement once as a fully correlated rank-one component. |
| Lepton resolution | Will implement as a distinct flavour-correlated signal-line-shape nuisance from the data/MC $Z$ width residual. |
| Lepton reconstruction, ID, isolation, trigger | Will implement from data/MC scale factors with flavour/kinematic correlations. |
| Luminosity | Phase 2 deadline: locate a calibration/manifest applicable to this exact supplied dataset. If absent, report $\mu$ explicitly conditional on the supplied 10 fb$^{-1}$ and do not invent a 2016 CMS uncertainty. |
| Pileup, JEC/JER, b tagging/mistag | Will implement from raw-object inputs as VBF/non-VBF and top-control migration nuisances after the raw-data route is proven. |
| Signal acceptance, production mixture, PDF, QCD scale | Phase 2 deadline: inventory in-file weights and alternate configurations for every signal process. Implement available variations; if none exist, document the search and issue a formally reviewed disposition rather than reporting zero. |
| Shower, hadronization, underlying event | Same Phase 2 inventory/deadline; propagate an independent available configuration/variation to category and NN shapes, or document a tested, evidence-backed limitation. |
| MC statistics and mass morphing | Will implement gamma constraints and leave-node morph closure. |
| $q\bar q$/$gg$ ZZ | Will implement cited theory/sideband rate and shape variations separately by process. |
| DY fakes | Will implement the transfer covariance, rate mode, and correlated mass/score shape modes from the contract above. |
| Top, conversion/$Z\gamma^{(*)}$, electroweak multileptons | Will inventory and either model/control each separately or quantify an explicitly propagated effect; no “other” catch-all. |
| NN modelling | Will implement input calibration, score-transfer, seed/architecture, score-bin, and decorrelation variations. |
| Fit discretization/model | Will implement the predeclared mass-bin and template alternatives with toy-backed envelope/discrete treatment. |

The Phase 2 deadlines occur before its artifact is closed, not at Phase 4:
raw-data provenance, luminosity applicability, and generator/weight coverage
must each have an evidence record. Failure of raw-data provenance leaves the
VBF/FSR data result blocked; it does not authorize a proxy or a reduction of
the user's VBF requirement.

## Reference comparison matrix

No comparison is called independent until the fake-data provenance establishes
whether it overlaps the 35.9 fb$^{-1}$ 2016 CMS dataset and provides a usable
covariance. “Difference” below means a descriptive signed difference with both
uncertainties reported, not a pull or p value.

| Result/distribution | Target and audit provenance | Comparability and overlap treatment | Statistic/reporting rule |
|---|---|---|---|
| $m_H$ | HIG-16-041 $125.26\pm0.20\,({\rm stat})\pm0.08\,({\rm syst})$ GeV, [R1, Table 6]; PDG Higgs listing [R5] | Same decay channel/energy and same fit window for HIG; fake-data overlap/covariance unresolved. PDG is a validation target, never a calibration input. | Use a covariance-aware pull only with documented covariance; otherwise report difference. A >3$\sigma$ or >30% deviation triggers the review investigation rule. |
| $\mu$ | HIG-16-041 $1.05^{+0.19}_{-0.17}$ at fixed $m_H=125.09$ GeV, [R1, Table 3] | Compare only after reproducing the fixed/conditional mass convention. No universal world $\mu$ is directly comparable because channel, production assumptions, theory normalization, and overlap differ. | HIG: covariance-aware pull only if overlap known; otherwise difference. World average: explicitly “not meaningfully comparable,” not an omitted target. |
| Calibrated $m_{4\ell}$ | HIG-16-041 mass spectra/selection, [R1, Sec. 4 and figures on the official result page] | Same channel/window after rebinning to a common selection. The calibration control is not independent of its external $M_Z$ input. | Shape comparison with covariance-aware $\chi^2$ only after common bins/covariance; otherwise overlay and documented qualitative selection difference. |
| VBF tag/jet observables | HIG-16-041 categorization and HEPData record 80189 [R1, Sec. 6; R1H] | Requires raw data and the reference-like jet definition; current flat tree is not comparable. | Compare binned distributions/efficiency only after data route and matching selections; otherwise state “not yet comparable: raw-data jet provenance absent.” |
| Angular NN score | No identical published observable; HIG uses matrix-element discriminants [R1, Sec. 10.3] | The five underlying angles are physically related but $D_{\rm NN}$ is not $\mathcal D^{\rm kin}_{\rm bkg}$. | No numerical reference pull. Report score controls and the reference-like-fit comparison as validation, with the non-equivalence stated. |
| $Z$ peaks and scale | PDG $M_Z$ [R4] | This is a calibration input, not an independent validation target. | Show pre/post fits and closure/injected-offset recovery; do not quote a pull to the input as a validation result. |

## Flagship figures, handoff, and sources

The binding Phase 5 flagship figures remain: (1) before/after $ee$ and
$\mu\mu$ peaks with residuals; (2) calibrated $m_{4\ell}$ by final state;
(3) VBF jet/discriminant validation and VBF mass category; (4) NN score in
fit, sideband, and fake controls; (5) simultaneous $(m_H,\mu)$ profiles; and
(6) uncertainty, pulls, and scale-covariance breakdown. Phase 5 also draws
the calibration-to-fit chain, five-angle definition, and exclusive
VBF/non-VBF/control topology. Each is tracked in `COMMITMENTS.md`.

Phase 2 must produce the raw-data, luminosity, generator, and composition
inventories; Phase 3 implements both selection approaches, raw-object
re-ntupling, FSR, calibration, fake controls, and independent closures; Phase
4 builds the workspace and all predeclared toy tests. A later phase may not
replace any decision here without a documented revision and review.

- **[S1] Supplied task specification:** `../../prompt.md` — inputs,
  cross sections, target luminosity, and Metadata normalization rule.
- **[S2] Raw-object feasibility:** `RAW_OBJECT_FEASIBILITY_v1.md` and
  `raw_object_pilot_vbf_v1.json` — MC branch coverage/event match and the
  unresolved fake-data raw-input absence.
- **[S3] Current-flat inventory:** `mc_preselection_inventory_v1.json` —
  supplied-normalization 10 fb$^{-1}$ preselection hierarchy only.
- **[S4] Numerical scale preflight:** `z_scale_sanity_preflight_v1.json` —
  control-statistics floor and explicitly non-nuisance alert thresholds.
- **[R1] CMS HIG-16-041:** [official CMS result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-16-041/), JHEP 11 (2017) 047, [open-access paper](https://cds.cern.ch/record/2272260/files/cms-hig-16-041-arxiv.pdf), especially Secs. 4, 6, 7.2, 8--10 and Tables 3 and 6.
- **[R1H] HIG-16-041 HEPData:** [record 80189](https://www.hepdata.net/record/ins1608162).
- **[R2] CMS HIG-13-002:** [official result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-13-002/index.html), Phys. Rev. D 89 (2014) 092007.
- **[R3] CMS HIG-21-019:** [official result page](https://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-21-019/), Phys. Rev. D 111 (2025) 092014.
- **[R4] PDG $Z$ listing:** [2025 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-z-boson.pdf).
- **[R5] PDG Higgs listing:** [2025 listing](https://pdg.lbl.gov/2025/listings/rpp2025-list-higgs-boson.pdf).
- **[M1] Local methodology:** `methodology/03-phases.md`,
  `methodology/06-review.md`, and `methodology/12-downscoping.md` — phase
  gates, toy/closure requirements, comparison investigations, and the rule
  against silent downscoping.
