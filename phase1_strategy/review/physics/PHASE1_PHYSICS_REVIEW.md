# Phase 1 physics review — calibrated CMS Open Data $H\to ZZ^*\to4\ell$ measurement

**Verdict: ITERATE — not approved for physics advancement.**

The strategy has a sound intended measurement: it uses the requested
$105<m_{4\ell}<140$ GeV window, an external-$Z$ rather than Higgs-mass
calibration, a simultaneous category fit, and explicitly refuses to fake a
VBF tag from $p_\mathrm{T}^{4\ell}$.  Those are important strengths.  It is
not yet publication-credible because the VBF/FSR prerequisite is unproven
and the required DY+jets fake template has no quantitatively predictive
control-to-signal-region validation.  The Category A items below must be
resolved, and every Category B item must be addressed before a PASS.

## Review scope and evidence

This is an independent physics review of only the stated physics goal in
`prompt.md` and the Phase 1 artifact `outputs/STRATEGY.md`.  No methodology
or convention document was used.  The strategy itself states that Phase 1
creates no analysis figure, and it identifies no compiled PDF.  Therefore
there are no figures or data/MC ratio panels to inspect at this phase; the
mandatory visual-normalisation check is deferred until such plots exist.

The evidence cited below is by strategy section, decision, limitation, or
table so that each issue can be checked without inference from later work.

## Physics strengths

- The physics scope is appropriate: calibrated $m_{4\ell}$, a common
  $(m_H,\mu)$ extraction, the three final states, VBF/non-VBF categories,
  and an angular discriminant all directly address the requested result.
- The calibration is non-circular in principle.  [D3] uses separate
  $Z\to ee$ and $Z\to\mu\mu$ peak fits and an external $M_Z$ target rather
  than a Higgs mass or rate.  The calibration section correctly requires
  lepton-four-vector rebuilding before recomputing $m_{4\ell}$ and the
  angles, and it defines flavour correlations in $4e$, $4\mu$, and
  $2e2\mu$.
- The proposed angular NN uses the five conventional four-lepton-rest-frame
  angles, excludes $m_{4\ell}$ and event identity, and has sensible
  prerequisites: held-out agreement, input data/MC checks, mass-sculpting
  checks, and predeclared architecture/binning replicas.  The cut-based
  analysis remains an independent cross-check.
- The irreducible $q\bar q\to ZZ^{(*)}$ and $gg\to ZZ^{(*)}$ components are
  identified separately, and the strategy recognizes DY fakes and top/nonprompt
  contamination as reducible backgrounds.  The likelihood is correctly
  formulated as a simultaneous extended fit with shared parameters and
  correlated nuisances.
- The comparison plan names CMS HIG-16-041, preserves its requested mass-fit
  window, and treats the PDG Higgs mass as a validation target rather than a
  calibration input.  It also correctly warns that 10 fb$^{-1}$ cannot have
  the precision of the 35.9 fb$^{-1}$ reference.

## Findings

### A1 — Unproven raw-object path makes the binding VBF and FSR programme non-executable

**Classification: A — must resolve.**

**Evidence.** [A1] states that the delivered H4l tree has no jets, b tags,
missing momentum, or FSR photons, then requires future re-ntupling from
corresponding NanoAOD.  [D5] makes a reconstructed-jet VBF category binding,
while the shared reconstruction also requires FSR recovery.  The supplied
physics goal provides flat-tree locations but does not itself establish that
the matching raw data and every MC process are available, match the event
identifiers, and can be re-ntuplized with the required objects.

Calling the missing inputs a future prerequisite is not a feasibility
demonstration.  Without them, neither the requested VBF category nor
reference-like FSR candidate reconstruction can be performed; this cannot be
replaced by a proxy tag or an unvalidated omission.

**Required resolution.** Before advancing, demonstrate a pilot re-ntupling
path for data and all relevant MC processes, including stable event identity,
selected jet/FSR inputs, and the intended object calibrations.  Define the
jet tag before looking at the signal region and give an expected 10 fb$^{-1}$
yield, VBF efficiency, background composition, and purity/uncertainty study.
Then demonstrate that a simultaneous VBF mass category has non-trivial
information rather than a zero- or prior-dominated template.  The VBF
category remains required even if this check is difficult.

### A2 — The DY+jets fake template is not yet a predictive four-lepton background model

**Classification: A — must resolve.**

**Evidence.** [D6] and the backgrounds table designate `DYJetsToLL` MC as the
primary reducible/instrumental template.  The only stated protection is that
loose/not-tight and same-sign regions will “constrain its rate and test
transfer of its shape.”  The document gives neither a control-to-signal
transfer definition nor a closure/coverage criterion.  It also provides no
relative expected yield or explicit treatment/justification for other small
but potentially different reducible compositions such as electroweak
multilepton, conversion/$Z\gamma^{(*)}$, and top-associated processes.

A hard-process DY simulation alone does not establish detector fake-lepton
probabilities or the rare four-lepton tail.  A rate-only control constraint
cannot protect $m_H$ or $\mu$ from a shape/composition bias in the fit region.
This is particularly material because the NN is trained against the same
fake template.

**Required resolution.** Keep DY+jets MC as the required nominal model, but
make it physically testable: predefine fake-enriched regions, the transfer
from each region to the signal selection, independent closure tests, and
rate-and-shape nuisance parameters with correlations in the simultaneous
fit.  Inventory all reducible sources and either include a template/control
constraint or quantitatively demonstrate their negligible impact.  Test the
NN-score as well as the $m_{4\ell}$ transfer; a visible control-region shape
failure must not be absorbed into an unconstrained systematic.

### B1 — Data/MC calibration treatment is ambiguous at the template level

**Classification: B — should address before PASS.**

**Evidence.** The calibration section defines
$c_f=M_Z^{\rm ref}/\widehat m_{Z,f}^{\rm data}$, fits MC separately only as a
validation, and then says to apply corrections before rebuilding “templates.”
It does not explicitly state whether the nominal and varied corrections are
applied to data only, to MC using its own fitted response, or through a
data--MC residual transformation.  These choices are not equivalent when the
reconstructed MC $Z$ peak is displaced from its generated/reference scale.

The intended calibration is independent of $m_H$, which is good, but an
ambiguous data/template convention can still bias a mass fit or make a scale
nuisance artificially cancel.

**Required resolution.** Specify one data and MC calibration convention,
including the role of the separately fitted MC peak, and persist the exact
nominal and shifted transformations applied to every fitted template.  Show
with injected known scale offsets that the procedure recovers the offset and
has valid coverage in each final state.

### B2 — The residual-scale programme lacks the requested quantitative 10 fb$^{-1}$ sanity test

**Classification: B — should address before PASS.**

**Evidence.** The scale plan appropriately lists peak-fit covariance,
line-shape, binning/domain, residual, and supported kinematic variations,
with sensible electron/muon correlations.  However, the only numerical
benchmark in [L1] is an extrapolation of the *Higgs mass statistical*
uncertainty.  The strategy says Phase 4 will compare the profiled scale
impact to [L1], but it defines no expected scale-impact calculation from the
available $Z$ control statistics or residual calibration precision.

Consequently a quoted post-fit scale error could be implausibly small, or a
large one could conceal a defective calibration, without a predeclared
physics check.  Treating several overlapping fit alternatives as independent
nuisances would also risk double counting unless their covariance/envelope
rule is defined.

**Required resolution.** Predeclare a numerical calibration-to-$m_H$
sanity calculation for 10 fb$^{-1}$: use the measured/expected $Z$-fit
precision and injected scale shifts to predict the propagated mass effect,
then verify it with pseudoexperiments and coverage.  Define how fit-model
alternatives become a covariance or a single justified envelope rather than
multiple double-counted errors.  Compare with the scaled reference systematic
only where the differing reconstruction and luminosity make that comparison
meaningful.

### B3 — The simultaneous fit needs an explicit sparse-category and resolving-power preflight

**Classification: B — should address before PASS.**

**Evidence.** The fit crosses three flavours with VBF/non-VBF selections and
adds NN-score bins, yet the exposure is only 10 fb$^{-1}$.  The strategy says
the VBF score dimension is used only if populated and promises category
occupancy in a decision table, but it gives no predeclared minimum expected
counts, mass/score binning relative to resolution, MC-statistical treatment
for sparse templates, or toy criterion for identifying both $m_H$ and $\mu$.

The likelihood form is appropriate, but a formally simultaneous likelihood
can still simply reproduce template priors if its categories and NN bins are
mostly empty.  This is a physics sensitivity issue, not an implementation
detail.

**Required resolution.** Before unblinding, publish the expected yield and
background hierarchy in every final-state/category/score bin; predeclare
binning at a scale demonstrated to retain mass information; and use
independent pseudo-data to show unbiased $(m_H,\mu)$ recovery, finite pulls,
and sensitivity to at least a ±20% signal-strength change.  Preserve the
required VBF category while making any score-binning decision reproducible.

### B4 — Comparison coverage is incomplete for $\mu$ and does not establish independence

**Classification: B — should address before PASS.**

**Evidence.** The reference table gives an HIG-16-041 value for both $m_H$
and $\mu$ and a PDG value only for $m_H$.  It promises a covariance-aware
compatibility but does not identify a global/world signal-strength comparator
for $\mu$, or explicitly document why none is meaningfully comparable.  It
also calls HIG-16-041 directly comparable without determining whether the
10 fb$^{-1}$ sample overlaps that 35.9 fb$^{-1}$ data set; an overlap would
make an apparent compatibility calculation correlated rather than an
independent validation.

**Required resolution.** For each reported quantity, name the comparison
target or state precisely why no target is comparable.  Establish data-set
overlap and include its covariance where known; otherwise report a descriptive
difference rather than an independence-implying pull/p-value.  Retain the
required HIG-16-041 comparison and clearly separate a world-average
validation target from calibration input.

### C1 — Complete the common reference-scale correlation

**Classification: C — suggestion.**

The external $M_Z$ uncertainty is tiny compared with the expected 10 fb$^{-1}$
mass uncertainty, but it is common to the electron and muon calibrations.
Include it as a fully correlated reference component (or explicitly quantify
and document its negligible effect) so the stated covariance is complete.

### C2 — Strengthen the NN score validation beyond individual inputs

**Classification: C — suggestion.**

The planned input-level data/MC checks and mass-sculpting studies are good.
Also validate the calibrated *score* in sideband and fake-enriched data with
the physical process mixture used for training, and record the class-weight
choice.  That makes a disagreement in the learned combination diagnosable
even when each one-dimensional angular input appears acceptable.

## Requested-topic assessment

| Topic | Physics assessment |
|---|---|
| Backgrounds | Irreducible components are correctly separated, but the reducible/fake programme is not yet predictive enough for a mass/rate fit (A2). |
| Calibration independence | Sound in intent: the $Z$ mass, not $m_H$ or $\mu$, anchors the correction.  The exact data/MC template convention must be fixed (B1). |
| VBF feasibility | Correctly refuses an invalid proxy, but its raw-object and 10 fb$^{-1}$ viability are unproven (A1, B3). |
| Angular NN | Physically sensible five-angle input and useful gates; score-level fake/control validation should be added (A2, C2). |
| DY fake model | Satisfies the requested nominal model only conditionally; it needs a defined transfer, closure, and composition treatment (A2). |
| Simultaneous $m_{4\ell}$ fit | Appropriate structure and shared correlated nuisances; its sparse-template resolving power remains unshown (B3). |
| Reference comparisons | HIG-16-041/window and PDG-$m_H$ treatment are good; $\mu$ coverage and sample correlation must be resolved (B4). |
| Scale nuisances | Broad source list and flavour correlations are credible; calibration/template convention and a numerical 10 fb$^{-1}$ validation are required (B1, B2). |

Permission to pass through a later Phase 4b human gate does not relax these
technical requirements.  A Phase 4b gate cannot validate a VBF category that
cannot be constructed or a fake model that has not demonstrated predictive
closure.
