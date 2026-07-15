# Phase 1 critical review — calibrated CMS Open Data $H\to ZZ^*\to4\ell$

**Verdict: ITERATE — do not advance Phase 1.**

This is an evidence-based review of the committed strategy only. It does not
ask for Phase 2/3 validation results prematurely: no Phase 1 figures, fit
outputs, or compiled PDF exist, so visual data/MC and closure *outcomes* are
not expected yet. The findings instead identify plans that are too undefined
to make binding Phase 1 commitments executable or auditable. All Category A
and B items require resolution and an independent re-review before PASS.

## Scope and evidence read

Read: prompt.md, root and Phase 1 instructions, plan.md, experiment_log.md,
retrieval_log.md, outputs/STRATEGY.md, both candidate conventions,
methodology/03-phases.md Phase 1, and methodology/06-review.md §§6.3–6.4.
I also checked the cited local ntuplizer claims. h4l_ntuplizer.cpp:123-148
defines the flat output as event/PV, candidate, and four lepton records, and
its reader at :528-602 requests event/PV, trigger, muon, and electron fields
but no jet, b-tag, FSR-photon, or MET collection. Conversely, the dilepton
output retains event identifiers and lepton quality fields
(dilep_ntuplizer.cpp:122-145,348-414).

The corpus connector limitation is correctly documented in
retrieval_log.md:3-8. I did not wait for unavailable corpus tools or undertake
new external literature research. The strategy has direct official source links
at STRATEGY.md:404-415; the findings distinguish citation presence from
unrecorded semantic verification of the cited claims.

## What is already sound

- **Technique and convention applicability.** [D1]/[D9] correctly define a
  detector-level profile-likelihood measurement rather than a closed-form
  double-tag result or unfolded particle-level spectrum
  (STRATEGY.md:64-67,105-107). The non-applicability reasons for core unfolding
  deliverables agree with conventions/unfolding.md:6-30. The extraction table
  maps analogous sources, which is a useful completeness check
  (STRATEGY.md:316-352).
- **Requested fit and calibration direction.** [D2] commits every nominal and
  systematic fit to $105<m_{4\ell}<140$ GeV (STRATEGY.md:69-72), and [D3] uses
  separate external-$Z$ electron and muon calibration rather than forcing the
  Higgs result (:74-77,201-240). The intended lepton-level recomputation before
  $m_{4\ell}$, angular variables, and templates is correct (:227-230).
- **Two genuine selection approaches.** The cut-only analysis and angular-NN
  fit are qualitatively different and both are retained, with expected
  precision and goodness of fit in the comparison (STRATEGY.md:165-197). The
  NN excludes $m_{4\ell}$ and specifies the five requested rest-frame angles
  (:141-154).
- **Input normalization and a numerical check.** The sample cross sections in
  STRATEGY.md:127-133 reproduce prompt.md:25-35, and the Metadata.nEvents
  denominator agrees with prompt.md:37-39 (STRATEGY.md:111-121). The [L1]
  arithmetic is correct: $0.20\sqrt{35.9/10}=0.378946$ GeV, rounded to
  0.38 GeV (:52-57). It is a mass-statistics benchmark, not the requested
  lepton-scale sanity calculation.

## Category A — must resolve

### A1 — [A1]/[A2]/[D5] do not define an executable raw-object re-ntupling gate

**Evidence.** The strategy accurately says that the delivered flat H4l tree
lacks objects needed for VBF and FSR, then defers a replacement NanoAOD
ntupling (STRATEGY.md:31-38). [D5] nevertheless makes a reconstructed-jet VBF
category binding (:84-88), while common reconstruction requires FSR recovery
before candidate masses (:157-163). The user supplied only flat H4l and
dilepton paths (prompt.md:17-23); plan.md:11-14 proposes only future raw-input
inspection. No raw file location, process-by-process availability check, pilot
output, event-identity comparison, jet/FSR object definition, or predeclared
VBF tag exists. There is also a timing ambiguity: [A1] says re-ntupling must
happen *before Phase 3* (:34-36), but the handoff assigns raw-object VBF/FSR
implementation to Phase 3 (:388-392).

This is not an implementation detail: without the stated objects, the
required VBF category and reference-like FSR reconstruction cannot be built. A
$p_{\mathrm T}^{4\ell}$ proxy would violate the strategy's own constraint.

**Required resolution.** Make [A1] a hard, ordered gate before selection
optimization: locate matching raw NanoAOD for data and every relevant MC
process; run a small re-ntupling pilot; persist run/lumi/event, selected jets,
b tags, FSR candidates, and calibration inputs; and demonstrate event/candidate
matching to the flat tree where applicable. Predeclare jet cleaning, VBF
variables/thresholds, FSR association/recovery, and a viability test (expected
10 fb$^{-1}$ yield, purity/background composition, and non-prior-dominated fit
information). If raw data cannot be demonstrated, formally revise the strategy;
do not silently omit VBF.

### A2 — [D6] is a nominal DY label, not a predictive fake-background and background-completeness plan

**Evidence.** The background table has correct broad classifications, but no
relative expected yields or hierarchy despite the Phase 1 requirement to
estimate importance (STRATEGY.md:242-250; phase1_strategy/CLAUDE.md,
"Enumerate backgrounds"). DYJetsToLL is the nominal fake template while
$t\bar t$ is grouped with an undefined "other nonprompt" component
(STRATEGY.md:249-250). The only transfer prescription is that loose/not-tight
and same-sign regions will constrain/test it (:252-256): there is no
loose/tight definition, region orthogonality, transfer observable/factor,
independent closure sample, acceptance criterion, or rate/shape correlation
model. This same DY model is used in NN training (:175-183), so an untested
fake shape can bias both score and mass fit.

conventions/extraction.md:85-104 requires explicit contamination and
composition treatment whenever its analogous ingredients are applicable; the
strategy's own convention audit marks both applicable (STRATEGY.md:327-330).
A generic "other" row is not a disposition for
electroweak/conversion/heavy-flavour/nonprompt contributions: they need not be
invented as samples, but must be inventoried and either constrained, shown
negligible with uncertainty, or formally unavailable.

**Required resolution.** Preserve DY+jets MC as the user-mandated nominal
model, but predeclare fake-enriched regions, loose/tight object definitions,
control-to-signal transfer, independent closure/coverage tests and pass
criteria, and separate correlated rate/shape nuisances in $m_{4\ell}$ and
$D_{\rm NN}$. Provide expected yields/background hierarchy by flavour/category
and a disposition for each reducible source. A control-region shape failure
must propagate to the template uncertainty rather than be absorbed by an
unconstrained normalization.

### A3 — the binding systematic and commitment ledger is incomplete

**Evidence.** Phase methodology requires COMMITMENTS.md at Phase 1 completion,
listing every systematic, cross-check, validation, comparison, and planned
figure with a machine-readable status (methodology/03-phases.md:771-784). The
executor plan promised it (plan.md:17-22), but it is absent from the worktree.
This removes the formal downstream trace for [A]/[L]/[D] decisions.

The systematic table has applicable sources without a final "will implement"
or justified "not applicable" disposition. Luminosity is "blocked by [A3]"
(STRATEGY.md:298) although $\mu$ is primary, and signal acceptance/PDF/QCD-scale
plus shower/hadronization/underlying-event rows are merely "inventory required
by [L2]" (:301-302). The extraction analogue calls relevant efficiency, model,
physics-parameter, and production-mixture sources applicable (:324-332). The
Phase 1 rule requires an actual disposition for every applicable source, not an
open inventory task; methodology/06-review.md:301-313 makes unjustified
completeness omissions Category A.

**Required resolution.** Create COMMITMENTS.md and list each [A], [L], [D],
systematic, validation, reference comparison, and flagship figure with explicit
status/action. Set a decision tree for luminosity: obtain a dataset-appropriate
prior, or define exactly how a conditional-on-10-fb$^{-1}$ $\mu$ is represented
and reported without a false luminosity systematic. For generator/model
sources, commit to specific in-file weights, alternative configurations, or an
evidenced formal downscope after inventory; do not leave "inventory required"
as the final strategy state.

### A4 — published-method parity is asserted but not demonstrated or made a cross-check

**Evidence.** Phase 1 requires comparison to the *statistical extraction* in
the references and either the same/better method or a justified simpler method
with the published method committed as a cross-check
(methodology/03-phases.md:211-231). The strategy chooses a binned template
likelihood (STRATEGY.md:262-284) and the reference table says only "match ...
likelihood treatment" (:356-360). It gives no source-backed account of the
HIG-16-041 likelihood dimensionality, per-event resolution treatment,
signal-shape treatment, or which feature makes the binned approximation
equivalent at 10 fb$^{-1}$. Nor is a reference-like fit listed as a cross-check.
The same table identifies per-event mass and angular information in the earlier
reference (:359), underscoring why a generic statement of parity is
insufficient.

**Required resolution.** Add a source-specific method-parity row: state the
published likelihood inputs, category/discriminant and resolution treatment,
then either reproduce each material feature or quantify a predeclared
binned-versus-reference-like pseudoexperiment cross-check. Test bias and
coverage for both $m_H$ and $\mu$, not merely nominal-template self-consistency.

## Category B — must fix before PASS

### B1 — the $Z$ control calibration is not fully specified at the data/MC template level

**Evidence.** The strategy calls for a clean control selection but leaves
"tight flavour-appropriate ID/isolation," primary-vertex compatibility, and a
"predeclared broad $Z$-pole fit domain" without numerical/source-selected
working points or a background model (STRATEGY.md:201-210). This matters because
the supplied dilepton ntuplizer deliberately has **no** $m_{\ell\ell}$
requirement and retains the pair closest to $M_Z$
(dilep_ntuplizer.cpp:242-269); its loose baseline is only at :548-552.
The strategy defines $c_f=M_Z^{\rm ref}/\hat m^{\rm data}_{Z,f}$ but says MC
peak fits merely validate response, then requires rebuilding all templates
(STRATEGY.md:211-230). It never fixes whether nominal/varied transformations
apply to data only, MC with its own response, or a defined data/MC residual
scheme. Those choices can change fitted mass or make a scale nuisance cancel.

**Required resolution.** Specify Z-control trigger/object/vertex/mass
selection, signal-plus-background line shape, fit range, and a predeclared
goodness-of-fit/closure criterion. Fix nominal and varied data/MC
transformations, retain fitted MC response, and test injected electron/muon
offsets in $4e$, $4\mu$, and $2e2\mu$. State quantitative occupancy and
nested-test criteria before permitting $\eta$/$p_T$ subdivisions.

### B2 — residual scale variations have no non-overlapping combination rule or 10 fb$^{-1}$ numerical sanity target

**Evidence.** Peak-fit covariance, model, binning/domain, residual, and
kinematic alternatives are listed (STRATEGY.md:232-240), but the document does
not define whether overlapping alternatives form one envelope, fitted covariance,
or independent nuisances. It also gives no expected Z-peak precision, injected
scale study, pseudoexperiment coverage target, or calibration-to-$m_H$ estimate.
[L1] instead scales the *published mass statistical* uncertainty (:52-57),
which cannot test the requested lepton-scale contribution at 10 fb$^{-1}$.

**Required resolution.** Before the fit is built, predeclare a non-double-
counted covariance/envelope construction and produce the numerical chain:
expected/measured Z fit precision $\to$ electron/muon scale variations
$\to m_{4\ell}$ shift $\to$ profiled $(m_H,\mu)$ contribution. Validate with
injected shifts and pseudoexperiment coverage, retaining stated flavour and
category correlations.

### B3 — angular-NN gates need quantitative decision criteria and complete background treatment

**Evidence.** The five-angle input set, split strategy, and qualitative
validation list are good (STRATEGY.md:175-197), but no threshold is given for
input data/MC agreement, score–mass correlation, overtraining, or acceptable
score-bin population. Primary training backgrounds are ZZ and DY only
(:177-180), whereas the fit also includes top/other nonprompt backgrounds
(:249-250); their classifier treatment is unstated. The promised formal
strategy revision if a gate fails ([D4], :79-82) is not reproducible.

**Required resolution.** Predeclare score validation statistics and acceptance
limits, class mixture/weights and disposition of all fit backgrounds, a
mass-sculpting metric, score-bin occupancy rule, and exact fallback/revision
trigger. Validate the *score* in sideband and fake-control mixtures, not only
each one-dimensional angular input.

### B4 — the simultaneous fit has no sparse-category, morphing, or joint-coverage preflight

**Evidence.** Six flavour×VBF categories and optional score bins share
$(m_H,\mu)$ and nuisance parameters (STRATEGY.md:258-284), yet the only
sparsity safeguard is "if expected population supports it" (:181-183). There
are no predeclared expected yields by bin, mass/score binning versus resolution,
mass-hypothesis grid/morphing closure metric, MC-statistical treatment for
sparse templates, or independent pseudo-data coverage test. The unfolding
analogue adopts split-MC and injected tests (:337-352) but does not operationalize
them for this likelihood.

**Required resolution.** Publish expected yields/background hierarchy per
final-state/category/score bin, define binning and morphing tests, and use
independent pseudo-data to demonstrate unbiased joint recovery, finite pulls,
healthy goodness of fit, and sensitivity to a predeclared $\pm20\%$ change in
$\mu$. Preserve VBF while making score-bin retention deterministic.

### B5 — comparison coverage and comparability are incomplete, despite direct primary citations

**Evidence.** The strategy gives the requested HIG-16-041 fit window and
$m_H$/$\mu$ targets, and labels PDG Higgs mass validation rather than
calibration input (STRATEGY.md:354-367). However, it calls HIG-16 directly
comparable without establishing whether supplied 10 fb$^{-1}$ data overlap the
35.9 fb$^{-1}$ reference, which determines whether a pull is independent. It
gives no world/comparable $\mu$ target or explicit explanation that none
exists, and names HEPData for later jet comparisons without mapping each
planned key distribution to a table, binning, covariance, or documented reason
it is not comparable. The retrieval record documents unavailable corpus tools
but contains no page/table or HEPData-resource provenance for external numeric
claims (retrieval_log.md:3-8; STRATEGY.md:404-415).

**Required resolution.** Add a result-by-result comparison matrix for $m_H$,
$\mu$, and material published distributions: target, observable/selection
equivalence, data-overlap/covariance treatment, comparison statistic, and
precise reason where no comparison is meaningful. Record publication
section/table or HEPData resource for every numerical target. Do not present
an overlap-correlated comparison as independent validation.

## Category C — suggestions

### C1 — retain the common external-$M_Z$ component explicitly

The PDG $M_Z$ uncertainty is likely negligible, but it is common to electron
and muon calibration equations (STRATEGY.md:74-77,211-219). Quantify it once
and either include it as a fully correlated component or state its demonstrated
negligible impact.

### C2 — make the fallback source record audit-friendly

Direct official URLs in STRATEGY.md:410-415 are a good fallback after the
connector limitation. Add one retrieval-log row per external number with source
page/table/HEPData table identifier and access date so a later reviewer can
verify [R1]–[R5] without reconstructing the literature search.

## Binding-label and requirement audit

| Commitment / requirement | Evidence | Review result |
|---|---|---|
| [D1]/[D9] template fit, not extraction/unfolding | STRATEGY.md:64-67,105-107,314-352 | Direction correct; published-method parity remains A4. |
| [D2] reference fit window | STRATEGY.md:69-72 | Present. |
| [D3] independent $Z$ calibration and full propagation | STRATEGY.md:74-77,201-240 | Intent matches user scope; operational detail/scale validation are B1–B2. |
| [D4] angular rest-frame NN plus cut cross-check | STRATEGY.md:141-154,165-197 | Two approaches present; reproducible NN gates are B3. |
| [A1]/[A2]/[D5] VBF/FSR re-ntupling | STRATEGY.md:31-44,84-88,157-163,388-392; ntuplizer evidence above | A1. |
| [D6] DY MC fake model and transfer/closure | STRATEGY.md:90-92,242-256,305 | A2. |
| [D7] staged pseudo-data/hashed tenth | STRATEGY.md:94-98 | Present; later validation must be substantive. |
| [D8] supplied cross section/Metadata normalization | prompt.md:25-39; STRATEGY.md:100-121 | Present and internally consistent. |
| Systematic enumeration and downstream traceability | STRATEGY.md:286-312; missing COMMITMENTS.md | A3. |
| Reference analyses, numerical targets, comparisons | STRATEGY.md:354-367,404-415 | Three references/direct links present; parity A4, comparability B5. |
| Six flagship figures / diagrams | STRATEGY.md:369-384 | Present as Phase 5 commitments; put them in the required ledger. |

## Required re-review evidence

The next critical review should receive: COMMITMENTS.md; a raw-object
re-ntupling pilot/gate record; formal fake-region/transfer/closure specification
and expected background yields; a method-parity table; complete data/MC
calibration and residual-scale validation plan; NN and simultaneous-fit
pseudoexperiment criteria; and a comparison matrix with source provenance. A
rewritten strategy alone is insufficient unless these artifacts supply the
listed evidence.
