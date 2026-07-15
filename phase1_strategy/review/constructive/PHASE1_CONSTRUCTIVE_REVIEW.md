# Phase 1 constructive review — calibrated CMS Open Data $H\to ZZ^*\to4\ell$ measurement

**Verdict: ITERATE — promising and substantially aligned, but not yet credible
enough to advance.**

The strategy makes the right high-level choices: it preserves the required
mass window, refuses a fake VBF proxy, uses an independent $Z$ calibration,
keeps a cut-based cross-check, and commits to a simultaneous likelihood.  The
central concern is feasibility rather than intent.  Four required outcomes are
only asserted as future work without a falsifiable gate: raw-object VBF/FSR
recovery, predictive DY-fake transfer, the numerical 10 fb$^{-1}$ scale
sanity test, and durable tracking of the binding commitments.  Resolving those
will make the otherwise strong plan executable rather than aspirational.

## Scope, evidence, and review limit

This is an evidence-based review of `prompt.md`, the committed
`outputs/STRATEGY.md`, Phase 1/review methodology, both candidate conventions,
the phase experiment/retrieval logs, and a bounded read of the two supplied
ntuplizers.  The corpus connector is unavailable as recorded in
`retrieval_log.md`; consequently, this review does **not** independently
verify the numerical values in the cited CMS/PDG sources.  It assesses whether
the strategy contains a credible and testable plan, not whether later
implementation has already succeeded.

The bounded source inspection supports two important claims in the strategy:
the H4l ntuplizer persists four selected lepton $p_T,\eta,\phi,m$, charge,
flavour, $Z$ assignment, and run/lumi/event identifiers
(`ntuplizer/h4l_ntuplizer.cpp:123-130, 412-478, 529-531`), so an angular
calculation is feasible from candidate-level leptons; it has no `Jet_*`,
b-tag, or photon/FSR input/output handling.  Thus the angular-NN plan is not
blocked by the current tree, whereas the VBF/FSR plan genuinely requires a
new raw-object path.

## Constructive strengths

- The endpoint is correctly constrained to calibrated detector-level
  $(m_H,\mu)$ extraction rather than an unsupported unfolded claim
  (strategy [D1], [D9], lines 64-67 and 105-107).  The fixed
  $105<m_{4\ell}<140$ GeV window is explicitly binding (lines 69-72), as the
  task requires.
- The NN is meaningfully part of the proposed inference, not decorative.  It
  uses five specified rest-frame angles, excludes $m_{4\ell}$ and event
  identity, retains a cut-based alternative, and requires held-out,
  input-modelling, sculpting, and replica tests (lines 141-155 and 175-197).
- The calibration is independent in principle: separate $ee$/$\mu\mu$ peaks
  anchor lepton corrections to $M_Z$, never to a fitted Higgs mass or rate
  (lines 199-225).  Rebuilding four-vectors, masses, angles, scores, and
  templates after the correction (lines 227-240) is the right propagation
  order.
- The strategy has a real simultaneous extended-likelihood structure with
  shared $m_H$, $\mu$, and correlated nuisances, rather than combining
  independently fitted channels after the fact (lines 258-284).
- It acknowledges the 10 fb$^{-1}$ statistical limitation numerically
  (0.38 GeV planning benchmark; lines 52-57) without presenting it as an
  observed result.  It also distinguishes the PDG Higgs mass as validation,
  rather than calibration input (lines 362-367).

## Requested-feature audit

| Requested feature | Evidence in the committed strategy | Assessment and required gate |
|---|---|---|
| Raw-object VBF category | [A1]/[D5] correctly reject a $p_T^{4\ell}$ proxy and require a jet-based, exclusive category (lines 31-38, 84-88). | **Blocked (A1).** A future re-ntupling assertion is not a demonstrated raw-data/MC path. |
| Four-lepton-rest-frame angular NN | Five input angles, train/test split, no mass input, a cut-based cross-check, and pre-fit validation are specified (lines 141-155, 175-197). | **Conditionally credible (B2).** Candidate lepton four-vectors exist, but score-level transfer criteria and a quantified decision rule are missing. |
| Flavour-separated $Z$ calibration | Clean-control intent, separate peak fits, a scale equation, recomputation, and flavour correlations appear at lines 199-240. | **Partly specified (B1, A4).** The data/MC template convention, fit specification, covariance/envelope rule, and required numerical scale sanity test need gates. |
| DY+jets fake model | DY is the required nominal template, with loose/not-tight and same-sign controls named (lines 242-256, 305). | **Blocked (A2).** No defined transfer factor, closure, composition inventory, or shape-coverage condition exists. |
| Simultaneous fit for $m_H$ and $\mu$ | Six nominal flavour/category cells, score bins when populated, extended likelihood, template morphing, and profile outputs are specified (lines 258-284). | **Conditionally credible (B3).** It needs sparse-bin rules and independent pseudoexperiment coverage/identifiability requirements. |
| CMS/PDG comparisons | Three references and several numerical targets are tabulated; the PDG $m_H$ comparison is named (lines 354-367). | **Partly specified (B4).** A result-by-result metric and data-overlap/correlation treatment are absent, especially for $\mu$. |
| 10 fb$^{-1}$ scale sanity | The only stated numeric benchmark is a scaled *statistical* mass uncertainty (lines 52-57); Phase 4 is asked only to compare the profiled scale contribution to it (lines 232-240). | **Blocked (A4).** This does not test whether the *lepton-scale* uncertainty is physically reasonable at 10 fb$^{-1}$. |

## Findings

### A1 — Raw-object recovery is a prerequisite, not yet an executable VBF/FSR plan

**Classification: A — must resolve.**

**Evidence.** The strategy explicitly says that the delivered flat tree has
no jets, b tags, MET, or FSR photons and requires a later re-ntupling from
“corresponding NanoAOD” (lines 31-38).  It nevertheless makes a reconstructed
jet VBF category binding (lines 84-88) and requires FSR before reference-like
candidate reconstruction (lines 157-163).  The local ntuplizer corroborates
the object deficit: it reads only muon/electron collections in the relevant
input setup (`ntuplizer/h4l_ntuplizer.cpp:549-602`), even though it preserves
event identifiers.

**Why this matters.** The required VBF result cannot be recovered from a
candidate-level proxy.  Without a verified raw source for data and each needed
MC process, the analysis cannot establish jet calibration, VBF/non-VBF
migration, b-veto/top-control behaviour, FSR recovery, or event matching.

**Required upgrade / exit evidence.** Before Phase 2 closes, provide a pilot
raw-object re-ntupling record that names accessible input locations for data
and all fit MC, preserves run/lumi/event matching, and writes selected jets,
JEC/JER inputs, b-tag inputs, FSR candidates, and their systematic branches.
Predeclare the VBF selection and demonstrate, on independent MC, its 10
fb$^{-1}$ expected yield, signal efficiency, background composition, purity,
and migration uncertainty.  Keep a VBF mass category even if it is sparse;
do not replace it with a proxy.

### A2 — The DY+jets fake template has no predictive transfer or closure contract

**Classification: A — must resolve.**

**Evidence.** DY+jets is appropriately retained as the user-mandated nominal
fake template (lines 249-256), but its controls are only named as
“loose/not-tight and same-sign” and described generically as rate/shape tests.
The systematic table lists “control-statistics, composition, and transfer
tests” (line 305), without defining a signal-selection transfer, a closure
sample, an expected fake hierarchy, a pass criterion, or a score-bin test.

**Why this matters.** A hard-process DY sample does not by itself validate the
rare four-lepton fake tail.  A rate constraint cannot protect the mass fit
from an $m_{4\ell}$, final-state, VBF, or NN-score shape bias.  The same model
is also used in NN training (lines 177-183), so an untested transfer can bias
both selection and fit.

**Required upgrade / exit evidence.** Define the tight/loose and same-sign
regions, disjoint validation region(s), and the transfer function from each
control selection to every final-state/category/score template before seeing
the signal region.  For each transfer, report expected composition and yield,
rate and shape closure with a predeclared statistic/p-value or coverage band,
and rate/shape nuisance correlations in the workspace.  Inventory
$t\bar t$, conversion/$Z\gamma^{(*)}$, and any electroweak multilepton
component; either model/control each or quantify a negligible propagated
impact.  A failed control-shape test must drive a resolved model change or a
bounded, validated shape variation—not an unconstrained catch-all error.

### A3 — The mandatory commitment-tracking artifact is absent

**Classification: A — must resolve.**

**Evidence.** `methodology/03-phases.md` requires a Phase 1
`COMMITMENTS.md` that tracks every systematic, cross-check, validation,
comparison target, and planned figure with status.  `phase1_strategy/plan.md`
also promised that file, yet it is absent from the phase directory.  The
strategy contains binding [A], [L], and [D] labels (lines 29-107), but prose
alone cannot provide the machine-readable phase-to-phase audit required by
the methodology.

**Why this matters.** The strategy intentionally makes VBF, the NN,
calibration, fake controls, fit validation, and comparison commitments
binding.  Without a tracked status for each, later phases can silently omit a
hard requirement while still appearing internally consistent.

**Required upgrade / exit evidence.** Create `phase1_strategy/COMMITMENTS.md`
before PASS.  Include each [A]/[L]/[D], every systematic row, both selection
approaches, every fake/fit/calibration validation, the reference comparison
matrix, and the six flagship figures.  Use the required `[ ]`, `[x]`, or `[D]`
status semantics and update it at every phase boundary.

### A4 — The strategy lacks the required numerical 10 fb$^{-1}$ lepton-scale sanity test

**Classification: A — must resolve.**

**Evidence.** [L1] scales HIG-16-041’s *statistical* $m_H$ uncertainty to
0.38 GeV (lines 52-57).  The residual-scale section lists sensible sources
and says that Phase 4 will report the profiled contribution (lines 232-240),
but it gives no numerical calculation linking the achievable $Z$-peak
precision to the propagated $(m_H,\mu)$ impact.  The user explicitly asks for
this check rather than a nominal $m_{4\ell}$ statement.

**Why this matters.** A mass-statistics extrapolation cannot establish that a
quoted lepton-scale uncertainty is plausible.  It could mask an
under-constrained peak fit, an over-large envelope, or an accidental
data/template cancellation.

**Required upgrade / exit evidence.** Predeclare a scale-only 10 fb$^{-1}$
test: obtain expected or measured separate $Z\to ee$ and $Z\to\mu\mu$ peak-fit
covariances, convert them to $c_e,c_\mu$ variations, rerun independent
pseudoexperiments/templates, and report the induced $\Delta m_H$ and
$\Delta\mu$ with pull/coverage.  State an expected numerical range before
the final fit and compare the observed profiled contribution with it.  Treat
fit-model/binning/domain alternatives as a documented covariance or one
justified envelope, rather than summing overlapping alternatives as
independent nuisances.

### B1 — Calibration control selection and data/MC template convention need a reproducible definition

**Classification: B — must fix before PASS.**

**Evidence.** The clean control selection is qualitative—trigger matching,
OS-SF, tight ID/isolation, vertex compatibility, and a “predeclared broad
$Z$-pole fit domain” (lines 201-206).  The nominal correction is defined from
the data peak (lines 208-225), MC peaks are called validation, and then all
“templates” are rebuilt (lines 227-230).  The artifact never states whether
data and MC receive separate absolute corrections, a data/MC residual
correction, or another common-scale convention.

**Required upgrade / exit evidence.** Specify the fit domain, background and
resolution models, fit-quality/closure condition, trigger and ID working
points, and a predeclared occupancy/nested-test rule for $\eta$ or $p_T$
subdivision.  Define the exact nominal and varied transformation for data and
for each MC template, then validate it by injecting known electron/muon scale
offsets and demonstrating recovery and coverage in $4e$, $4\mu$, and
$2e2\mu$.  Retain the common $M_Z$ reference uncertainty as correlated or
show quantitatively that it is negligible.

### B2 — The angular-NN plan needs a quantitative score-transfer and decision contract

**Classification: B — must fix before PASS.**

**Evidence.** The strategy specifies the five angles and good qualitative
checks (lines 145-155 and 185-197), and the lepton four-vectors needed to
calculate them do exist in the current output source.  However, it gives no
predeclared acceptance criterion for input/score data--MC agreement,
mass-sculpting, or meaningful separation, and it does not require a
mixture-aware score validation in a fake-enriched region.  “B is primary if
these gates pass” is not reproducible without thresholds (line 193).

**Required upgrade / exit evidence.** Fix the pairing/angle convention test
against a reference implementation and establish pre-unblinding gates for
input and score agreement, mass correlation/sculpting, held-out performance,
and architecture/seed spread.  Validate the *score*, not only one-dimensional
inputs, in sideband and fake-control mixtures with the same class weights and
process composition as training.  The selection table should state which
metric decides B versus A and show that a smaller expected uncertainty cannot
win if the modelling/GoF gate fails.

### B3 — The simultaneous fit needs sparse-category, identifiability, and coverage gates

**Classification: B — must fix before PASS.**

**Evidence.** Six flavour-by-VBF categories and optional NN-score bins are
planned (lines 258-260), while the only numerical 10 fb$^{-1}$ check is a
global mass-statistics extrapolation.  The likelihood and morphing plan are
sound (lines 262-284), but there are no expected per-bin yields, minimum
effective template-statistics rule, mass/score binning choice, independent
pseudoexperiment pull/coverage requirement, boundary rule, or test of
sensitivity to a non-nominal $\mu$.  The convention audit calls split-MC and
injected tests “mandatory” (lines 337-348) but does not define their
acceptance evidence.

**Required upgrade / exit evidence.** Before data fitting, publish an Asimov
and independent-pseudoexperiment preflight table for each final-state,
VBF/non-VBF, and score bin: expected signal/background yields, MC-statistical
treatment, chosen binning relative to mass resolution, and nuisance impact.
Inject non-nominal $m_H$ and at least a ±20% $\mu$ change, demonstrate unbiased
recovery with finite pull widths/coverage and no boundary saturation, and use
toy-based GoF.  Keep the VBF category binding while using score bins only when
the predeclared information/occupancy criterion is met.

### B4 — Comparison commitments need a quantity-by-quantity metric and independence statement

**Classification: B — must fix before PASS.**

**Evidence.** The reference table gives HIG-16-041, earlier CMS, and a modern
CMS mass result with useful comparability prose (lines 354-367).  It promises
a covariance-aware compatibility for $m_H$, but does not define the statistic
or reference covariance.  It gives no world/combined $\mu$ target or explicit
statement that no such comparison is meaningful.  Nor does it establish
whether the supplied 10 fb$^{-1}$ data overlap the 35.9 fb$^{-1}$ reference,
which determines whether a pull can be described as independent.

**Required upgrade / exit evidence.** Add a comparison matrix covering each
reported quantity and observable: target, observable/energy/data-overlap
status, correlation assumption, metric (difference, pull, or
covariance-aware $\chi^2$), and the reason when no comparison is meaningful.
For $m_H$, distinguish PDG validation from a calibration input; for $\mu$,
retain HIG-16-041 and either identify a valid reference or explicitly state
why a world average is not comparable.  Use descriptive differences rather
than independence-implying p-values if overlap covariance cannot be obtained.

### C1 — Turn unresolved luminosity and generator coverage into an early feasibility deadline

**Classification: C — suggested strengthening.**

[A3] correctly refuses to import the wrong luminosity uncertainty, and [L2]
correctly refuses to omit generator tests (lines 46-63).  Add a Phase 2
decision date: either identify a dataset-appropriate luminosity uncertainty
and independent/varied generator inputs, or record why $\mu$ is conditional
and which shape/rate conclusions cannot be made.  This keeps the uncertainty
discussion honest without inventing a borrowed flat uncertainty.

### C2 — Make raw-object recovery preserve the already feasible angular cross-check

**Classification: C — suggested strengthening.**

Because current lepton four-vectors and event identifiers are retained, run
the angle calculator both before and after re-ntupling on matched events.
Document the fraction whose pairing/angles change after FSR or calibrated
reconstruction.  This protects the valid angular-NN programme from becoming
an accidental casualty of the required VBF/FSR reconstruction upgrade.

## Overall constructive path to a credible Phase 1 PASS

The artifact should be revised, not downscoped.  Its best features are already
the correct ones: a raw-jet VBF category, an angle-only NN with a cut-based
cross-check, independent flavour calibrations, DY as the nominal fake model,
and a shared-nuisance mass/rate fit.  A credible revision must add the
feasibility and validation contracts above, then create the missing commitment
ledger.  The strategy can then carry forward an honest statement of what is
known (candidate-level angular inputs and the supplied normalisations) and
what is still to be proven (raw-object VBF/FSR, fake transfer, calibrated fit
coverage, and 10 fb$^{-1}$ scale reach).

**PASS condition for this reviewer:** all A and B findings have explicit,
testable commitments in a revised strategy and `COMMITMENTS.md`; the revision
must preserve the user-required VBF, angular NN, DY nominal model, calibrated
simultaneous fit, reference window, and scale propagation rather than replacing
them with easier proxies.
