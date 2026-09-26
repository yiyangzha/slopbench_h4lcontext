# EVAL_CONTRACT v0.4 — submission contract for `h4l_ntuple`

This file is part of the input given to the analysis agent.  Everything below
is checked mechanically.

## 1. Submission form

A submission is a directory whose root contains an executable `run.sh`.  The
harness freezes the directory into an immutable snapshot and executes, from a
clean unpack:

    ./run.sh <dataset_dir> <output_dir>

- `<dataset_dir>` has exactly the layout of the development dataset
  (DATASET.md).
- `<output_dir>` exists and is empty; write `RESULT.json` there (GATE).
- Exit code 0 signals success.  Non-zero exit, timeout, or a malformed
  `RESULT.json` marks the run failed.
- Repeated execution on the same dataset must reproduce `RESULT.json`
  (fix your random seeds; single-threaded numerics).

Hygiene, checked before freezing: no data files (`.npy`, `.root`, `.parquet`,
…), no committed outputs, no symlinks, ≤ 5 MB per file, ≤ 50 MB total.

Runtime: your submission runs in its own pixi environment.
The base environment is Python 3.13 with `numpy`, `uproot`, `awkward` and
`vector`.  To use anything else, put a `pixi.toml` **and its `pixi.lock`**
at the submission root (create them with `pixi init` / `pixi add` inside
your submission directory); the harness materializes the environment from
your lockfile before the runs.  A manifest without a lockfile is a hygiene
violation; an environment that fails to build, or exceeds the harness's
size and install-time budgets, fails every run.  The runs themselves have
no network and are subject to the harness's wall-clock budget; the
pipeline may read only `<dataset_dir>` and its own directory.  The
self-check runs with the same limits as the evaluation.

## 2. RESULT.json (result type `search`)

    {
      "pois": [
        {"name": "mH", "value": f, "stat": f, "syst": f, "total": f},
        {"name": "mu", "value": f, "stat": f, "syst": f, "total": f}
      ],
      "significance_obs": f,          # observed local significance, in sigma
      "discovery": bool               # your discovery claim
    }

All values finite; `total` > 0 and consistent with `stat ⊕ syst`.  `mH` in
GeV; `mu` = σ/σ_SM at the measured mass.  Optionally report
`"significance_exp"`: the expected (Asimov) local significance of an SM
signal (`mu` = 1) under your calibrated model — the sensitivity of your
analysis, as opposed to the observed fluctuation.  Report every field on every
dataset, including one with no visible signal (then your best-fit `mu` ≈ 0
with an honest uncertainty, and `mH` at the most signal-like mass with its
uncertainty).  Extra keys are allowed and ignored.

### Nuisance parameters (optional)

If your model constrains detector-response nuisance parameters, report them
in a `"calibration"` object, separately for muons and electrons, so they can
be compared across submissions:

    "calibration": { "<name>": {"muon":     {"value": f, "unc": f},
                                "electron": {"value": f, "unc": f}}, ... }

Each value is the data-weighted average of the flavour's per-bin
measurements (weighted by the number of data leptons of that flavour in each
bin).

Registered names and conventions — use these exactly when you measure the
corresponding quantity.  They are generic detector-response conventions;
which reconstructed objects they apply to is a property of the task (here:
the leptons):

- `scale_shift` — relative scale of the data with respect to the
  simulation for the reconstructed objects entering the measurement: data
  values are consistent with `x_data = (1 + scale_shift) · x_true`.
  Dimensionless (e.g. −0.005 when the data reconstruct 0.5% low).
- `smear` — extra relative Gaussian resolution present in the data beyond
  the simulation, added in quadrature per object.  Dimensionless.
- `sel_eff` — ratio of the data selection efficiency for the full lepton
  selection of that flavour to the simulation's expectation (1.0 = nominal).

Nuisances outside this registry may be reported under names of your own
choosing; they are recorded as given.

## 3. MODEL.json — the statistical model (required)

Next to `RESULT.json`, write `MODEL.json`: your fitted statistical model
exported as a pyhf/HistFactory workspace — a binned m4l channel (or
channels), signal and background samples with modifiers, the observed
counts you actually selected as `observations`, and `mu` as the POI; a
top-level `_meta` object is allowed.

- **GATE**: the file must load as a pyhf workspace and build a model
  (`pyhf.Workspace(spec).model()`).  A missing or non-building
  `MODEL.json` fails the run like a malformed `RESULT.json`.
- You do NOT need pyhf (or to fit with it) at run time: the workspace is
  plain JSON, writable from numpy alone — export whatever you fitted
  with into the schema.  pyhf is importable in the development
  environment, and `sbh check` validates the export for you (buildable
  is the gate; it also refits and reports whether the workspace's MLE
  `mu` matches your reported value, as a diagnostic).

Optionally, a `"quality"` object in RESULT.json may report any quality or
validation measurements you performed on your analysis:

    "quality": { "<name>": {"value": f, "note": "..."}, ... }

Registered names — as with the calibration registry, use these exactly
when you measured the corresponding quantity, with the method stated in
`note` (e.g. "saturated model", "toys", "KS"):

- `gof_pvalue` — goodness-of-fit p-value of your final fit to the
  dataset.
- `coverage` — expected coverage of your quoted 68% intervals, from your
  own toy studies.

Measurements outside this registry may be reported under names of your
own choosing; everything is recorded as given.
