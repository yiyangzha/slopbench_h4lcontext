"""On-shell width of the Higgs boson from the m4l line shape (JHEP 11 (2017) 047 10.4).

    pixi run py -- analysis_v3/inference/scripts/fit_width.py --model v5/model_v1 --yr4 v1 --label width_v1

Unbinned fit of the 1D model (m4l with the Z1 refit, 105-140 GeV) with the resonant line shape the
per-final-state DCB convolved with a Breit-Wigner of width Gamma_H (h4l_shapes.dcb_bw_pdf; no
signal-background interference: substitution), m_H floating, the fermion- and vector-induced signal
strengths mu_F and mu_V free (the categories separate them), every nuisance profiled.  A profile
scan of Gamma_H in [0, 10] GeV gives -2 Delta ln L, observed and expected (Asimov at the SM: mu = 1,
m_H = 125 GeV, Gamma_H = 4.1 MeV); the 95 % CL upper limit is the crossing of 3.84 (the 68 % one of
1.0).  Writes production_v3/inference/<model>/<label>/width.json and plots/.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/inference"))
sys.path.insert(0, str(REPO / "analysis_v3/inference/scripts"))
import fit_model as fm  # noqa: E402
import h4l_likelihood as lk  # noqa: E402

GRID = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
SM_WIDTH = 0.0041


def crossing(x, y, level):
    """The first x above the minimum where y crosses the level (linear interpolation); None if it does not."""
    i0 = int(np.argmin(y))
    for i in range(i0, len(x) - 1):
        if y[i] <= level < y[i + 1]:
            return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return None


def scan(lik: lk.Model, events: dict) -> dict:
    fitter = fm.Fitter(lik, events)
    # The likelihood in Gamma_H can have several local minima: the best of fits started at 4.1 MeV, 0.5 and 1.5 GeV.
    best = min((fitter.fit(start={"GammaH": g}) for g in (SM_WIDTH, 0.5, 1.5)), key=lambda m: m.fval)
    values = {n: float(best.values[n]) for n in fitter.names}
    rows = []
    for g in GRID:
        # An invalid minimum is refitted from other starting points (m_H shifted, the strengths at 1); the lowest valid
        # minimum is kept (final review 2026-09-25: several observed points were flagged invalid).
        fits = [fitter.fit(start=dict(values, GammaH=g), fixed={"GammaH": g})]
        if not fits[0].valid:
            for dm in (-0.5, 0.5, -1.0, 1.0):
                fits.append(fitter.fit(start=dict(values, GammaH=g, mH=values["mH"] + dm, mu_F=1.0, mu_V=1.0), fixed={"GammaH": g}))
        valid = [f for f in fits if f.valid]
        m = min(valid or fits, key=lambda f: f.fval)
        rows.append({"GammaH": g, "nll": float(m.fval), "valid": bool(m.valid), "mH": float(m.values["mH"]), "attempts": len(fits)})
    nll0 = min(float(best.fval), min(r["nll"] for r in rows))
    for r in rows:
        r["q"] = 2.0 * (r["nll"] - nll0)
    x = np.array([r["GammaH"] for r in rows])
    y = np.array([r["q"] for r in rows])
    return {"best": {"GammaH": values["GammaH"], "mH": values["mH"], "mu_F": values["mu_F"], "mu_V": values["mu_V"], "nll": float(best.fval),
                     "valid": bool(best.valid)},
            "scan": rows, "limit_68": crossing(x, y, 1.0), "limit_95": crossing(x, y, 3.84)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--set", default="b", choices=["a", "b"])
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    out_dir = model_dir / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    model, events, signal_model, yr4 = fm.load(model_dir, args.yr4)
    syst = json.loads(args.systematics.read_text(encoding="utf-8"))
    lik = lk.Model(model, signal_model, yr4, "1D", True, args.set, syst, "fv", width=True)
    report = {"schema": "h4l_v3_width/1", "model": str(model_dir), "dimension": "1D", "refit": True, "result_set": args.set,
              "poi_scheme": "fv + GammaH", "systematics": syst, "grid": GRID,
              "substitutions": ["no signal-background interference term (no interference sample)"]}
    observed = scan(lik, fm.data_events(model, events, lik.variable, lik.error))
    p = lik.default_parameters(1.0, 125.0)
    p[lik.index["GammaH"]] = SM_WIDTH
    expected = scan(lik, fm.asimov_events(lik, p))
    report["observed"], report["expected"] = observed, expected
    out_dir.mkdir(parents=True)
    (out_dir / "plots").mkdir()
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for res, label, style in ((observed, "observed", "k-"), (expected, "expected", "b--")):
        ax.plot([r["GammaH"] for r in res["scan"]], [r["q"] for r in res["scan"]], style,
                label=f"{label}: Gamma_H < {res['limit_95']:.2f} GeV (95% CL)" if res["limit_95"] else label)
    for level in (1.0, 3.84):
        ax.axhline(level, color="grey", lw=0.6, ls=":")
    ax.set_xlabel("Gamma_H [GeV]")
    ax.set_ylabel("-2 Delta ln L")
    ax.set_ylim(0, 10)
    ax.set_title(f"On-shell width, 1D m4l (Z1 refit), {model['luminosity_fb']} fb-1", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "plots" / "width_scan.png", dpi=120)
    fig.savefig(out_dir / "plots" / "width_scan.pdf")
    plt.close(fig)
    (out_dir / "width.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[width] observed best Gamma_H {observed['best']['GammaH']:.3f} GeV, < {observed['limit_95']} GeV (95% CL); expected < "
          f"{expected['limit_95']} GeV\n[width] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
