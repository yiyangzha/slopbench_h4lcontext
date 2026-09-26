"""Profile-likelihood scans of the H -> 4l fit in one or two parameters (stage 8).

    pixi run py -- analysis_v3/inference/scripts/scan_likelihood.py --model v5/model_v1 --yr4 v1 --label scan_mu_v1 \
        --x mu 0 2.5 26 [--y mH 123 127 21] [--poi-scheme fv] [--dimension 3D] [--no-refit] [--set b] \
        [--dataset data|asimov] [--fix-mh 125.09] [--stat-only] [--workers 8]

At every grid point the scanned parameters are fixed and every other parameter (m_H unless scanned or
fixed, the other POIs of the scheme, all nuisances) is profiled; q = 2 [NLL(point) - NLL(best)], with the
best fit the global one (and never above the lowest grid point).  --stat-only fixes every nuisance at its
best-fit value (the statistical curve).  1D: the 68 % and 95 % intervals are the crossings of q = 1 and
3.84 (linear interpolation); 2D: the contours q = 2.30 and 5.99.  The Asimov dataset is the one of
fit_model.py (mu = 1, m_H = 125 GeV or --fix-mh, nuisances at 0).  Writes
production_v3/inference/<model>/<label>/scan.json and plots/scan.{png,pdf}.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import multiprocessing
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

_STATE = {}


def make_likelihood(args) -> tuple:
    model_dir = PRODUCTION / "inference" / args.model
    model, events, signal_model, yr4 = fm.load(model_dir, args.yr4)
    syst = json.loads(args.systematics.read_text(encoding="utf-8"))
    lik = lk.Model(model, signal_model, yr4, args.dimension, not args.no_refit, args.set, syst, args.poi_scheme)
    if args.dataset == "asimov":
        ev = fm.asimov_events(lik, lik.default_parameters(1.0, args.fix_mh if args.fix_mh is not None else 125.0))
    else:
        ev = fm.data_events(model, events, lik.variable, lik.error)
    return lik, ev, syst


def init_worker(args):
    lik, ev, _ = make_likelihood(args)
    _STATE["fitter"] = fm.Fitter(lik, ev, args.fix_mh)


def one_point(task):
    fixed, start = task
    m = _STATE["fitter"].fit(start=start, fixed=fixed)
    return {"point": fixed, "nll": float(m.fval), "valid": bool(m.valid), "values": {n: float(m.values[n]) for n in _STATE["fitter"].names}}


def crossings(x: np.ndarray, q: np.ndarray, level: float) -> list:
    out = []
    for i in range(len(x) - 1):
        if (q[i] - level) * (q[i + 1] - level) < 0:
            out.append(float(x[i] + (level - q[i]) * (x[i + 1] - x[i]) / (q[i + 1] - q[i])))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--x", nargs=4, required=True, metavar=("NAME", "LOW", "HIGH", "N"))
    parser.add_argument("--y", nargs=4, default=None, metavar=("NAME", "LOW", "HIGH", "N"))
    parser.add_argument("--poi-scheme", default="inclusive")
    parser.add_argument("--dimension", default="3D", choices=["1D", "2D", "2Dmass", "3D"])
    parser.add_argument("--no-refit", action="store_true")
    parser.add_argument("--set", default="b", choices=["a", "b"])
    parser.add_argument("--dataset", default="data", choices=["data", "asimov"])
    parser.add_argument("--fix-mh", type=float, default=None)
    parser.add_argument("--stat-only", action="store_true")
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    out_dir = model_dir / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    lik, ev, syst = make_likelihood(args)
    fitter = fm.Fitter(lik, ev, args.fix_mh)
    best = fitter.fit()
    best_values = {n: float(best.values[n]) for n in fitter.names}
    axes = [(args.x[0], np.linspace(float(args.x[1]), float(args.x[2]), int(args.x[3])))]
    if args.y:
        axes.append((args.y[0], np.linspace(float(args.y[1]), float(args.y[2]), int(args.y[3]))))
    for name, _ in axes:
        if name not in fitter.names:
            raise SystemExit(f"{name} is not a parameter of the scheme {args.poi_scheme}: {fitter.names}")
    frozen = {n: best_values[n] for n in lik.nuisances} if args.stat_only else {}
    grid = [dict(frozen, **{axes[0][0]: float(a)}) for a in axes[0][1]] if len(axes) == 1 else \
        [dict(frozen, **{axes[0][0]: float(a), axes[1][0]: float(b)}) for b in axes[1][1] for a in axes[0][1]]
    tasks = [(point, best_values) for point in grid]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("fork"),
                                                initializer=init_worker, initargs=(args,)) as pool:
        rows = list(pool.map(one_point, tasks, chunksize=1))
    nll0 = min(float(best.fval), min(r["nll"] for r in rows))
    for r in rows:
        r["q"] = 2.0 * (r["nll"] - nll0)
        r["point"] = {k: v for k, v in r["point"].items() if k in (axes[0][0], axes[-1][0])}
    report = {"schema": "h4l_v3_scan/1", "model": str(model_dir), "dimension": args.dimension, "refit": not args.no_refit,
              "result_set": args.set, "dataset": args.dataset, "poi_scheme": args.poi_scheme, "fixed_mh": args.fix_mh,
              "stat_only": args.stat_only, "systematics": str(args.systematics), "best": {"nll": float(best.fval),
                                                                                           "valid": bool(best.valid),
                                                                                           "values": best_values},
              "axes": [{"name": n, "values": v.tolist()} for n, v in axes], "points": rows}
    out_dir.mkdir(parents=True)
    (out_dir / "plots").mkdir()
    fig, ax = plt.subplots(figsize=(6, 4.4))
    if len(axes) == 1:
        x = axes[0][1]
        q = np.array([r["q"] for r in rows])
        report["intervals"] = {"68": crossings(x, q, 1.0), "95": crossings(x, q, 3.84)}
        ax.plot(x, q, "k-", marker=".")
        for level in (1.0, 3.84):
            ax.axhline(level, color="grey", lw=0.6, ls=":")
        ax.set_xlabel(axes[0][0])
        ax.set_ylabel("-2 Delta ln L")
        ax.set_ylim(0, max(6.0, min(float(np.max(q)), 12.0)))
    else:
        x, y = axes[0][1], axes[1][1]
        q = np.array([r["q"] for r in rows]).reshape(len(y), len(x))
        contour = ax.contour(x, y, q, levels=[2.30, 5.99], colors=["k", "k"], linestyles=["-", "--"])
        ax.clabel(contour, fmt={2.30: "68%", 5.99: "95%"}, fontsize=7)
        ax.plot(best_values[axes[0][0]], best_values[axes[1][0]], "k+", ms=10, label="best fit")
        ax.plot(1.0, 1.0, "r*", ms=9, label="SM")
        ax.set_xlabel(axes[0][0])
        ax.set_ylabel(axes[1][0])
        ax.legend(fontsize=8)
        report["q_grid"] = q.tolist()
    ax.set_title(f"{args.dataset}, {args.dimension} {'refit' if not args.no_refit else 'no refit'}, set {args.set}"
                 f"{', stat only' if args.stat_only else ''}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "plots" / "scan.png", dpi=120)
    fig.savefig(out_dir / "plots" / "scan.pdf")
    plt.close(fig)
    (out_dir / "scan.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[scan] best {', '.join(f'{n} {best_values[n]:.4f}' for n, _ in axes)}; {len(rows)} points\n[scan] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
