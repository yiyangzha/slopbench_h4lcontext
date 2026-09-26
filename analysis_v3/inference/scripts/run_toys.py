"""Pseudo-experiment validation of the H -> 4l fit (stage 9).

    pixi run py -- analysis_v3/inference/scripts/run_toys.py --model v5/model_v1 --yr4 v1 --label toys_v1 \
        --mu 0 0.5 1 2 3 --mh 121 123 125 127 129 --ntoys 500 --workers 8

For every (mu, m_H) point, unbinned pseudo-experiments are generated from the model (fit_model.toy_events:
Poisson counts per component and channel, m4l from the component densities, D_bkg^kin and the relative
mass error from the templates; nuisances at 0, the global observables drawn from N(0, 1): frequentist
toys, so that the total intervals are tested) and fitted with the same likelihood as the data (mu
unbounded, m_H floating in 110-140 GeV, all nuisances profiled).  Per toy: the best fit, the MINOS
intervals of mu and m_H and, for mu_true = 0 and 1, the local significance.  The summary gives per point
the mean and width of the pulls (measured - true) / (interval on the side of the truth), the slopes of
measured against true mu and m_H, the fraction of toys with Z > 3 and Z > 5 at mu = 0 (false discovery)
and the coverage of the 68 % intervals.  Seeds: stream = hash(point) + toy index (reproducible).
Writes production_v3/inference/<model>/<label>/toys.json and plots.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import multiprocessing
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/inference"))
sys.path.insert(0, str(REPO / "analysis_v3/inference/scripts"))
import fit_model as fm  # noqa: E402
import h4l_likelihood as lk  # noqa: E402

_STATE = {}


def init_worker(model_path: str, yr4: str, dimension: str, refit: bool, syst_path: str):
    model_dir = Path(model_path)
    model, events, signal_model, yr4_table = fm.load(model_dir, yr4)
    syst = json.loads(Path(syst_path).read_text(encoding="utf-8"))
    _STATE["lik"] = lk.Model(model, signal_model, yr4_table, dimension, refit, "b", syst)


def one_toy(args):
    mu, mh, index, seed = args
    lik = _STATE["lik"]
    rng = np.random.default_rng([seed, index])
    p = np.zeros(len(lik.parameters))
    p[0], p[1] = mu, mh
    ev = fm.toy_events(lik, p, rng)
    # Frequentist pseudo-experiment: the global observables drawn around the generated nuisances (0).
    lik.global_observables = rng.normal(size=len(lik.nuisances))
    fitter = fm.Fitter(lik, ev)
    best = fitter.fit()
    out = {"mu_true": mu, "mh_true": mh, "index": index, "valid": bool(best.valid), "nll": float(best.fval),
           "mu": float(best.values["mu"]), "mH": float(best.values["mH"]), "n_events": int(sum(len(v["m"]) for v in ev.values()))}
    try:
        best.minos("mu", "mH")
        out["mu_err"] = [float(best.merrors["mu"].lower), float(best.merrors["mu"].upper)]
        out["mH_err"] = [float(best.merrors["mH"].lower), float(best.merrors["mH"].upper)]
    except Exception as error:  # noqa: BLE001
        out["minos_error"] = str(error)
        out["mu_err"] = [-float(best.errors["mu"]), float(best.errors["mu"])]
        out["mH_err"] = [-float(best.errors["mH"]), float(best.errors["mH"])]
    null = fitter.fit(start={n: float(best.values[n]) for n in fitter.names}, fixed={"mu": 0.0, "mH": float(best.values["mH"])})
    q0 = max(2.0 * (null.fval - best.fval), 0.0)
    out["Z"] = math.sqrt(q0) if out["mu"] > 0 else 0.0
    return out


def one_pair(args):
    """A paired probe (the benchmark's paired-injection test): dataset A at (mu, m_H), dataset B = A plus an
    independent signal-only pseudo-experiment of mu = 1 (every background event shared); the same global
    observables for both fits.  Returns both best fits and dmu = mu_B - mu_A (expected 1)."""
    mu, mh, index, seed = args
    lik = _STATE["lik"]
    rng = np.random.default_rng([seed, index])
    ev_a = fm.toy_events(lik, lik.default_parameters(mu, mh), rng)
    extra = fm.toy_events(lik, lik.default_parameters(1.0, mh), rng, components=("resonant", "nonres"))
    ev_b = {ch: {k: np.concatenate([ev_a[ch][k], extra[ch][k]]) for k in ev_a[ch]} for ch in ev_a}
    lik.global_observables = rng.normal(size=len(lik.nuisances))
    out = {"mu_true": mu, "mh_true": mh, "index": index}
    for tag, ev in (("a", ev_a), ("b", ev_b)):
        best = fm.Fitter(lik, ev).fit()
        out[f"mu_{tag}"], out[f"mH_{tag}"], out[f"valid_{tag}"] = float(best.values["mu"]), float(best.values["mH"]), bool(best.valid)
        out[f"mu_err_{tag}"] = float(best.errors["mu"])
    out["dmu"] = out["mu_b"] - out["mu_a"]
    return out


def pull(measured, true, err):
    lower, upper = err
    sigma = upper if measured < true else -lower
    return (measured - true) / sigma if sigma > 0 else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--mu", type=float, nargs="+", default=[1.0])
    parser.add_argument("--mh", type=float, nargs="+", default=[125.0])
    parser.add_argument("--ntoys", type=int, default=100)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--dimension", default="3D")
    parser.add_argument("--no-refit", action="store_true")
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    parser.add_argument("--paired", action="store_true", help="paired probes: A at mu, B = A + an extra mu = 1 signal")
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    out_dir = model_dir / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    tasks = []
    for mu in args.mu:
        for mh in args.mh:
            point_seed = args.seed + int(round(mu * 1000)) * 1000 + int(round(mh * 10))
            tasks += [(mu, mh, i, point_seed) for i in range(args.ntoys)]
    results = []
    # The fork start method: the default forkserver binds a Unix socket in TMPDIR, which the EOS FUSE mount refuses.
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("fork"),
                                                initializer=init_worker,
                                                initargs=(str(model_dir), args.yr4, args.dimension, not args.no_refit,
                                                          str(args.systematics))) as pool:
        for k, r in enumerate(pool.map(one_pair if args.paired else one_toy, tasks, chunksize=4)):
            results.append(r)
            if (k + 1) % 100 == 0:
                print(f"[toys] {k + 1}/{len(tasks)}", flush=True)
    if args.paired:
        summary = {}
        for mu in args.mu:
            for mh in args.mh:
                rows = [r for r in results if r["mu_true"] == mu and r["mh_true"] == mh and r["valid_a"] and r["valid_b"]]
                d = np.array([r["dmu"] for r in rows])
                summary[f"mu{mu:g}_mh{mh:g}"] = {"mu_true": mu, "mh_true": mh, "n_valid": len(rows), "dmu_mean": float(np.mean(d)),
                                                 "dmu_median": float(np.median(d)), "dmu_width": float(np.std(d, ddof=1)),
                                                 "dmu_mean_error": float(np.std(d, ddof=1) / math.sqrt(len(d)))}
        out_dir.mkdir(parents=True)
        report = {"schema": "h4l_v3_paired_toys/1", "model": str(model_dir), "dimension": args.dimension, "refit": not args.no_refit,
                  "seed": args.seed, "ntoys": args.ntoys, "summary": summary, "toys": results}
        (out_dir / "toys.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
        for key, s in summary.items():
            print(f"[paired] {key}: dmu mean {s['dmu_mean']:.3f} +- {s['dmu_mean_error']:.3f}, median {s['dmu_median']:.3f}, "
                  f"width {s['dmu_width']:.3f}", flush=True)
        print(f"[paired] {out_dir}")
        return 0
    summary = {}
    for mu in args.mu:
        for mh in args.mh:
            rows = [r for r in results if r["mu_true"] == mu and r["mh_true"] == mh and r["valid"]]
            pm = np.array([pull(r["mu"], mu, r["mu_err"]) for r in rows])
            ph = np.array([pull(r["mH"], mh, r["mH_err"]) for r in rows])
            z = np.array([r["Z"] for r in rows])
            key = f"mu{mu:g}_mh{mh:g}"
            summary[key] = {"mu_true": mu, "mh_true": mh, "n_valid": len(rows), "n_total": args.ntoys,
                            "mu_mean": float(np.mean([r["mu"] for r in rows])), "mu_median": float(np.median([r["mu"] for r in rows])),
                            "mH_mean": float(np.mean([r["mH"] for r in rows])), "mH_median": float(np.median([r["mH"] for r in rows])),
                            "pull_mu_mean": float(np.nanmean(pm)), "pull_mu_width": float(np.nanstd(pm, ddof=1)),
                            "pull_mH_mean": float(np.nanmean(ph)), "pull_mH_width": float(np.nanstd(ph, ddof=1)),
                            "coverage_mu": float(np.nanmean(np.abs(pm) <= 1)), "coverage_mH": float(np.nanmean(np.abs(ph) <= 1)),
                            "frac_Z_gt_3": float(np.mean(z > 3)), "frac_Z_gt_5": float(np.mean(z > 5)), "median_Z": float(np.median(z))}
    # Slopes of measured against true (medians), over the scanned points.
    slopes = {}
    if len(args.mu) > 1:
        for mh in args.mh:
            x = np.array(args.mu)
            y = np.array([summary[f"mu{m:g}_mh{mh:g}"]["mu_median"] for m in args.mu])
            slopes[f"mu_at_mh{mh:g}"] = float(np.polyfit(x, y, 1)[0])
    if len(args.mh) > 1:
        for mu in args.mu:
            if mu <= 0:
                continue
            x = np.array(args.mh)
            y = np.array([summary[f"mu{mu:g}_mh{m:g}"]["mH_median"] for m in args.mh])
            slopes[f"mH_at_mu{mu:g}"] = float(np.polyfit(x, y, 1)[0])
    out_dir.mkdir(parents=True)
    report = {"schema": "h4l_v3_toys/1", "model": str(model_dir), "dimension": args.dimension, "refit": not args.no_refit,
              "seed": args.seed, "ntoys": args.ntoys, "summary": summary, "slopes": slopes, "toys": results}
    (out_dir / "toys.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    for key, s in summary.items():
        print(f"[toys] {key}: mu median {s['mu_median']:.3f} pull {s['pull_mu_mean']:+.3f} +- width {s['pull_mu_width']:.3f}; "
              f"mH median {s['mH_median']:.3f} pull {s['pull_mH_mean']:+.3f} width {s['pull_mH_width']:.3f}; Z>3 {s['frac_Z_gt_3']:.3f}",
              flush=True)
    print(f"[toys] slopes {slopes}\n[toys] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
