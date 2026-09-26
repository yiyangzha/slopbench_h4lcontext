"""Saturated goodness-of-fit test of the H -> 4l likelihood with pseudo-experiments (stage 9).

    pixi run py -- analysis_v3/inference/scripts/gof.py --model v5/model_v1 --fit fit_v1 --yr4 v1 --label gof_v1 \
        --toys 300 --workers 8

The fitted unbinned model (the fit label's dimension, refit and result set) is binned per channel
(final state x category) in m4l (5 GeV bins in the window) x the D_bkg^kin template bins; the
expected content of a bin is the model's expected-count grid (fit_model.asimov_events at the
parameters) summed over the bin.  The saturated statistic
    q = 2 sum_bins [nu - n + n ln(n / nu)]
is computed for the data at the best fit of the data.  Pseudo-experiments are generated from the
best fit (fit_model.toy_events: Poisson counts, unbinned values; the global observables drawn from
N(theta_hat, 1), frequentist toys), each fitted with the same
likelihood (mu unbounded, m_H floating, nuisances profiled) and binned the same way; the p-value is
the fraction of pseudo-experiments with q at least the observed one.  Writes
production_v3/inference/<model>/<label>/gof.json and plots/gof.png.
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

M_STEP = 5.0
_STATE = {}


def bins_of(lik: lk.Model):
    m_edges = np.arange(lik.lo, lik.hi + 1e-9, M_STEP)
    return m_edges, lik.d_edges


def binned(lik: lk.Model, events: dict, weighted: bool) -> np.ndarray:
    """Counts per channel x m4l bin x D bin (flattened in the channel order of the likelihood)."""
    m_edges, d_edges = bins_of(lik)
    out = []
    for ch in lik.channels:
        ev = events.get(ch)
        if ev is None or len(ev["m"]) == 0:
            out.append(np.zeros((len(m_edges) - 1) * (len(d_edges) - 1)))
            continue
        w = ev["w"] if weighted else None
        h, _, _ = np.histogram2d(ev["m"], ev["d"], bins=[m_edges, d_edges], weights=w)
        out.append(h.ravel())
    return np.concatenate(out)


def expected(lik: lk.Model, p: np.ndarray) -> np.ndarray:
    return binned(lik, fm.asimov_events(lik, p), weighted=True)


def q_saturated(n: np.ndarray, nu: np.ndarray) -> float:
    nu = np.maximum(nu, 1e-12)
    terms = nu - n
    positive = n > 0
    terms[positive] += n[positive] * np.log(n[positive] / nu[positive])
    return float(2.0 * np.sum(terms))


def parameters_of(lik: lk.Model, values: dict) -> np.ndarray:
    p = np.zeros(len(lik.parameters))
    for name, value in values.items():
        if name in lik.index:
            p[lik.index[name]] = value
    return p


def init_worker(model_dir: str, yr4: str, dimension: str, refit: bool, result_set: str, syst: dict):
    model, events, signal_model, yr4_table = fm.load(Path(model_dir), yr4)
    _STATE["lik"] = lk.Model(model, signal_model, yr4_table, dimension, refit, result_set, syst)


def one_toy(args):
    p_best, seed, index = args
    lik = _STATE["lik"]
    rng = np.random.default_rng([seed, index])
    p_best = np.asarray(p_best)
    ev = fm.toy_events(lik, p_best, rng)
    # Frequentist pseudo-experiment: the global observables drawn around the fitted nuisances.
    lik.global_observables = np.array([p_best[lik.index[n]] for n in lik.nuisances]) + rng.normal(size=len(lik.nuisances))
    fitter = fm.Fitter(lik, ev)
    best = fitter.fit()
    p_toy = parameters_of(lik, {n: float(best.values[n]) for n in fitter.names})
    n = binned(lik, {ch: dict(v, w=np.ones(len(v["m"]))) for ch, v in ev.items()}, weighted=False)
    return {"index": index, "valid": bool(best.valid), "q": q_saturated(n, expected(lik, p_toy)), "mu": float(best.values["mu"]),
            "mH": float(best.values["mH"])}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--fit", required=True, help="label of the data fit whose best fit is tested")
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--toys", type=int, default=300)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    out_dir = model_dir / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    fit = json.loads((model_dir / args.fit / "fit.json").read_text(encoding="utf-8"))
    if fit.get("poi_scheme", "inclusive") != "inclusive" or fit.get("dataset") != "data":
        raise SystemExit("the goodness of fit tests the inclusive fit of the data")
    model, events, signal_model, yr4 = fm.load(model_dir, args.yr4)
    lik = lk.Model(model, signal_model, yr4, fit["dimension"], fit["refit"], fit["result_set"], fit["systematics"])
    p_best = parameters_of(lik, fit["fit"]["values"])
    data = fm.data_events(model, events, lik.variable, lik.error)
    n_obs = binned(lik, {ch: dict(v, w=np.ones(len(v["m"]))) for ch, v in data.items()}, weighted=False)
    nu_best = expected(lik, p_best)
    q_obs = q_saturated(n_obs, nu_best)
    print(f"[gof] observed q_saturated {q_obs:.2f} over {len(n_obs)} bins ({int(n_obs.sum())} events, expected {nu_best.sum():.1f})",
          flush=True)
    tasks = [(p_best.tolist(), args.seed, i) for i in range(args.toys)]
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("fork"),
                                                initializer=init_worker,
                                                initargs=(str(model_dir), args.yr4, fit["dimension"], fit["refit"], fit["result_set"],
                                                          fit["systematics"])) as pool:
        for k, r in enumerate(pool.map(one_toy, tasks, chunksize=2)):
            results.append(r)
            if (k + 1) % 25 == 0:
                print(f"[gof] {k + 1}/{len(tasks)} toys", flush=True)
    q_toys = np.array([r["q"] for r in results if r["valid"]])
    p_value = float(np.mean(q_toys >= q_obs)) if len(q_toys) else float("nan")
    p_error = float(np.sqrt(p_value * (1 - p_value) / len(q_toys))) if len(q_toys) else float("nan")
    out_dir.mkdir(parents=True)
    (out_dir / "plots").mkdir()
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.hist(q_toys, bins=40, histtype="step", color="#5790fc", label=f"pseudo-experiments ({len(q_toys)})")
    ax.axvline(q_obs, color="#e42536", label=f"data: q = {q_obs:.1f}, p = {p_value:.3f}")
    ax.set_xlabel("saturated q")
    ax.set_ylabel("pseudo-experiments")
    ax.set_title(f"Saturated GoF, {fit['dimension']} {'refit' if fit['refit'] else 'no refit'}, {len(n_obs)} bins", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "plots" / "gof.png", dpi=120)
    fig.savefig(out_dir / "plots" / "gof.pdf")
    plt.close(fig)
    report = {"schema": "h4l_v3_gof/1", "model": str(model_dir), "fit": str(model_dir / args.fit / "fit.json"),
              "method": f"saturated, binned per channel in m4l ({M_STEP:g} GeV) x D_bkg^kin template bins, toys from the best fit "
                        f"refitted with the full unbinned likelihood",
              "n_bins": int(len(n_obs)), "q_saturated": q_obs, "p_value": p_value, "p_value_error": p_error, "n_toys": len(q_toys),
              "n_toys_invalid": int(sum(not r["valid"] for r in results)), "seed": args.seed, "toys": results}
    (out_dir / "gof.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[gof] p-value {p_value:.3f} +- {p_error:.3f} from {len(q_toys)} toys\n[gof] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
