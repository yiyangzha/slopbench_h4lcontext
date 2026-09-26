"""Fits of the binned likelihood: best fit with m_H floating (grid profile then free fit), MINOS intervals, statistical
intervals with the nuisances fixed at their best fit, fixed-nuisance impacts for the systematic part, local
significances (observed and post-fit Asimov with mu = 1), the saturated goodness of fit (the main analysis's gof.py:
q = 2 sum [nu - n + n ln(n / nu)] over 5 GeV m4l bins x the template bins of every channel, at the best fit; toys from
the best fit with the global observables drawn, each refitted with m_H floating) and the coverage of the mu interval."""

from __future__ import annotations

import math

import numpy as np

from . import config as C
from .model import NUISANCES, PARAMS, rebin_gof, saturated_q

GOF_GROUP = int(round(C.GOF_M_STEP / C.BIN_WIDTH))


def interval(m, name):
    me = m.merrors[name] if name in m.merrors else None
    if me is not None and me.is_valid:
        return [float(me.lower), float(me.upper)]
    e = float(m.errors[name])
    return [-e, e]


_TOY = {}


def _toy(i):
    """One pseudo-experiment (its own seeded generator, so the result does not depend on the process that runs it)."""
    lik, fitter, p_hat, nu_hat, seed = _TOY["lik"], _TOY["fitter"], _TOY["p_hat"], _TOY["nu_hat"], _TOY["seed"]
    rng = np.random.default_rng([seed, i])
    obs = [rng.poisson(nu).astype(float) for nu in nu_hat]
    aux = p_hat[2:] + rng.normal(size=len(NUISANCES))
    m = fitter.fit(p_hat, observed=obs, aux=aux)
    if not m.valid:
        return None
    q = saturated_q(rebin_gof(lik.expected(np.array(m.values)), GOF_GROUP), rebin_gof(obs, GOF_GROUP))
    try:
        m.minos("mu")
        lo, hi = interval(m, "mu")
    except Exception:  # noqa: BLE001
        lo, hi = -float(m.errors["mu"]), float(m.errors["mu"])
    mu_t = float(m.values["mu"])
    return {"q": q, "mu": mu_t, "covered": bool(mu_t + lo <= p_hat[0] <= mu_t + hi)}


def run_inference(lik, fitter, log, n_toys=200, seed=1) -> dict:
    best, scan = fitter.global_fit()
    p_hat = np.array(best.values)
    try:
        best.minos("mu", "mH")
    except Exception as error:  # noqa: BLE001
        log(f"[fit] MINOS failed: {error}")
    total = {n: interval(best, n) for n in ("mu", "mH")}
    fixed = {n: float(best.values[n]) for n in NUISANCES}
    stat_fit = fitter.fit(p_hat, fixed=fixed)
    try:
        stat_fit.minos("mu", "mH")
    except Exception as error:  # noqa: BLE001
        log(f"[fit] stat MINOS failed: {error}")
    stat = {n: interval(stat_fit, n) for n in ("mu", "mH")}
    impacts = {}
    for n in NUISANCES:
        err = float(best.errors[n])
        shifts = {}
        for sign in (1, -1):
            m = fitter.fit(p_hat, fixed={n: float(best.values[n]) + sign * err})
            shifts[sign] = {poi: float(m.values[poi]) - float(best.values[poi]) for poi in ("mu", "mH")}
        impacts[n] = {"up": shifts[1], "down": shifts[-1], "post_fit_error": err, "pull": float(best.values[n])}
    syst = {poi: math.sqrt(sum((0.5 * (abs(v["up"][poi]) + abs(v["down"][poi]))) ** 2 for v in impacts.values())) for poi in ("mu", "mH")}
    null = fitter.fit(p_hat, fixed={"mu": 0.0, "mH": float(best.values["mH"])})
    q0 = max(2.0 * (null.fval - best.fval), 0.0)
    z_obs = math.sqrt(q0) if best.values["mu"] > 0 else 0.0
    # Expected significance: the post-fit Asimov dataset with mu = 1 at the fitted m_H (global observables at the fitted
    # nuisances, so that the Asimov fit returns them).
    p_exp = p_hat.copy()
    p_exp[0] = 1.0
    asimov = lik.expected(p_exp)
    aux_exp = p_hat[2:].copy()
    fa = fitter.fit(p_exp, fixed={"mH": float(p_hat[1])}, observed=asimov, aux=aux_exp)
    fa0 = fitter.fit(p_exp, fixed={"mu": 0.0, "mH": float(p_hat[1])}, observed=asimov, aux=aux_exp)
    z_exp = math.sqrt(max(2.0 * (fa0.fval - fa.fval), 0.0))
    log(f"[fit] mu {best.values['mu']:.4f} [{total['mu'][0]:+.4f}, {total['mu'][1]:+.4f}] (stat [{stat['mu'][0]:+.4f}, "
        f"{stat['mu'][1]:+.4f}], syst {syst['mu']:.4f}); mH {best.values['mH']:.3f} [{total['mH'][0]:+.3f}, {total['mH'][1]:+.3f}] "
        f"(syst {syst['mH']:.3f}); Z {z_obs:.2f} (expected {z_exp:.2f}); valid {best.valid}")
    # Toys from the best fit: Poisson counts and global observables drawn around the fitted nuisances; each toy fitted from
    # the generating point (m_H free).
    from .parallel import pmap
    nu_hat = lik.expected(p_hat)
    q_obs = saturated_q(rebin_gof(nu_hat, GOF_GROUP), rebin_gof(lik.observed, GOF_GROUP))
    _TOY.update({"lik": lik, "fitter": fitter, "p_hat": p_hat, "nu_hat": nu_hat, "seed": seed})
    toys = [t for t in pmap(_toy, list(range(n_toys)), chunksize=4) if t is not None]
    valid = len(toys)
    q_toys = np.array([t["q"] for t in toys])
    covered = [t["covered"] for t in toys]
    mu_toys = [t["mu"] for t in toys]
    p_gof = float(np.mean(q_toys >= q_obs)) if len(q_toys) else float("nan")
    coverage = float(np.mean(covered)) if covered else float("nan")
    log(f"[gof] saturated q {q_obs:.1f}; toys {valid}/{n_toys} valid: p = {p_gof:.3f}; coverage of the mu interval {coverage:.3f}")
    n_bins = int(len(rebin_gof(lik.observed, GOF_GROUP)))
    return {"best": {n: float(best.values[n]) for n in PARAMS}, "errors": {n: float(best.errors[n]) for n in PARAMS},
            "valid": bool(best.valid), "nll": float(best.fval), "total": total, "stat": stat, "syst": syst, "impacts": impacts,
            "q0": q0, "z_obs": z_obs, "z_exp": z_exp, "scan": scan,
            "gof": {"q_saturated": q_obs, "p_value": p_gof, "n_toys": int(n_toys), "n_valid": valid, "n_bins": n_bins,
                    "q_toys_mean": float(np.mean(q_toys)) if len(q_toys) else None},
            "coverage": {"value": coverage, "n": len(covered), "mu_generated": float(p_hat[0]),
                         "mu_toys_mean": float(np.mean(mu_toys)) if mu_toys else None,
                         "mu_toys_std": float(np.std(mu_toys)) if mu_toys else None}}
