"""Fits of the H -> 4l likelihood (stage 7): mu and m_H with profile-likelihood intervals.

    pixi run py -- analysis_v3/inference/scripts/fit_model.py --model v5/model_v1 --yr4 v1 --label fit_v1 \
        [--dimension 3D] [--no-refit] [--set b] [--dataset data|asimov] [--toys N] \
        [--poi-scheme inclusive|final_state|category|fv|mode|stxs0|fid_fs|fid_int|fid_<obs>] [--fix-mh 125.09] \
        [--fiducial production_v3/signal_model/v5/gen_v1/fiducial_xsec.json]

The likelihood is analysis_v3/inference/h4l_likelihood.py.  mu is unbounded; m_H floats in
[110, 140] GeV (or is fixed with --fix-mh).  With a POI scheme other than inclusive the global mu
is fixed at 1 and the scheme's POIs (signal strengths or cross-section ratios, in [0, 20] as in the
paper: user decision 2026-09-25) are fitted.
Reported per fit:
  * the best-fit values and the MINOS (profile-likelihood) intervals of mu and m_H (total);
  * the statistical intervals: MINOS of the same likelihood with every nuisance fixed at its
    best-fit value;
  * the systematic part from explicit fixed-nuisance impacts: for every nuisance, the shift of
    the POI when it is fixed at its best fit +- its post-fit error and the rest re-minimized;
    syst = quadrature sum of the symmetrized impacts (never a subtraction of intervals);
  * the local significance with m_H floating: sqrt(2 [NLL(mu = 0) - NLL(best)]);
  * the pulls of the nuisances.
The Asimov dataset is the expected-count grid of the model at mu = 1, m_H = 125 GeV (nuisances at 0):
m4l cells of 0.25 GeV x the D and relative-error template bins, as weighted events (the pre-fit
expectation); --dataset asimov_postfit --postfit-from <fit> builds it at the best fit of that fit (the
post-fit expectation of the paper).
Toys: unbinned pseudo-experiments generated from the model (Poisson counts per component, m4l from
the component density, D and the relative error from the templates) with the global observables
drawn from N(0, 1) (frequentist toys), fitted like the data.
Writes production_v3/inference/<model>/<label>/fit.json.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from iminuit import Minuit

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/inference"))
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_likelihood as lk  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

MH_LIMITS = (110.0, 140.0)
E_SUB = 5  # D_mass points per template bin of the Asimov grid


def load(model_dir: Path, yr4_version: str):
    model = json.loads((model_dir / "model.json").read_text(encoding="utf-8"))
    events = dict(np.load(model_dir / "events.npz"))
    signal_model = json.loads(Path(model["signal_model"]).read_text(encoding="utf-8"))
    yr4 = json.loads((PRODUCTION / "signal_model" / "yr4" / yr4_version / "yr4.json").read_text(encoding="utf-8"))
    return model, events, signal_model, yr4


def data_events(model: dict, events: dict, variable: str, error: str) -> dict:
    out = {}
    for ch in model["channels"]:
        out[ch] = {"m": events[f"{ch}__{variable}"], "dm": events[f"{ch}__{error}"], "d": events[f"{ch}__d_bkg_kin"]}
    return out


def asimov_events(lik: lk.Model, p: np.ndarray) -> dict:
    """Expected-count grid as weighted events."""
    out = {}
    step = 0.25
    m = np.arange(lik.lo + step / 2, lik.hi, step)
    d_rep = 0.5 * (lik.d_edges[1:] + lik.d_edges[:-1])
    for ch in lik.channels:
        fs = lik.model["channels"][ch]["final_state"]
        fse = lik.model["final_states"][fs]
        # D_mass: E_SUB points per template bin, evenly inside the bin's raw quantile range (the range the events
        # occupy, as in the toys), each carrying 1/E_SUB of the bin probability (final review 2026-09-25: the bin
        # midpoints, e1/2 for the first bin, made the per-event widths too small).
        raw = np.array(fse["e_edges_raw"][lik.variable])
        e_points = [(raw[i] + (raw[i + 1] - raw[i]) * (np.arange(E_SUB) + 0.5) / E_SUB) for i in range(len(raw) - 1)]
        e_rep = np.concatenate(e_points)
        y = lik.yields(p, ch)
        ms, ds, es, ws = [], [], [], []
        for di, dval in enumerate(d_rep if lik.use_d else d_rep[:1]):
            for ei, eval_ in enumerate(e_rep if lik.use_e else [np.nan]):
                mm = m
                dd = np.full(len(m), dval)
                ee = np.full(len(m), eval_ if lik.use_e else 0.01)
                dm = ee * mm
                total = np.zeros(len(m))
                for mode in lk.MODES:
                    dens = lik.signal_m_density(p, fs, mm, dm, mode)
                    if lik.use_d:
                        dens = dens * lik.d_density("signal", fs, mm, dd)
                    if lik.use_e:
                        dens = dens * lik.e_density("signal", fs, ee)
                    total += y[mode] * dens
                for name in ("qqZZ", "ggZZ", "zx"):
                    if name == "zx":
                        dens = shapes.zx_pdf(mm, fse["zx"]["shape"], lik.lo, lik.hi)
                    else:
                        dens = shapes.bernstein_pdf(mm, fse["backgrounds"][name][f"bernstein_{lik.variable}"], lik.lo, lik.hi)
                    if lik.use_d:
                        dens = dens * lik.d_density(name, fs, mm, dd)
                    if lik.use_e:
                        dens = dens * lik.e_density(name, fs, ee)
                    total += y[name] * dens
                # Without D the D loop has one pass (the D density absent); each D_mass sub-point carries 1/E_SUB of
                # its bin probability.
                ms.append(mm)
                ds.append(dd)
                es.append(dm)
                ws.append(total * step / (E_SUB if lik.use_e else 1.0))
        out[ch] = {"m": np.concatenate(ms), "d": np.concatenate(ds), "dm": np.concatenate(es), "w": np.concatenate(ws)}
    return out


def standard_dcb_sampler(v: dict):
    """Inverse-CDF sampler of the standard DCB (tails in units of sigma) on a wide t grid."""
    t = np.linspace(-60.0, 60.0, 240001)
    cdf = shapes.dcb_primitive(t, v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"])
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    return lambda u: np.interp(u, cdf, t)


def toy_events(lik: lk.Model, p: np.ndarray, rng: np.random.Generator,
               components: tuple = ("resonant", "nonres", "qqZZ", "ggZZ", "zx")) -> dict:
    """Unbinned pseudo-experiment from the model: Poisson counts per component and channel; the relative mass
    error from its template (uniform within the raw quantile bins); signal m4l from the DCB with the event's own
    width (t from the standard DCB, m = mean + sigma_i t, events outside the window regenerated), the VH
    non-resonant part and the backgrounds by inverse CDF of their densities; D_bkg^kin from P(D | m4l bin)."""
    out = {}
    grid = np.linspace(lik.lo, lik.hi, 3501)
    centres = 0.5 * (grid[1:] + grid[:-1])
    th = {n: p[lik.index[n]] for n in ("scale_mu", "scale_e", "res_mu", "res_e")}
    d_scale = lik.syst.get("scale_uncertainty", {"mu": 0.0005, "e": 0.0015})
    d_res = lik.syst.get("resolution_uncertainty", {"mu": 0.10, "e": 0.10})
    for ch in lik.channels:
        fs = lik.model["channels"][ch]["final_state"]
        fse = lik.model["final_states"][fs]
        mh = lik.mh_of(p, fs)
        k = mh / 125.0
        y = lik.yields(p, ch)
        raw_edges = np.array(fse["e_edges_raw"][lik.variable])
        sm = fse["signal"][lik.variable]
        vh = sm["vh_nonresonant"]
        f_mu, f_e = lik.flavour_fraction[fs]
        g_mu, g_e = lik.resolution_fraction[fs]
        scale = 1.0 + f_mu * d_scale["mu"] * th["scale_mu"] + f_e * d_scale["e"] * th["scale_e"]
        res = max(1.0 + g_mu * d_res["mu"] * th["res_mu"] + g_e * d_res["e"] * th["res_e"], 0.2)
        v = sm["dcb_per_event_width"]["values"] if lik.use_e else sm["dcb"]["values"]
        sampler = standard_dcb_sampler(v)
        mean = v["mean"] * k * scale
        ms, ds, dms = [], [], []

        def sample_e(key, n):
            probs = np.array(fse["e_templates"][key][lik.variable])
            idx = rng.choice(len(probs), size=n, p=probs / probs.sum())
            return rng.uniform(raw_edges[idx], raw_edges[idx + 1])

        def inverse_cdf(dens, n):
            cdf = np.cumsum(dens)
            cdf /= cdf[-1]
            return np.interp(rng.uniform(size=n), cdf, centres)

        n_res = rng.poisson(max(y["ggH"] + y["VBF"] + vh["f_res"] * y["VH"], 0.0))
        n_nonres = rng.poisson(max((1 - vh["f_res"]) * y["VH"], 0.0))
        generated = [("resonant", n_res, "signal"), ("nonres", n_nonres, "signal")] + \
                    [(name, rng.poisson(max(y[name], 0.0)), name) for name in ("qqZZ", "ggZZ", "zx")]
        for comp, n, key in generated:
            if n == 0 or comp not in components:
                continue
            if lik.use_e:
                e = sample_e(key, n)
            else:
                e = np.full(n, 0.01)
            if comp == "resonant":
                m = np.full(n, np.nan)
                todo = np.arange(n)
                while len(todo):
                    # In the 3D model the width is s x D_mass x m_H at the event's relative error D_mass = e.
                    width = v["width"] * (e[todo] * mh if lik.use_e else k) * res
                    cand = mean + width * sampler(rng.uniform(size=len(todo)))
                    ok = (cand > lik.lo) & (cand < lik.hi)
                    m[todo[ok]] = cand[ok]
                    todo = todo[~ok]
            elif comp == "nonres":
                m = inverse_cdf(shapes.landau_pdf(centres, vh["landau_mpv"], vh["landau_width"], lik.lo, lik.hi), n)
            elif comp == "zx":
                m = inverse_cdf(shapes.zx_pdf(centres, fse["zx"]["shape"], lik.lo, lik.hi), n)
            else:
                m = inverse_cdf(shapes.bernstein_pdf(centres, fse["backgrounds"][comp][f"bernstein_{lik.variable}"], lik.lo, lik.hi), n)
            table = np.array(fse["d_templates"][key])
            mi = np.clip(np.searchsorted(lik.m_edges, m, side="right") - 1, 0, len(lik.m_edges) - 2)
            cum = np.cumsum(table[mi], axis=1)
            cum /= cum[:, -1:]
            kbin = (rng.uniform(size=(n, 1)) > cum).sum(axis=1)
            lo_d = np.maximum(lik.d_edges[kbin], 0.0)
            hi_d = np.minimum(lik.d_edges[kbin + 1], 1.0)
            d = rng.uniform(lo_d, np.maximum(hi_d, lo_d))
            ms.append(m)
            ds.append(d)
            dms.append(e * m)
        if ms:
            out[ch] = {"m": np.concatenate(ms), "d": np.concatenate(ds), "dm": np.concatenate(dms)}
        else:
            out[ch] = {"m": np.zeros(0), "d": np.zeros(0), "dm": np.zeros(0)}
    return out


class Fitter:
    def __init__(self, lik: lk.Model, events: dict, fixed_mh: float | None = None):
        self.lik = lik
        self.prepared = lik.prepare(events)
        self.names = lik.parameters
        self.fixed_mh = fixed_mh

    def minuit(self, start: dict | None = None, fixed: dict | None = None) -> Minuit:
        x0 = np.zeros(len(self.names))
        x0[0], x0[1] = 1.0, 125.0
        if start:
            for k, v in start.items():
                x0[self.names.index(k)] = v
        for k in self.lik.extra_pois:
            if not start or k not in start:
                x0[self.names.index(k)] = 125.0 if k.startswith("mH_") else 1.0
        m = Minuit(lambda x: self.lik.nll(np.asarray(x), self.prepared), x0, name=self.names)
        m.errordef = Minuit.LIKELIHOOD
        m.limits["mH"] = MH_LIMITS
        # mu unbounded in practice: limits far outside any reachable value (a negative signal is stopped by the density
        # at the events).
        m.limits["mu"] = (-20.0, 50.0)
        for n in self.lik.extra_pois:
            # The scheme's signal strengths and cross-section ratios are non-negative, as in the paper (user decision
            # 2026-09-25); the inclusive mu stays unbounded.
            m.limits[n] = MH_LIMITS if n.startswith("mH_") else (0.0, 20.0)
        if "GammaH" in self.names:
            # The width is physical only at or above zero; start at the SM value (4.1 MeV).
            m.limits["GammaH"] = (0.0, 20.0)
            if not start or "GammaH" not in start:
                m.values["GammaH"] = 0.0041
        for n in self.lik.nuisances:
            m.limits[n] = (-6.0, 6.0)
        if [n for n in self.lik.extra_pois if n != "GammaH" and not n.startswith("mH_")]:
            # The scheme's POIs carry the signal normalization: the global mu is fixed at 1.
            m.values["mu"] = 1.0
            m.fixed["mu"] = True
        if self.fixed_mh is not None:
            m.values["mH"] = self.fixed_mh
            m.fixed["mH"] = True
        if fixed:
            for k, v in fixed.items():
                m.values[k] = v
                m.fixed[k] = True
        m.strategy = 1
        return m

    def fit(self, start=None, fixed=None) -> Minuit:
        m = self.minuit(start, fixed)
        m.migrad(ncall=20000)
        if not m.valid:
            m.migrad(ncall=20000)
        if not m.valid:
            # A minimum above the EDM goal (many POIs at bounds, set (a)): simplex, then Migrad with strategy 2.
            m.simplex(ncall=20000)
            m.strategy = 2
            m.migrad(ncall=40000)
        # A minimum on the boundary of the physical region (a strength at 0 in a category without signal-like
        # events): the internal transformation of a bounded parameter at its limit makes the Hessian there singular
        # (HESSE then raises the EDM above the goal) and MINOS unusable.  The constrained minimum: every parameter at
        # a limit where the NLL rises into the allowed region is fixed at that limit and the others are minimized;
        # the one-sided interval of such a parameter comes from a profile scan (run_fit).
        self.bound = self.at_limit(m)
        if self.bound:
            for n, edge in self.bound:
                m.values[n] = edge
                m.fixed[n] = True
            m.strategy = 1
            m.migrad(ncall=40000)
        m.hesse()
        return m

    def at_limit(self, m: Minuit) -> list:
        """The free parameters at a limit (within 1e-4 of their range) whose NLL rises when moved inside."""
        out = []
        x = np.array(m.values)
        f0 = self.lik.nll(x, self.prepared)
        for i, n in enumerate(self.names):
            if m.fixed[n] or m.limits[n] is None:
                continue
            lo, hi = m.limits[n]
            if not (math.isfinite(lo) and math.isfinite(hi)):
                continue
            tol = 1e-4 * (hi - lo)
            for edge, inward in ((lo, 1.0), (hi, -1.0)):
                if abs(x[i] - edge) < tol:
                    y = x.copy()
                    y[i] = edge
                    at_edge = self.lik.nll(y, self.prepared)
                    y[i] = edge + inward * 10.0 * tol
                    if self.lik.nll(y, self.prepared) > min(at_edge, f0):
                        out.append((n, edge))
        return out


def one_sided_interval(fitter: Fitter, values: dict, poi: str, edge: float, nll_min: float, fixed: dict) -> list:
    """[lower, upper] offsets of a POI at its limit `edge`: zero on the limit side, the profile edge (NLL + 0.5) on the
    other, found by bracketing and bisection of the profiled NLL (the other parameters minimized, `fixed` held)."""
    lo, hi = fitter.minuit().limits[poi]
    inward = 1.0 if edge == lo else -1.0
    far = hi if inward > 0 else lo

    def excess(x: float) -> float:
        return fitter.fit(start=values, fixed={**fixed, poi: x}).fval - nll_min - 0.5

    inside, step = edge, 0.25
    outside = edge + inward * step
    while excess(outside) < 0.0:
        inside = outside
        step *= 2.0
        if abs(outside - far) < 1e-9:
            return [0.0, far - edge] if inward > 0 else [far - edge, 0.0]
        outside = edge + inward * step
        if (outside - far) * inward > 0:
            outside = far
    for _ in range(30):
        mid = 0.5 * (inside + outside)
        if excess(mid) < 0.0:
            inside = mid
        else:
            outside = mid
        if abs(outside - inside) < 1e-3:
            break
    offset = 0.5 * (inside + outside) - edge
    return [0.0, offset] if inward > 0 else [offset, 0.0]


def default_pois(lik: lk.Model, fixed_mh: float | None) -> tuple:
    return tuple(lik.extra_pois or ["mu"]) + (() if fixed_mh is not None else ("mH",))


def run_fit(lik: lk.Model, events: dict, pois=None, impacts: bool = True, significance: bool = True,
            fixed_mh: float | None = None) -> dict:
    pois = tuple(pois) if pois is not None else default_pois(lik, fixed_mh)
    fitter = Fitter(lik, events, fixed_mh)
    best = fitter.fit()
    bound = dict(fitter.bound)
    values = {n: float(best.values[n]) for n in fitter.names}
    errors = {n: float(best.errors[n]) for n in fitter.names}
    result = {"valid": bool(best.valid), "nll": float(best.fval), "values": values, "hesse": errors, "pois": {}}
    try:
        best.minos(*[poi for poi in pois if poi not in bound])
    except Exception as error:  # noqa: BLE001
        result["minos_error"] = str(error)
    # Statistical: every nuisance fixed at its best fit.
    stat = fitter.fit(start=values, fixed={n: values[n] for n in lik.nuisances})
    stat_bound = dict(fitter.bound)
    try:
        stat.minos(*[poi for poi in pois if poi not in stat_bound])
    except Exception as error:  # noqa: BLE001
        result["minos_stat_error"] = str(error)
    for poi in pois:
        me = best.merrors[poi] if poi in best.merrors else None
        ms = stat.merrors[poi] if poi in stat.merrors else None
        total = [float(me.lower), float(me.upper)] if me else [-errors[poi], errors[poi]]
        stat_interval = [float(ms.lower), float(ms.upper)] if ms else [-float(stat.errors[poi]), float(stat.errors[poi])]
        if poi in bound:
            # One-sided interval of a POI at its limit: the profile-likelihood edge where the NLL rises by 0.5.
            total = one_sided_interval(fitter, values, poi, bound[poi], best.fval, {})
            result.setdefault("at_limit", {})[poi] = bound[poi]
        if poi in stat_bound or poi in bound:
            stat_interval = one_sided_interval(fitter, values, poi, bound.get(poi, stat_bound.get(poi)), stat.fval,
                                               {n: values[n] for n in lik.nuisances})
        result["pois"][poi] = {"value": values[poi], "total": total, "stat": stat_interval}
    if impacts:
        result["impacts"] = {}
        for n in lik.nuisances:
            shifts = {}
            for sign in (1, -1):
                target = values[n] + sign * errors[n]
                m = fitter.fit(start=values, fixed={n: target})
                shifts[sign] = {poi: float(m.values[poi]) - values[poi] for poi in pois}
            result["impacts"][n] = {poi: [shifts[1][poi], shifts[-1][poi]] for poi in pois}
        for poi in pois:
            syst = math.sqrt(sum((0.5 * (abs(v[poi][0]) + abs(v[poi][1]))) ** 2 for v in result["impacts"].values()))
            result["pois"][poi]["syst"] = syst
    if significance and not lik.extra_pois:
        null = fitter.fit(start=values, fixed={"mu": 0.0, "mH": values["mH"]})
        q0 = max(2.0 * (null.fval - best.fval), 0.0)
        result["significance"] = {"q0": q0, "Z": math.sqrt(q0) if values["mu"] > 0 else 0.0}
    result["pulls"] = {n: values[n] for n in lik.nuisances}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="<select>/<model label> under production_v3/inference")
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--dimension", default="3D", choices=["1D", "2D", "2Dmass", "3D"])
    parser.add_argument("--no-refit", action="store_true")
    parser.add_argument("--set", default="b", choices=["a", "b"])
    parser.add_argument("--dataset", default="data", choices=["data", "asimov", "asimov_postfit"])
    parser.add_argument("--postfit-from", default=None, help="fit label whose best fit defines the post-fit Asimov dataset")
    parser.add_argument("--toys", type=int, default=0)
    parser.add_argument("--toy-mu", type=float, default=1.0)
    parser.add_argument("--toy-mh", type=float, default=125.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    parser.add_argument("--poi-scheme", default="inclusive")
    parser.add_argument("--fix-mh", type=float, default=None, help="fix m_H (GeV) instead of profiling it")
    parser.add_argument("--final-state", default=None, help="fit only the channels of this final state (4mu, 4e, 2e2mu)")
    parser.add_argument("--fiducial", type=Path, default=None,
                        help="fiducial_xsec.json of fiducial_xsec.py: the SM fiducial cross sections of the fid_int scheme")
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    out_dir = model_dir / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    model, events, signal_model, yr4 = load(model_dir, args.yr4)
    syst = json.loads(args.systematics.read_text(encoding="utf-8")) if args.systematics.exists() else {}
    fiducial_sm = None
    if args.fiducial:
        fiducial = json.loads((args.fiducial if args.fiducial.is_absolute() else REPO / args.fiducial).read_text(encoding="utf-8"))
        fiducial_sm = fiducial["total"]["sigma_fid_final_state"]
    lik = lk.Model(model, signal_model, yr4, args.dimension, not args.no_refit, args.set, syst, args.poi_scheme,
                   fiducial_sm=fiducial_sm)
    if args.final_state:
        # The likelihood of one final state: its channels only (the other final states' nuisances keep their constraints).
        lik.channels = [c for c in lik.channels if model["channels"][c]["final_state"] == args.final_state]
        if not lik.channels:
            raise SystemExit(f"no channel of final state {args.final_state}")
    variable, error = lik.variable, lik.error
    started = time.time()
    report = {"schema": "h4l_v3_fit/1", "model": str(model_dir), "dimension": args.dimension, "refit": not args.no_refit,
              "result_set": args.set, "dataset": args.dataset, "systematics": syst, "parameters": lik.parameters,
              "poi_scheme": args.poi_scheme, "fixed_mh": args.fix_mh, "final_state": args.final_state, "channels": lik.channels,
              "fiducial": str(args.fiducial) if args.fiducial else None, "fiducial_sm": fiducial_sm}
    if args.toys > 0:
        rng = np.random.default_rng(args.seed)
        p = np.zeros(len(lik.parameters))
        p[0], p[1] = args.toy_mu, args.toy_mh
        toys = []
        for t in range(args.toys):
            ev = toy_events(lik, p, rng)
            # Frequentist pseudo-experiment: the global observables drawn around the generated nuisances (0).
            lik.global_observables = rng.normal(size=len(lik.nuisances))
            r = run_fit(lik, ev, impacts=False, significance=False, fixed_mh=args.fix_mh)
            toys.append({"mu": r["pois"]["mu"], "mH": r["pois"]["mH"], "valid": r["valid"]})
            print(f"[toy {t}] mu {r['pois']['mu']['value']:.3f} mH {r['pois']['mH']['value']:.3f}", flush=True)
        report["toys"] = {"mu_true": args.toy_mu, "mh_true": args.toy_mh, "results": toys}
    else:
        if args.dataset == "asimov":
            p = lik.default_parameters(1.0, args.fix_mh if args.fix_mh is not None else 125.0)
            ev = asimov_events(lik, p)
        elif args.dataset == "asimov_postfit":
            # The post-fit expectation: the Asimov dataset at every parameter's best fit of the given fit.
            if not args.postfit_from:
                raise SystemExit("--dataset asimov_postfit needs --postfit-from")
            source = json.loads((model_dir / args.postfit_from / "fit.json").read_text(encoding="utf-8"))
            if source["parameters"] != lik.parameters:
                raise SystemExit(f"{args.postfit_from} has other parameters than this likelihood")
            p = np.array([source["fit"]["values"][n] for n in lik.parameters])
            report["postfit_from"] = str(model_dir / args.postfit_from / "fit.json")
            ev = asimov_events(lik, p)
            # The constraint centres at the fitted nuisance values, so that the post-fit Asimov returns the best fit.
            lik.global_observables = np.array([p[lik.index[n]] for n in lik.nuisances])
        else:
            ev = data_events(model, events, variable, error)
        report["fit"] = run_fit(lik, ev, fixed_mh=args.fix_mh)
        f = report["fit"]
        for poi, r in f["pois"].items():
            print(f"[fit] {poi} = {r['value']:.4f} total [{r['total'][0]:+.4f}, {r['total'][1]:+.4f}] stat [{r['stat'][0]:+.4f}, "
                  f"{r['stat'][1]:+.4f}] syst {r.get('syst', float('nan')):.4f}", flush=True)
        if "significance" in f:
            print(f"[fit] significance Z = {f['significance']['Z']:.2f}; valid {f['valid']}", flush=True)
    report["elapsed_s"] = time.time() - started
    out_dir.mkdir(parents=True)
    (out_dir / "fit.json").write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    print(f"[fit] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
