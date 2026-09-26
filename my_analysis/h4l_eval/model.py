"""The binned H -> 4l likelihood: channels = final state x D_mass bin, m4l (Z1-refitted) bins of 0.5 GeV in 105-140 GeV.

Samples per channel (templates.py): the signal (normfactor mu; per mode the yield at 125 GeV times sigma_eff(m_H) /
sigma_eff(125) (YR4, the main analysis's table and extrapolation) times the acceptance ratio A(m_H) of the scale morphing,
moved towards the shift morphing by the morphing nuisance; the resonant DCB of the channel with mean and width scaled by
m_H / 125, VH with its non-resonant Landau), qqZZ, ggZZ and Z + X.  Nuisances (result set b): lepton scale and resolution
per flavour (signal shape, histosys), the A x eff(m_H) morphing (histosys), lepton efficiency per flavour (normsys on
signal, qqZZ, ggZZ), Z + X per final state (asymmetric normsys).  The likelihood evaluates exactly the pyhf semantics
(histosys code 0, normsys code 1, unit Gaussian constraints) so that the exported MODEL.json at the fitted m_H is the
model that was fitted; m_H enters through the signal templates.
"""

from __future__ import annotations

import math

import numpy as np
from iminuit import Minuit

from . import config as C
from . import shapes as SH

NUISANCES = ["scale_mu", "scale_e", "res_mu", "res_e", "eff_mu", "eff_e", "morph", "zx_4mu", "zx_4e", "zx_2e2mu"]
PARAMS = ["mu", "mH"] + NUISANCES
SHAPE_NUISANCES = ("scale_mu", "scale_e", "res_mu", "res_e", "morph")


class SigmaEff:
    """sigma_eff(mode, m_H) / sigma_eff(mode, 125) from the YR4 table (log-quadratic outside it; the main analysis's
    h4l_likelihood.SigmaEff)."""

    def __init__(self):
        self.funcs = {}
        for mode, t in C.CONSTANTS["sigma_eff"].items():
            self.funcs[mode] = (np.array(t["mass"]), np.array(t["ratio"]), np.array(t["log_quadratic"]))

    def ratio(self, mode: str, mh: float) -> float:
        masses, ratios, coeff = self.funcs[mode]
        if masses[0] <= mh <= masses[-1]:
            return float(np.interp(mh, masses, ratios))
        return float(np.exp(np.polyval(coeff, mh - 125.0)))


def signal_template(ch: dict, mh: float, xs: SigmaEff, scale: float = 1.0, res: float = 1.0, morph: float = 0.0) -> np.ndarray:
    lo, hi = C.WINDOW
    s = ch["signal"]
    v = s["dcb"]
    k = mh / 125.0
    resonant = SH.dcb_fractions(ch["edges"], v, v["mean"] * k * scale, v["width"] * k * max(res, 0.2), lo, hi)
    total = np.zeros(len(ch["edges"]) - 1)
    for mode, y in s["yields_125"].items():
        a_sc = float(np.interp(mh, C.MH_GRID, s["a_scale"][mode]))
        a_sh = float(np.interp(mh, C.MH_GRID, s["a_shift"][mode]))
        yy = y * xs.ratio(mode, mh) * (a_sc + morph * (a_sh - a_sc))
        if mode == "VH":
            total += yy * (s["f_res"] * resonant + (1.0 - s["f_res"]) * s["landau_fractions"])
        else:
            total += yy * resonant
    return total


class Likelihood:
    def __init__(self, channels: list, observed: list):
        self.channels = channels
        self.observed = [np.asarray(o, float) for o in observed]
        self.xs = SigmaEff()
        self.aux = np.zeros(len(NUISANCES))
        self.index = {n: i for i, n in enumerate(PARAMS)}
        self._cache_key = None
        self._cache = None
        # Lower limit of mu: the expected count of every bin stays positive for every m_H of the range (the Poisson model is
        # defined for positive expectations only); 90 % of the most restrictive -B/S.
        worst = -20.0
        for mh in np.arange(C.MH_LIMITS[0], C.MH_LIMITS[1] + 1e-9, 0.5):
            for ch, (nom, _) in zip(self.channels, self.signal_templates(float(mh))):
                bkg = ch["qqZZ"]["template"] + ch["ggZZ"]["template"] + ch["zx"]["template"]
                ok = nom > 1e-12
                if np.any(ok):
                    worst = max(worst, float(np.max(-bkg[ok] / nom[ok])))
        self.mu_low = max(-20.0, 0.9 * worst)

    def signal_templates(self, mh: float):
        """Per channel at m_H: the nominal signal template and the +-1 sigma templates of the shape nuisances."""
        key = float(mh)
        if key == self._cache_key:
            return self._cache
        out = []
        for ch in self.channels:
            s = ch["signal"]
            nom = signal_template(ch, mh, self.xs)
            var = {}
            for fl in ("mu", "e"):
                ks, kr = s["scale_kappa"][fl], s["res_kappa"][fl]
                var[f"scale_{fl}"] = (signal_template(ch, mh, self.xs, scale=1 + ks), signal_template(ch, mh, self.xs, scale=1 - ks)) if ks > 0 else None
                var[f"res_{fl}"] = (signal_template(ch, mh, self.xs, res=1 + kr), signal_template(ch, mh, self.xs, res=1 - kr)) if kr > 0 else None
            up, down = signal_template(ch, mh, self.xs, morph=1.0), signal_template(ch, mh, self.xs, morph=-1.0)
            var["morph"] = (up, down) if np.any(np.abs(up - nom) > 1e-12 * np.maximum(nom, 1e-30)) else None
            out.append((nom, var))
        self._cache_key, self._cache = key, out
        return out

    def expected(self, p: np.ndarray) -> list:
        mu, mh = p[0], p[1]
        th = dict(zip(NUISANCES, p[2:]))
        sig = self.signal_templates(mh)
        out = []
        for ch, (nom, var) in zip(self.channels, sig):
            s = nom.copy()
            for name, v in var.items():
                if v is None:
                    continue
                t = th[name]
                up, down = v
                s += t * (up - nom) if t >= 0 else t * (nom - down)
            f_sig = 1.0
            for fl in ("mu", "e"):
                k = ch["signal"]["kappa_eff"][fl]
                t = th[f"eff_{fl}"]
                f_sig *= (1 + k) ** t if t >= 0 else (1 - k) ** (-t)
            total = mu * s * f_sig
            for name in ("qqZZ", "ggZZ"):
                b = ch[name]
                f = 1.0
                for fl in ("mu", "e"):
                    k = b["kappa_eff"][fl]
                    t = th[f"eff_{fl}"]
                    f *= (1 + k) ** t if t >= 0 else (1 - k) ** (-t)
                total = total + b["template"] * f
            z = ch["zx"]
            t = th[z["nuisance"]]
            k_lo, k_hi = z["kappa"]
            # The main analysis's asymmetric log-normal: +1 sigma at central (1 + kappa_high), -1 sigma at central (1 - kappa_low).
            total = total + z["template"] * ((1 + k_hi) ** t if t >= 0 else (1 - k_lo) ** (-t))
            out.append(total)
        return out

    def nll(self, p, observed=None, aux=None) -> float:
        observed = self.observed if observed is None else observed
        aux = self.aux if aux is None else aux
        total = 0.0
        for nu, n in zip(self.expected(np.asarray(p, float)), observed):
            if np.any(nu <= 0) or not np.all(np.isfinite(nu)):
                return 1e30
            total += float(np.sum(nu - n * np.log(nu)))
        th = np.asarray(p[2:], float)
        return total + 0.5 * float(np.sum((th - aux) ** 2))


def rebin_gof(values: list, step_bins: int) -> np.ndarray:
    """Counts per channel summed into groups of step_bins consecutive m4l bins (the goodness-of-fit binning)."""
    return np.concatenate([np.add.reduceat(np.asarray(v, float), np.arange(0, len(v), step_bins)) for v in values])


def saturated_q(nu: np.ndarray, n: np.ndarray) -> float:
    return float(2.0 * np.sum(nu - n + np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / nu), 0.0)))


# --------------------------------------------------------------------------------------------------------- fitting
class Fitter:
    def __init__(self, lik: Likelihood):
        self.lik = lik

    def minuit(self, start, fixed=None, observed=None, aux=None):
        m = Minuit(lambda x: self.lik.nll(x, observed, aux), np.asarray(start, float), name=PARAMS)
        m.errordef = Minuit.LIKELIHOOD
        m.limits["mu"] = (self.lik.mu_low, 50.0)
        m.limits["mH"] = C.MH_LIMITS
        for n in NUISANCES:
            m.limits[n] = (-5.0, 5.0)
        for n, v in (fixed or {}).items():
            m.values[n] = v
            m.fixed[n] = True
        m.strategy = 1
        m.tol = C.PRECISION["minuit_tolerance"]
        return m

    def fit(self, start, fixed=None, observed=None, aux=None):
        m = self.minuit(start, fixed, observed, aux)
        m.migrad(ncall=20000)
        if not m.valid:
            m.migrad(ncall=20000)
        if not m.valid:
            m.simplex(ncall=20000)
            m.migrad(ncall=20000)
        m.hesse()
        return m

    def global_fit(self, observed=None, aux=None, step=0.5):
        """Profile scan over the m_H grid (mu and the nuisances profiled), then a free fit from the best grid point."""
        start = np.zeros(len(PARAMS))
        start[0], start[1] = 1.0, 125.0
        best, scan = None, []
        values = start.copy()
        for mh in np.arange(C.MH_LIMITS[0], C.MH_LIMITS[1] + 1e-9, step):
            values[1] = mh
            m = self.fit(values, fixed={"mH": float(mh)}, observed=observed, aux=aux)
            scan.append((float(mh), float(m.fval), float(m.values["mu"])))
            if best is None or m.fval < best.fval:
                best = m
            values = np.array(m.values)
        free = self.fit(np.array(best.values), observed=observed, aux=aux)
        return free, scan
