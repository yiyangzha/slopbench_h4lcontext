"""The binned H -> 4l likelihood: channels = final state x D_mass bin, m4l (Z1-refitted) bins of 0.5 GeV in 105-140 GeV.

Samples per channel: the signal (normfactor mu; its m4l template a double-sided Crystal Ball fitted to the 125 GeV
signal MC of the channel, peak position and width tracking m_H; yield sigma x BR(m_H) / sigma x BR(125) per mode from
YR4, the acceptance constant in m_H as the task prescribes, the window fraction of the line shape followed), qqZZ and
ggZZ (Bernstein polynomials fitted to their MC, normalized as 1000 L sigma_eff n_sel / n_presel with the lepton SFs)
and Z + X (the OS-method yield, Landau + exponential shape).  Nuisances (result set b, the uncertainties real for the
pseudo-data): lepton scale and resolution per flavour (signal shape, histosys), lepton efficiency per flavour (normsys on
signal, qqZZ, ggZZ), Z + X per final state (normsys).  The likelihood evaluates exactly the pyhf semantics
(histosys code 0, normsys code 1, unit Gaussian constraints) so that the exported MODEL.json at the fitted m_H is the
model that was fitted; m_H enters through the signal templates.
"""

from __future__ import annotations

import math

import numpy as np
from iminuit import Minuit
from scipy import special

from . import config as C

SQRT_HALF_PI = math.sqrt(0.5 * math.pi)
NUISANCES = ["scale_mu", "scale_e", "res_mu", "res_e", "eff_mu", "eff_e", "zx_4mu", "zx_4e", "zx_2e2mu"]
PARAMS = ["mu", "mH"] + NUISANCES
NORM_RANGE = (70.0, 200.0)


# ---------------------------------------------------------------------------------------------------------- shapes
def dcb_shape(t, al, nl, ar, nr):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.exp(-0.5 * t * t)
    left, right = t < -al, t > ar
    if np.any(left):
        out[left] = (nl / al) ** nl * math.exp(-0.5 * al * al) * np.power(nl / al - al - t[left], -nl)
    if np.any(right):
        out[right] = (nr / ar) ** nr * math.exp(-0.5 * ar * ar) * np.power(nr / ar - ar + t[right], -nr)
    return out


def dcb_primitive(t, al, nl, ar, nr):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    a_l = (nl / al) ** nl * math.exp(-0.5 * al * al)
    b_l = nl / al - al
    a_r = (nr / ar) ** nr * math.exp(-0.5 * ar * ar)
    b_r = nr / ar - ar
    f_left = a_l * (b_l + al) ** (1 - nl) / (nl - 1)
    erf_l = special.erf(-al / math.sqrt(2))
    f_right = f_left + SQRT_HALF_PI * (special.erf(ar / math.sqrt(2)) - erf_l)
    out = f_left + SQRT_HALF_PI * (special.erf(np.clip(t, -al, ar) / math.sqrt(2)) - erf_l)
    left, right = t < -al, t > ar
    if np.any(left):
        out[left] = a_l * np.power(b_l - t[left], 1 - nl) / (nl - 1)
    if np.any(right):
        out[right] = f_right + a_r / (nr - 1) * ((b_r + ar) ** (1 - nr) - np.power(b_r + t[right], 1 - nr))
    return out


def dcb_fractions(edges, mean, sigma, al, nl, ar, nr):
    """Probability of each bin under the DCB normalized over NORM_RANGE."""
    f = dcb_primitive((np.asarray(edges) - mean) / sigma, al, nl, ar, nr)
    lo, hi = dcb_primitive((np.array(NORM_RANGE) - mean) / sigma, al, nl, ar, nr)
    return np.diff(f) / (hi - lo)


def fit_dcb(m, w, lo, hi):
    """Weighted unbinned ML fit of a DCB normalized in [lo, hi]."""
    m, w = np.asarray(m, float), np.asarray(w, float)
    med = float(np.median(m))
    q16, q84 = np.percentile(m, [16, 84])

    def nll(mean, sigma, al, nl, ar, nr):
        norm = sigma * (dcb_primitive((hi - mean) / sigma, al, nl, ar, nr) - dcb_primitive((lo - mean) / sigma, al, nl, ar, nr))[0]
        if not norm > 0:
            return 1e30
        d = dcb_shape((m - mean) / sigma, al, nl, ar, nr) / norm
        if np.any(d <= 0):
            return 1e30
        return float(-np.sum(w * np.log(d)))

    mt = Minuit(nll, mean=med, sigma=max(0.5 * (q84 - q16) * 0.7, 0.5), al=1.0, nl=3.0, ar=1.5, nr=5.0)
    mt.errordef = Minuit.LIKELIHOOD
    mt.limits["mean"] = (med - 3, med + 3)
    mt.limits["sigma"] = (0.2, 6.0)
    mt.limits["al"] = (0.3, 4.0)
    mt.limits["ar"] = (0.3, 4.0)
    mt.limits["nl"] = (1.1, 50.0)
    mt.limits["nr"] = (1.1, 50.0)
    mt.migrad(ncall=10000)
    if not mt.valid:
        mt.migrad(ncall=10000)
    return {k: float(mt.values[k]) for k in ("mean", "sigma", "al", "nl", "ar", "nr")} | {"valid": bool(mt.valid)}


def bernstein_basis(u, degree):
    return np.stack([special.comb(degree, k) * u ** k * (1 - u) ** (degree - k) for k in range(degree + 1)], axis=-1)


def fit_bernstein(m, w, lo, hi, degree=3):
    u = (np.asarray(m, float) - lo) / (hi - lo)
    basis = bernstein_basis(u, degree)
    w = np.asarray(w, float)

    def nll(*c):
        c = np.abs(np.array(c))
        norm = c.sum() / (degree + 1)
        d = basis @ c / norm
        if np.any(d <= 0):
            return 1e30
        return float(-np.sum(w * np.log(d)))

    mt = Minuit(nll, *([1.0] * (degree + 1)), name=[f"c{k}" for k in range(degree + 1)])
    mt.errordef = Minuit.LIKELIHOOD
    mt.fixed["c0"] = True
    for k in range(1, degree + 1):
        mt.limits[f"c{k}"] = (0.0, 100.0)
    mt.migrad(ncall=10000)
    return [abs(float(mt.values[f"c{k}"])) for k in range(degree + 1)]


def bernstein_fractions(edges, coeffs, lo, hi, sub=20):
    """Bin probabilities of the Bernstein density normalized in [lo, hi] (numerical, fine sub-grid)."""
    edges = np.asarray(edges)
    fine = np.linspace(lo, hi, (len(edges) - 1) * sub + 1)
    centres = 0.5 * (fine[1:] + fine[:-1])
    d = bernstein_basis((centres - lo) / (hi - lo), len(coeffs) - 1) @ np.asarray(coeffs)
    per = d.reshape(len(edges) - 1, sub).sum(axis=1)
    return per / per.sum()


# ----------------------------------------------------------------------------------------------------- theory ratio
class XSBR:
    def __init__(self):
        y = C.CONSTANTS["yr4"]
        self.m = np.array(y["mass"])
        self.log = {mode: np.log(np.array(y["sigma_pb"][mode]) * np.array(y["br4l"])) for mode in C.SIGNAL_MODES}
        self.ref = {mode: float(np.interp(125.0, self.m, self.log[mode])) for mode in C.SIGNAL_MODES}

    def ratio(self, mode, mh):
        return math.exp(float(np.interp(mh, self.m, self.log[mode])) - self.ref[mode])


# ------------------------------------------------------------------------------------------------------ likelihood
class Likelihood:
    def __init__(self, channels: list, observed: list):
        self.channels = channels
        self.observed = [np.asarray(o, float) for o in observed]
        self.xsbr = XSBR()
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

    # signal templates of every channel at m_H: nominal and the +-1 sigma variations of the four shape nuisances
    def signal_templates(self, mh: float):
        key = float(mh)
        if key == self._cache_key:
            return self._cache
        out = []
        for ch in self.channels:
            s = ch["signal"]
            p = s["dcb"]
            scale = mh / 125.0
            mean0 = p["mean"] + (mh - 125.0)
            sigma0 = p["sigma"] * scale
            yield_total = sum(y * self.xsbr.ratio(mode, mh) for mode, y in s["yields_125"].items()) / s["window_fraction_125"]

            def tmpl(mean, sigma):
                return yield_total * dcb_fractions(ch["edges"], mean, sigma, p["al"], p["nl"], p["ar"], p["nr"])

            nom = tmpl(mean0, sigma0)
            var = {}
            for fl in ("mu", "e"):
                ks = s["kappa_scale"][fl]
                kr = s["kappa_res"][fl]
                var[f"scale_{fl}"] = (tmpl(mean0 * (1 + ks), sigma0), tmpl(mean0 * (1 - ks), sigma0)) if ks > 0 else None
                var[f"res_{fl}"] = (tmpl(mean0, sigma0 * (1 + kr)), tmpl(mean0, sigma0 * (1 - kr))) if kr > 0 else None
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
            k = z["kappa"]
            total = total + z["template"] * ((1 + k) ** t if t >= 0 else (1 - k) ** (-t))
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

    def saturated_q(self, p, observed=None, aux=None) -> float:
        observed = self.observed if observed is None else observed
        aux = self.aux if aux is None else aux
        q = 0.0
        for nu, n in zip(self.expected(np.asarray(p, float)), observed):
            q += float(np.sum(nu - n + np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / nu), 0.0)))
        return 2.0 * q + float(np.sum((np.asarray(p[2:], float) - aux) ** 2))


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
