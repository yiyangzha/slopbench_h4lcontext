"""Reducible background (Z + X) by the OS method of AN-16-442 7.2.1: fake rates from Z + 1 loose lepton in data (the
ZZ MC prompt contribution subtracted), applied to the 2P2F and 3P1F control regions:
    N = (1 - N_ZZ(3P1F) / N(3P1F)) sum_3P1F f / (1 - f) - sum_2P2F f3 f4 / ((1 - f3)(1 - f4)),
per final state inside the fit window.  The m4l shape per final state is the main analysis's Landau + exponential; the
D_mass fractions come from the control-region events.  Uncertainty: the control-region statistics in quadrature with
the relative MC-closure systematic of the method (main analysis)."""

from __future__ import annotations

import math

import numpy as np
from scipy import stats

from . import config as C
from .selection import fake_rate_bin

FS_NAMES = ("4mu", "4e", "2e2mu")


def fake_rates(data_zl: dict, mc_zl: list) -> dict:
    """Per flavour code: the fake rate per bin (pT bin x 2 + eta bin), numerator and denominator after the prompt
    subtraction.  mc_zl: [(rows, weight per row)]."""
    out = {}
    for code in (13, 11):
        edges = [e for e in C.CONSTANTS["zx"]["fake_rate_pt_edges"][str(code)] if e is not None]
        nb = 2 * len(edges)
        num, den = np.zeros(nb), np.zeros(nb)
        sel = (data_zl["flavour"] == code) & data_zl["sip_ok"]
        k = fake_rate_bin(np.full(sel.sum(), code), data_zl["pt"][sel].astype(float), data_zl["abs_eta"][sel].astype(float))
        np.add.at(den, k, 1.0)
        np.add.at(num, k[data_zl["pass"][sel]], 1.0)
        num_mc, den_mc = np.zeros(nb), np.zeros(nb)
        for rows, w in mc_zl:
            s = (rows["flavour"] == code) & rows["sip_ok"]
            if not np.any(s):
                continue
            kk = fake_rate_bin(np.full(s.sum(), code), rows["pt"][s].astype(float), rows["abs_eta"][s].astype(float))
            ws = np.asarray(w)[s]
            np.add.at(den_mc, kk, ws)
            passing = rows["pass"][s]
            np.add.at(num_mc, kk[passing], ws[passing])
        n_sub, d_sub = np.clip(num - num_mc, 0.0, None), np.clip(den - den_mc, 1e-9, None)
        rate = np.clip(np.where(den > 0, n_sub / d_sub, 0.0), 1e-4, 0.9)
        out[code] = {"rate": rate, "num": num, "den": den, "num_prompt_mc": num_mc, "den_prompt_mc": den_mc}
    return out


def weights(rates: dict, flavour, pt, abs_eta) -> np.ndarray:
    f = np.zeros(len(pt))
    for code, r in rates.items():
        sel = flavour == code
        f[sel] = r["rate"][fake_rate_bin(flavour[sel], pt[sel], abs_eta[sel])]
    return f


def landau_exp_fraction(shape: dict, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Probability of [a, b] under the Landau + exponential shape (components normalized in its fit range)."""
    r_lo, r_hi = shape["range"]
    mpv, width, f, slope = shape["mpv"], shape["width"], shape["exp_fraction"], shape["exp_slope"]
    n_l = stats.landau.cdf(r_hi, loc=mpv, scale=width) - stats.landau.cdf(r_lo, loc=mpv, scale=width)
    n_e = (math.exp(-slope * (r_lo - 70.0)) - math.exp(-slope * (r_hi - 70.0))) / slope
    land = (stats.landau.cdf(b, loc=mpv, scale=width) - stats.landau.cdf(a, loc=mpv, scale=width)) / n_l
    expo = (np.exp(-slope * (a - 70.0)) - np.exp(-slope * (b - 70.0))) / slope / n_e
    return (1 - f) * land + f * expo


def os_estimate(cr3: dict, cr2: dict, zz3: list, rates: dict) -> dict:
    """cr3 / cr2: data 3P1F / 2P2F candidates in the window (final_state, fail_flavour, fail_pt, fail_eta for the failing
    legs); zz3: [(3P1F candidates of a ZZ MC sample, per-event weight)].  Returns the yield and the relative statistical
    uncertainty per final state."""
    out = {}
    for code, name in enumerate(FS_NAMES):
        s3 = cr3["final_state"] == code
        f3 = weights(rates, cr3["fail_flavour"][s3], cr3["fail_pt"][s3], cr3["fail_eta"][s3])
        w3 = f3 / (1 - f3)
        s2 = cr2["final_state"] == code
        fa = weights(rates, cr2["fail_flavour"][s2][:, 0], cr2["fail_pt"][s2][:, 0], cr2["fail_eta"][s2][:, 0])
        fb = weights(rates, cr2["fail_flavour"][s2][:, 1], cr2["fail_pt"][s2][:, 1], cr2["fail_eta"][s2][:, 1])
        w2 = fa * fb / ((1 - fa) * (1 - fb))
        n3 = float(s3.sum())
        n_zz = float(sum(np.sum(w[c["final_state"] == code]) for c, w in zz3))
        zz_factor = max(1.0 - n_zz / n3, 0.0) if n3 > 0 else 0.0
        term3 = zz_factor * float(np.sum(w3))
        term2 = -float(np.sum(w2))
        value = term3 + term2
        var = zz_factor ** 2 * float(np.sum(w3 ** 2)) + float(np.sum(w2 ** 2))
        out[name] = {"yield_window": value, "stat": math.sqrt(var), "term_3p1f": term3, "term_2p2f": term2, "n_3p1f": n3,
                     "n_2p2f": float(s2.sum()), "n_zz_3p1f": n_zz, "zz_factor": zz_factor}
    return out
