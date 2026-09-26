"""Pass 2 of the evaluation run: calibrated event selection, templates, fit, validation and outputs."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from . import config as C
from . import model as M
from . import selection as S
from . import zx as Z
from .tnp import ScaleFactors

FS_NAMES = ("4mu", "4e", "2e2mu")


def lepton_table(rec: dict):
    lep = {k[2:]: v for k, v in rec.items() if k.startswith("l_")}
    for k in ("pt_raw", "eta", "eta_sc", "phi", "rel_err", "iso_fsr", "sip", "fsr_pt", "fsr_eta", "fsr_phi", "g"):
        lep[k] = lep[k].astype(np.float64)
    counts = rec["n_lep"].astype(np.int64)
    return lep, counts


def abs_eta_of(lep):
    return np.where(lep["flavour"] == 13, np.abs(lep["eta"]), np.abs(lep["eta_sc"]))


def calibrated_pt(lep, models, is_mc, smear_shift=None):
    """Data: pT / (1 + s); MC: pT (1 + r N) (r optionally shifted per flavour: {code: delta_r})."""
    pt = lep["pt_raw"].copy()
    ae = abs_eta_of(lep)
    for code, model in models.items():
        sel = lep["flavour"] == code
        if is_mc:
            r = model.smear(pt[sel], ae[sel])
            if smear_shift and code in smear_shift:
                r = np.clip(r + smear_shift[code], 0.0, None)
            pt[sel] = pt[sel] * (1.0 + r * lep["g"][sel])
        else:
            pt[sel] = model.corrected_data_pt(pt[sel], ae[sel])
    return np.where(pt > 0, pt, 1e-6)


def select_events(rec: dict, models: dict, is_mc: bool, control: bool, zl: bool, smear_shift=None) -> dict:
    lep, counts = lepton_table(rec)
    pt = calibrated_pt(lep, models, is_mc, smear_shift)
    selected = S.selected_mask(lep, pt)
    keep = S.cross_clean(lep, counts, selected) & S.loose_threshold_mask(lep, pt)
    selected &= keep
    out = {"lep": lep, "counts": counts, "pt": pt}
    sr = S.quad_candidates(lep, counts, pt, selected, selected, selected, C.SELECT["z2_low"])
    out["sr"] = sr
    if control:
        pool = keep & (lep["sip"] < C.SELECT["max_sip"])
        in_sr = np.zeros(len(counts), bool)
        in_sr[sr["event"]] = True
        for name, nfail in (("cr3", 1), ("cr2", 2)):
            c = S.quad_candidates(lep, counts, pt, pool, selected, pool, C.SELECT["z2_low"], n_fail_z2=nfail, fail=~selected)
            keep_ev = ~in_sr[c["event"]]
            out[name] = {k: v[keep_ev] for k, v in c.items()}
    if zl:
        out["zl"] = S.z_plus_lepton(lep, counts, pt, keep, selected, rec["e_met"].astype(float))
    return out


def with_kinematics(sel: dict, key: str) -> dict:
    c = sel[key]
    kin = S.candidate_kinematics(sel["lep"], c["legs"], sel["pt"])
    c = dict(c)
    c.update(kin)
    c["final_state"] = S.final_state(sel["lep"], c["legs"])
    return c


def failing_legs(sel: dict, c: dict, n: int):
    """Flavour, pT, |eta| of the failing Z2 legs of CR candidates (first n failing legs among legs 2, 3)."""
    lep, pt = sel["lep"], sel["pt"]
    selected = S.selected_mask(lep, pt)
    legs = c["legs"][:, 2:]
    fails = ~selected[legs]
    fl, fp, fe = [], [], []
    ae = abs_eta_of(lep)
    order = np.argsort(~fails, axis=1, kind="stable")
    for j in range(n):
        idx = legs[np.arange(len(legs)), order[:, j]]
        fl.append(lep["flavour"][idx])
        fp.append(pt[idx])
        fe.append(ae[idx])
    if n == 1:
        return fl[0], fp[0], fe[0]
    return np.stack(fl, 1), np.stack(fp, 1), np.stack(fe, 1)


def event_sf(sel: dict, c: dict, sfs: ScaleFactors):
    lep, pt = sel["lep"], sel["pt"]
    legs = c["legs"]
    flat = legs.ravel()
    sf, rs, rm, kk = sfs.lookup(lep["flavour"][flat], pt[flat], abs_eta_of(lep)[flat])
    n = len(legs)
    return {"sf": sf.reshape(n, 4).prod(axis=1), "relstat": rs.reshape(n, 4), "relmodel": rm.reshape(n, 4),
            "bin": kk.reshape(n, 4), "flavour": lep["flavour"][flat].reshape(n, 4)}


def efficiency_kappa(w, sfinfo, sel_mask):
    """Relative yield uncertainty per flavour from the SF uncertainties: statistical parts independent between bins,
    fit-model parts coherent."""
    out = {}
    for fl, code in (("mu", 13), ("e", 11)):
        y = float(np.sum(w[sel_mask]))
        if y <= 0:
            out[fl] = 0.0
            continue
        ww = w[sel_mask][:, None]
        is_f = sfinfo["flavour"][sel_mask] == code
        model = float(np.sum(ww * np.where(is_f, sfinfo["relmodel"][sel_mask], 0.0)))
        bins = sfinfo["bin"][sel_mask]
        contrib = (ww * np.where(is_f, sfinfo["relstat"][sel_mask], 0.0)).ravel()
        keyb = np.where(is_f, bins, -1).ravel()
        ok = keyb >= 0
        per_bin = np.bincount(keyb[ok], weights=contrib[ok]) if np.any(ok) else np.zeros(1)
        stat = float(np.sqrt(np.sum(per_bin ** 2)))
        out[fl] = math.hypot(stat, model) / y
    return out


def robust_width(m, w):
    order = np.argsort(m)
    cw = np.cumsum(w[order])
    cw /= cw[-1]
    q16 = m[order][np.searchsorted(cw, 0.16)]
    q84 = m[order][np.searchsorted(cw, 0.84)]
    return 0.5 * (q84 - q16)
