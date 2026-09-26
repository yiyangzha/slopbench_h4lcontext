"""Channels of the binned likelihood: D_mass bins, signal line shapes and yields, backgrounds, nuisance magnitudes."""

from __future__ import annotations

import math

import numpy as np

from . import analysis as A
from . import config as C
from . import model as M
from . import selection as S
from . import zx as Z

FLAVOUR_CODES = {"mu": 13, "e": 11}


def weighted_quantiles(x, w, qs):
    order = np.argsort(x)
    cw = np.cumsum(w[order])
    cw /= cw[-1]
    return [float(x[order][min(np.searchsorted(cw, q), len(x) - 1)]) for q in qs]


def scale_fractions(c: dict) -> dict:
    """Per final state and flavour: d ln m4l / d ln pT of that flavour's leptons (weighted signal average)."""
    sel = c["sel"]
    lep, pt, legs = sel["lep"], sel["pt"], c["legs"]
    flat = legs.ravel()
    base = S._mass(S.dressed_vectors(lep, flat, pt[flat]).reshape(len(legs), 4, 4).sum(axis=1))
    out = {}
    eps = 1e-3
    for fl, code in FLAVOUR_CODES.items():
        moved = pt[flat] * np.where(lep["flavour"][flat] == code, 1 + eps, 1.0)
        m = S._mass(S.dressed_vectors(lep, flat, moved).reshape(len(legs), 4, 4).sum(axis=1))
        out[fl] = (m / base - 1.0) / eps
    return out


def build_channels(ds, samples, mc_sr, data_sr, zx_est, cr3, cr2, models, cal_summary, sfs, norm, log):
    lo, hi = C.WINDOW
    edges = np.round(np.arange(lo, hi + 1e-9, C.BIN_WIDTH), 6)
    roles = {n: p["role"] for n, p in ds["processes"].items()}
    modes = {n: p["mode"] for n, p in ds["processes"].items()}
    # Window selection and D_mass of every candidate.
    for c in list(mc_sr.values()) + [data_sr, cr3, cr2]:
        c["in_window"] = (c["m4l_refit"] > lo) & (c["m4l_refit"] < hi)
        c["dmass"] = c["m4l_refit_err"] / c["m4l_refit"]
    sig = [n for n in mc_sr if roles[n] == "signal"]
    for n in sig:
        if modes[n] is None:
            raise SystemExit(f"signal sample {n}: unknown production mode")
    # Per final state D_mass edges: weighted terciles of the signal MC.
    dedges = {}
    for k, fs in enumerate(A.FS_NAMES):
        d = np.concatenate([mc_sr[n]["dmass"][mc_sr[n]["in_window"] & (mc_sr[n]["final_state"] == k)] for n in sig])
        w = np.concatenate([mc_sr[n]["w"][mc_sr[n]["in_window"] & (mc_sr[n]["final_state"] == k)] for n in sig])
        q = weighted_quantiles(d, w, [(i + 1) / C.N_DMASS_BINS for i in range(C.N_DMASS_BINS - 1)])
        dedges[fs] = [0.0] + q + [math.inf]

    def dbin(c, k, fs):
        return np.clip(np.searchsorted(np.array(dedges[fs]), c["dmass"], side="right") - 1, 0, C.N_DMASS_BINS - 1)

    # Scale fractions and resolution magnitudes per final state (signal MC).
    frac = {fs: {"mu": [], "e": [], "w": []} for fs in A.FS_NAMES}
    width_nom = {fs: {} for fs in A.FS_NAMES}
    kappa_res = {fs: {"mu": 0.0, "e": 0.0} for fs in A.FS_NAMES}
    for n in sig:
        c = mc_sr[n]
        f = scale_fractions(c)
        for k, fs in enumerate(A.FS_NAMES):
            m = c["in_window"] & (c["final_state"] == k)
            frac[fs]["mu"].append(f["mu"][m])
            frac[fs]["e"].append(f["e"][m])
            frac[fs]["w"].append(c["w"][m])
    scale_frac = {}
    for fs in A.FS_NAMES:
        w = np.concatenate(frac[fs]["w"])
        scale_frac[fs] = {fl: float(np.sum(np.concatenate(frac[fs][fl]) * w) / np.sum(w)) for fl in ("mu", "e")}
    for fl, code in FLAVOUR_CODES.items():
        name = "muon" if fl == "mu" else "electron"
        dr = cal_summary[name]["smear"]["unc"]
        delta = max(3.0 * dr, 0.004)
        m_nom, m_var, w_all, fs_all = [], [], [], []
        for n in sig:
            c = mc_sr[n]
            sel = c["sel"]
            pt_var = A.calibrated_pt(sel["lep"], models, True, smear_shift={code: delta})
            kin = S.candidate_kinematics(sel["lep"], c["legs"], pt_var)
            keep = c["in_window"]
            m_nom.append(c["m4l_refit"][keep])
            m_var.append(kin["m4l_refit"][keep])
            w_all.append(c["w"][keep])
            fs_all.append(c["final_state"][keep])
        m_nom, m_var, w_all, fs_all = map(np.concatenate, (m_nom, m_var, w_all, fs_all))
        for k, fs in enumerate(A.FS_NAMES):
            s = fs_all == k
            if s.sum() < 50:
                continue
            w0 = A.robust_width(m_nom[s], w_all[s])
            w1 = A.robust_width(m_var[s], w_all[s])
            kappa_res[fs][fl] = abs(w1 / w0 - 1.0) * dr / delta
            width_nom[fs][fl] = (w0, w1)
    scale_unc = {"mu": cal_summary["muon"]["scale_shift"]["unc"], "e": cal_summary["electron"]["scale_shift"]["unc"]}

    channels, observed, meta = [], [], {"dmass_edges": dedges, "scale_fractions": scale_frac, "kappa_res": kappa_res,
                                        "bin_edges": edges.tolist(), "channels": []}
    zx_cfg = C.CONSTANTS["zx"]
    # Z+X D_mass fractions per final state: the OS control-region events (2P2F + 3P1F) with m4l_refit in 70-400 GeV.
    zx_dfrac = {}
    for k, fs in enumerate(A.FS_NAMES):
        counts = np.zeros(C.N_DMASS_BINS)
        for c in (cr3, cr2):
            m = (c["final_state"] == k) & (c["m4l_refit"] > 70.0) & (c["m4l_refit"] < 400.0)
            if np.any(m):
                counts += np.bincount(dbin({"dmass": c["dmass"][m]}, k, fs), minlength=C.N_DMASS_BINS)
        if counts.sum() < 30:
            # Too few control events: the qqZZ D_mass fractions.
            qq = [n for n in mc_sr if roles[n] == "qqZZ"]
            counts = np.zeros(C.N_DMASS_BINS)
            for n in qq:
                c = mc_sr[n]
                m = c["in_window"] & (c["final_state"] == k)
                counts += np.bincount(dbin({"dmass": c["dmass"][m]}, k, fs), weights=c["w"][m], minlength=C.N_DMASS_BINS)
        zx_dfrac[fs] = counts / counts.sum()
    meta["zx_dmass_fractions"] = {fs: v.tolist() for fs, v in zx_dfrac.items()}
    for k, fs in enumerate(A.FS_NAMES):
        shape_frac = Z.landau_exp_fraction(zx_cfg["shapes"][fs], edges[:-1], edges[1:])
        shape_frac = shape_frac / shape_frac.sum()
        zy = zx_est[fs]
        z_yield = zy["yield_window"]
        z_kappa = math.hypot(zy["stat"] / z_yield, zx_cfg["relative_systematic"]) if z_yield > 0 else 1.0
        if z_yield <= 0:
            z_yield = 1e-3
        z_kappa = min(z_kappa, 0.99)
        for b in range(C.N_DMASS_BINS):
            ch = {"name": f"ch_{fs}_dmass{b}", "fs": fs, "dbin": b, "edges": edges}
            # signal
            ms, ws, md = [], [], []
            yields = {mode: 0.0 for mode in C.SIGNAL_MODES}
            eff_rows = []
            for n in sig:
                c = mc_sr[n]
                m = c["in_window"] & (c["final_state"] == k) & (dbin(c, k, fs) == b)
                ms.append(c["m4l_refit"][m])
                ws.append(c["w"][m])
                yields[modes[n]] += float(np.sum(c["w"][m]))
                eff_rows.append((c["w"], c["sfinfo"], m))
            ms, ws = np.concatenate(ms), np.concatenate(ws)
            dcb = M.fit_dcb(ms, ws, lo, hi)
            wf = float(np.sum(M.dcb_fractions(edges, dcb["mean"], dcb["sigma"], dcb["al"], dcb["nl"], dcb["ar"], dcb["nr"])))
            ke = _merge_eff(eff_rows)
            ch["signal"] = {"dcb": dcb, "yields_125": yields, "window_fraction_125": wf, "kappa_eff": ke,
                            "kappa_scale": {fl: abs(scale_frac[fs][fl]) * scale_unc[fl] for fl in ("mu", "e")},
                            "kappa_res": dict(kappa_res[fs])}
            for role in ("qqZZ", "ggZZ"):
                names = [n for n in mc_sr if roles[n] == role]
                mb, wb, rows = [], [], []
                for n in names:
                    c = mc_sr[n]
                    m = c["in_window"] & (c["final_state"] == k) & (dbin(c, k, fs) == b)
                    mb.append(c["m4l_refit"][m])
                    wb.append(c["w"][m])
                    rows.append((c["w"], c["sfinfo"], m))
                mb = np.concatenate(mb) if mb else np.zeros(0)
                wb = np.concatenate(wb) if wb else np.zeros(0)
                total = float(np.sum(wb))
                if len(mb) >= 30:
                    coeff = M.fit_bernstein(mb, wb, lo, hi)
                    tmpl = total * M.bernstein_fractions(edges, coeff, lo, hi)
                else:
                    coeff = None
                    tmpl = np.full(len(edges) - 1, max(total, 1e-6) / (len(edges) - 1))
                ch[role] = {"template": np.clip(tmpl, 1e-9, None), "kappa_eff": _merge_eff(rows) if rows else {"mu": 0.0, "e": 0.0},
                            "bernstein": coeff, "mc_events": int(len(mb))}
            ch["zx"] = {"template": np.clip(z_yield * zx_dfrac[fs][b] * shape_frac, 1e-9, None), "kappa": z_kappa,
                        "nuisance": f"zx_{fs}"}
            dm = (data_sr["final_state"] == k) & data_sr["in_window"] & (dbin(data_sr, k, fs) == b)
            obs, _ = np.histogram(data_sr["m4l_refit"][dm], bins=edges)
            channels.append(ch)
            observed.append(obs.astype(float))
            meta["channels"].append({"name": ch["name"], "observed": int(obs.sum()), "signal_125": sum(yields.values()),
                                     "qqZZ": float(ch["qqZZ"]["template"].sum()), "ggZZ": float(ch["ggZZ"]["template"].sum()),
                                     "zx": float(ch["zx"]["template"].sum()), "dcb": dcb, "kappa_eff": ke})
            log(f"[model] {ch['name']}: observed {int(obs.sum())}, signal {sum(yields.values()):.2f}, qqZZ "
                f"{ch['qqZZ']['template'].sum():.2f}, ggZZ {ch['ggZZ']['template'].sum():.2f}, Z+X {ch['zx']['template'].sum():.2f}; "
                f"DCB mean {dcb['mean']:.2f} sigma {dcb['sigma']:.2f}")
    return channels, observed, meta


def _merge_eff(rows):
    """kappa_eff per flavour of the union of several samples' selected events."""
    w = np.concatenate([r[0][r[2]] for r in rows])
    info = {k: np.concatenate([r[1][k][r[2]] for r in rows]) for k in ("relstat", "relmodel", "bin", "flavour")}
    if len(w) == 0:
        return {"mu": 0.0, "e": 0.0}
    return A.efficiency_kappa(w, info, np.ones(len(w), bool))
