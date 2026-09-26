"""Stage 4b: optimization of the lepton selection on MC only (never on data).

    pixi run py -- analysis_v3/reconstruction/scripts/optimize_selection.py --select v3 --label scan_i1

Objective (AGENTS.md, PLAN.md 4b): the joint expected precision of mu and m_H
at the luminosity of the data sample, mu = 1, m_H = 125 GeV,
    J = (sigma_mu / sigma_mu,ref)^2 + (sigma_mH / sigma_mH,ref)^2,
the reference being the AN starting point.  The precisions come from the
Fisher information of a binned Poisson model per final state (4mu, 4e,
2e2mu) in m4l (the Z1-refitted mass, 105-140 GeV, 0.25 GeV bins) x D_bkg^kin
(five bins, the quintiles of the reference signal) x the relative per-event
mass error (three bins, the tertiles of the reference signal per final state):
    nu = mu s(m4l - (m_H - 125)) + b_ZZ + b_ZX,
    I = sum_bins (d nu / d theta)(d nu / d theta)^T / nu,  theta = (mu, m_H),
    sigma = sqrt(diag(I^-1))  (each parameter profiled over the other).
s is the ggH + VBF + VH MC at 125 GeV smoothed with a Gaussian kernel of
0.3 GeV (its m_H derivative is minus its m4l derivative); b_ZZ the qqZZ and
ggZZ MC (kernel 2 GeV); b_ZX the fake-rate method applied to the DY and
TTbar MC: fake rates from the Z + 1 loose lepton rows (|m_Z1 - m_Z| < 7 GeV,
MET < 25 GeV, denominator loose + SIP, the user decision), in bins of probe
pT x barrel/endcap per flavour, applied to the 3P1F and 2P2F rows of the
same samples, N = sum_3P1F f/(1-f) - sum_2P2F f3 f4/((1-f3)(1-f4)); its m4l
shape is a 5 GeV histogram of the weighted rows interpolated linearly, its
D_kin and mass-error composition that of the rows in the window.
Options act on the stored candidate leptons and can only tighten the
h4l_select selection (candidate choice unchanged, a documented
approximation): electron MVA working point, isolation, SIP, the Z1 and Z2
minima, the lepton pT thresholds; the Z2 minima below 12 GeV use the SRZ4l
candidates (Z2 > 4 GeV), the others the SR candidates.  One-at-a-time scans
around the reference; every option records the yields, sigma_mu, sigma_mH
and J.  Writes production_v3/optimization/<select>/<label>/scan.json and the
plots scan_<variable>.png.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_select_io as io  # noqa: E402

Z_MASS = 91.1876
SIGNAL = ["GluGluToHToZZ_M125", "VBF_HToZZ_M125", "VHToZZ_M125"]
ZZ = ["ZZTo4L", "GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"]
FAKES = ["DYJetsToLL", "TTBar"]
ELECTRON_BITS = {"mvaFall17V2noIso_WPL": 0, "mvaFall17V2noIso_WP90": 1, "mvaFall17V2noIso_WP80": 2,
                 "mvaFall17V2Iso_WPL": 3, "mvaFall17V2Iso_WP90": 4, "mvaFall17V2Iso_WP80": 5}
REFERENCE = {"electron_id": "mvaFall17V2noIso_WPL", "iso": 0.35, "sip": 4.0, "z1_min": 40.0, "z2_min": 12.0,
             "lead_pt": 20.0, "sublead_pt": 10.0}
WINDOW = (105.0, 140.0)
BIN = 0.25
FR_PT = {13: [5, 7, 10, 15, 20, 30, 40, 50, 80, np.inf], 11: [7, 10, 15, 20, 30, 40, 50, 80, np.inf]}
FR_ETA_SPLIT = {13: 1.2, 11: 1.479}
CANDIDATE_BRANCHES = ["final_state", "m4l", "m4l_refit", "m4l_refit_err", "mz1", "mz2", "d_bkg_kin", "l_pdg", "l_pt", "l_eta",
                      "l_eta_sc", "l_iso", "l_sip", "l_flags", "l_pass", "region", "cr_type"]
ZL_BRANCHES = ["mz1", "met", "probe_pt", "probe_eta", "probe_eta_sc", "probe_iso", "probe_sip", "probe_flags", "probe_id",
               "probe_pdg"]


def lepton_pass(pdg, flags, iso, sip, opt) -> np.ndarray:
    """Selected-lepton test of the options on stored leptons (loose and AN muon ID from h4l_select)."""
    electron = np.abs(pdg) == 11
    bit = ELECTRON_BITS[opt["electron_id"]]
    id_ok = np.where(electron, (flags >> bit) & 1, 1).astype(bool)
    return id_ok & (iso < opt["iso"]) & (sip < opt["sip"])


def candidate_mask(c: dict, opt: dict) -> np.ndarray:
    passes = lepton_pass(c["l_pdg"], c["l_flags"], c["l_iso"], c["l_sip"], opt)
    all_pass = passes.all(axis=1)
    pts = np.sort(c["l_pt"], axis=1)
    kin = (pts[:, 3] > opt["lead_pt"]) & (pts[:, 2] > opt["sublead_pt"]) & (c["mz1"] > opt["z1_min"]) & (c["mz2"] > opt["z2_min"])
    return all_pass & kin


def fake_rates(zl: dict, opt: dict) -> dict:
    """Fake rates of the MC (DY, TTbar) Z + 1 loose lepton rows: denominator loose + SIP < sip."""
    base = (np.abs(zl["mz1"] - Z_MASS) < 7.0) & (zl["met"] < 25.0) & (zl["probe_sip"] < opt["sip"])
    electron = np.abs(zl["probe_pdg"]) == 11
    bit = ELECTRON_BITS[opt["electron_id"]]
    id_ok = np.where(electron, (zl["probe_flags"] >> bit) & 1, zl["probe_id"]).astype(bool)
    numerator = id_ok & (zl["probe_iso"] < opt["iso"])
    rates = {}
    for pdg in (13, 11):
        flavour = base & (np.abs(zl["probe_pdg"]) == pdg)
        eta = np.abs(zl["probe_eta"] if pdg == 13 else zl["probe_eta_sc"])
        table = np.zeros((len(FR_PT[pdg]) - 1, 2))
        for ie, region in enumerate((eta < FR_ETA_SPLIT[pdg], eta >= FR_ETA_SPLIT[pdg])):
            for ip in range(len(FR_PT[pdg]) - 1):
                sel = flavour & region & (zl["probe_pt"] >= FR_PT[pdg][ip]) & (zl["probe_pt"] < FR_PT[pdg][ip + 1])
                den = zl["w"][sel].sum()
                num = zl["w"][sel & numerator].sum()
                table[ip, ie] = min(max(num / den, 0.0), 0.95) if den > 0 else 0.0
        rates[pdg] = table
    return rates


def rate_of(rates: dict, pdg, pt, eta, eta_sc) -> np.ndarray:
    out = np.zeros(len(pdg))
    for flavour in (13, 11):
        sel = np.abs(pdg) == flavour
        edges = FR_PT[flavour]
        ip = np.clip(np.searchsorted(edges, pt[sel], side="right") - 1, 0, len(edges) - 2)
        ae = np.abs(eta[sel] if flavour == 13 else eta_sc[sel])
        ie = (ae >= FR_ETA_SPLIT[flavour]).astype(int)
        out[sel] = rates[flavour][ip, ie]
    return out


def zx_weights(cr: dict, rates: dict, opt: dict, region: int) -> tuple[np.ndarray, np.ndarray]:
    """OS-method weights of the control-region rows re-classified with the options."""
    passes = lepton_pass(cr["l_pdg"], cr["l_flags"], cr["l_iso"], cr["l_sip"], opt)
    sip_ok = (cr["l_sip"] < opt["sip"]).all(axis=1)
    pts = np.sort(cr["l_pt"], axis=1)
    kin = (pts[:, 3] > opt["lead_pt"]) & (pts[:, 2] > opt["sublead_pt"]) & (cr["mz1"] > opt["z1_min"]) & (cr["mz2"] > opt["z2_min"])
    base = (cr["region"] == region) & (cr["cr_type"] <= 1) & passes[:, 0] & passes[:, 1] & sip_ok & kin
    fails = (~passes[:, 2]).astype(int) + (~passes[:, 3]).astype(int)
    f3 = rate_of(rates, cr["l_pdg"][:, 2], cr["l_pt"][:, 2], cr["l_eta"][:, 2], cr["l_eta_sc"][:, 2])
    f4 = rate_of(rates, cr["l_pdg"][:, 3], cr["l_pt"][:, 3], cr["l_eta"][:, 3], cr["l_eta_sc"][:, 3])
    o3, o4 = f3 / (1 - f3), f4 / (1 - f4)
    weight = np.zeros(len(f3))
    one = base & (fails == 1)
    weight[one] = np.where(~passes[one, 2], o3[one], o4[one])
    two = base & (fails == 2)
    weight[two] = -o3[two] * o4[two]
    return weight * cr["w"], base & (fails >= 1)


def smooth_density(x: np.ndarray, w: np.ndarray, edges: np.ndarray, h: float) -> np.ndarray:
    """Weighted Gaussian-kernel density (events per GeV) at the given points."""
    if len(x) == 0:
        return np.zeros(len(edges))
    out = np.zeros(len(edges))
    for start in range(0, len(x), 20000):
        xs, ws = x[start:start + 20000], w[start:start + 20000]
        out += (ws[None, :] * np.exp(-0.5 * ((edges[:, None] - xs[None, :]) / h) ** 2)).sum(axis=1) / (h * np.sqrt(2 * np.pi))
    return out


class Model:
    """Binning (D quintiles, mass-error tertiles) frozen from the reference signal."""

    def __init__(self, signal: dict, mask: np.ndarray):
        s = {k: v[mask] for k, v in signal.items()}
        in_window = (s["m4l_refit"] > WINDOW[0]) & (s["m4l_refit"] < WINDOW[1])
        d = s["d_bkg_kin"][in_window]
        self.d_edges = np.quantile(d, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
        self.d_edges[0], self.d_edges[-1] = -1e-9, 1.0 + 1e-9
        self.err_edges = {}
        for fs in (0, 1, 2):
            rel = (s["m4l_refit_err"] / s["m4l_refit"])[in_window & (s["final_state"] == fs)]
            q = np.quantile(rel, [0, 1 / 3, 2 / 3, 1])
            q[0], q[-1] = 0.0, np.inf
            self.err_edges[fs] = q
        self.m_edges = np.arange(WINDOW[0], WINDOW[1] + 1e-9, BIN)

    def cells(self, rows: dict, fs: int):
        """Indices (D bin, error bin) of the rows of final state fs."""
        di = np.clip(np.searchsorted(self.d_edges, rows["d_bkg_kin"], side="right") - 1, 0, len(self.d_edges) - 2)
        rel = rows["m4l_refit_err"] / rows["m4l_refit"]
        ei = np.clip(np.searchsorted(self.err_edges[fs], rel, side="right") - 1, 0, 2)
        return di, ei

    def fisher(self, signal: dict, smask: np.ndarray, zz: dict, zmask: np.ndarray, zx: dict, zx_w: np.ndarray,
               zx_mask: np.ndarray) -> dict:
        info = np.zeros((2, 2))
        yields = {}
        centres = 0.5 * (self.m_edges[1:] + self.m_edges[:-1])
        for fs in (0, 1, 2):
            s = {k: v[smask & (signal["final_state"] == fs)] for k, v in signal.items()}
            b = {k: v[zmask & (zz["final_state"] == fs)] for k, v in zz.items()}
            x = {k: v[zx_mask & (zx["final_state"] == fs)] for k, v in zx.items()}
            xw = zx_w[zx_mask & (zx["final_state"] == fs)]
            sdi, sei = self.cells(s, fs)
            bdi, bei = self.cells(b, fs)
            xdi, xei = self.cells(x, fs)
            # Z+X: yield in the window from the 5 GeV shape, split by the cell composition of the rows.
            coarse = np.arange(70.0, 200.0 + 1e-9, 5.0)
            hist, _ = np.histogram(x["m4l_refit"], bins=coarse, weights=xw)
            density = np.clip(np.interp(centres, 0.5 * (coarse[1:] + coarse[:-1]), hist / 5.0), 0.0, None)
            in_w = (x["m4l_refit"] > WINDOW[0]) & (x["m4l_refit"] < WINDOW[1])
            total_zx = max(xw[in_w].sum(), 0.0)
            frac = np.zeros((5, 3))
            for di in range(5):
                for ei in range(3):
                    frac[di, ei] = max(xw[in_w & (xdi == di) & (xei == ei)].sum(), 0.0)
            frac = frac / frac.sum() if frac.sum() > 0 else np.full((5, 3), 1 / 15)
            y = {"signal": 0.0, "zz": 0.0, "zx": total_zx}
            for di in range(5):
                for ei in range(3):
                    sc = (sdi == di) & (sei == ei)
                    bc = (bdi == di) & (bei == ei)
                    s_dens = smooth_density(s["m4l_refit"][sc], s["w"][sc], self.m_edges, 0.3)
                    s_bins = 0.5 * (s_dens[1:] + s_dens[:-1]) * BIN
                    ds_bins = -(s_dens[1:] - s_dens[:-1])  # d/d m_H of the bin content (shift)
                    b_dens = smooth_density(b["m4l_refit"][bc], b["w"][bc], centres, 2.0) * BIN
                    zx_bins = density * BIN * frac[di, ei] * (total_zx / max(density.sum() * BIN, 1e-12))
                    nu = s_bins + b_dens + zx_bins
                    good = nu > 1e-12
                    grad = np.vstack([s_bins, ds_bins])[:, good]
                    info += grad @ (grad / nu[good]).T
                    y["signal"] += s_bins.sum()
                    y["zz"] += b_dens.sum()
            yields[["4mu", "4e", "2e2mu"][fs]] = y
        cov = np.linalg.inv(info)
        return {"sigma_mu": float(np.sqrt(cov[0, 0])), "sigma_mH": float(np.sqrt(cov[1, 1])),
                "correlation": float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])), "yields": yields}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True, help="h4l_select version with a full scan")
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    out_dir = PRODUCTION / "optimization" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    scan = io.load_scan(args.select)
    print(f"[opt] reading {args.select} (luminosity {io.lumi_fb(scan)} fb-1)", flush=True)
    trees = {}
    near_window = lambda a: (a["m4l_refit"] > WINDOW[0] - 10.0) & (a["m4l_refit"] < WINDOW[1] + 10.0)  # noqa: E731
    for tree in ("SR", "SRZ4l"):
        trees[("signal", tree)] = io.read(scan, SIGNAL, tree, CANDIDATE_BRANCHES, cut=near_window)
        trees[("zz", tree)] = io.read(scan, ZZ, tree, CANDIDATE_BRANCHES, cut=near_window)
    cr = io.read(scan, FAKES, "CR", CANDIDATE_BRANCHES)
    zl = io.read(scan, FAKES, "ZL", ZL_BRANCHES)
    print(f"[opt] rows: signal SR {len(trees[('signal', 'SR')]['w'])}, ZZ SR {len(trees[('zz', 'SR')]['w'])}, "
          f"fake CR {len(cr['w'])}, fake Z+l {len(zl['w'])}", flush=True)

    def evaluate(opt: dict, model: Model | None = None) -> tuple[dict, Model]:
        tree = "SR" if opt["z2_min"] >= 12.0 else "SRZ4l"
        region = 0 if tree == "SR" else 1
        signal, zz = trees[("signal", tree)], trees[("zz", tree)]
        smask = candidate_mask(signal, opt)
        zmask = candidate_mask(zz, opt)
        rates = fake_rates(zl, opt)
        zx_w, zx_mask = zx_weights(cr, rates, opt, region)
        if model is None:
            model = Model(signal, smask)
        result = model.fisher(signal, smask, zz, zmask, cr, zx_w, zx_mask)
        result["fake_rates"] = {str(k): v.tolist() for k, v in rates.items()}
        return result, model

    reference, model = evaluate(REFERENCE)
    print(f"[opt] reference: sigma_mu {reference['sigma_mu']:.4f}, sigma_mH {reference['sigma_mH']:.4f} GeV, "
          f"yields {json.dumps(reference['yields'])}", flush=True)
    scans = {"electron_id": list(ELECTRON_BITS), "iso": [0.15, 0.20, 0.25, 0.30, 0.35], "sip": [2.5, 3.0, 3.5, 4.0],
             "z2_min": [4.0, 6.0, 8.0, 10.0, 12.0, 14.0], "z1_min": [40.0, 45.0, 50.0],
             "lead_pt": [20.0, 22.0, 25.0], "sublead_pt": [10.0, 12.0, 15.0]}
    results = {"reference": {"options": REFERENCE, **reference}, "scans": {}}
    for variable, values in scans.items():
        rows = []
        for value in values:
            opt = dict(REFERENCE, **{variable: value})
            res, _ = evaluate(opt, model)
            objective = (res["sigma_mu"] / reference["sigma_mu"]) ** 2 + (res["sigma_mH"] / reference["sigma_mH"]) ** 2
            rows.append({"value": value, "objective": objective, **{k: res[k] for k in ("sigma_mu", "sigma_mH", "correlation", "yields")}})
            print(f"[opt] {variable} = {value}: J {objective:.4f}, sigma_mu {res['sigma_mu']:.4f}, "
                  f"sigma_mH {res['sigma_mH']:.4f}", flush=True)
        results["scans"][variable] = rows
    results["model"] = {"window": WINDOW, "bin": BIN, "d_edges": model.d_edges.tolist(),
                        "err_edges": {str(k): v.tolist() for k, v in model.err_edges.items()},
                        "signal": SIGNAL, "zz": ZZ, "fakes": FAKES, "fake_rate_bins": {str(k): list(map(float, v)) for k, v in FR_PT.items()}}
    results["select_scan"] = str(PRODUCTION / "h4l_select" / args.select / "scan.json")
    out_dir.mkdir(parents=True)
    (out_dir / "scan.json").write_text(json.dumps(results, indent=1, default=float) + "\n", encoding="utf-8")
    for variable, rows in results["scans"].items():
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        labels = [str(r["value"]) for r in rows]
        xs = np.arange(len(rows))
        for ax, key, title in zip(axes, ("objective", "sigma_mu", "sigma_mH"), ("J", "sigma(mu)", "sigma(m_H) [GeV]")):
            ax.plot(xs, [r[key] for r in rows], "o-")
            ax.set_xticks(xs, labels, rotation=30, fontsize=7)
            ax.set_title(title)
            ax.grid(alpha=0.3)
        fig.suptitle(f"stage 4b scan of {variable} (MC, {io.lumi_fb(scan)} fb-1; reference {REFERENCE[variable]})")
        fig.tight_layout()
        fig.savefig(out_dir / f"scan_{variable}.png", dpi=110)
        plt.close(fig)
    print(f"[opt] {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
