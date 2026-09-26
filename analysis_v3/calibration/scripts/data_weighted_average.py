"""Data-weighted global averages of the per-flavour lepton calibration and efficiency SF (user decision 2026-09-25).

    pixi run py -- analysis_v3/calibration/scripts/data_weighted_average.py --calibration v4/nominal \
        --sf production_v3/tnp/v2/v7/fits_i1c/sf.json \
        --systematics analysis_v3/inference/config/systematics_ul16_v3.json --label average_v1 [--workers 8]

RESULT.json reports one value per flavour for scale_shift, smear and sel_eff: the average of the per-lepton
calibration (the factorized model ln(1 + s) = a(|eta|) + b_R(pT), r^2 = c(|eta|) + d_R(pT) of the payload) and of
the per-bin tag-and-probe SF of the full single-lepton selection over the data leptons of the Z -> ll calibration
sample (user choice): both legs of the data control pairs of the calibration extract passing the calibration leg
selection (muons AN tight, electrons mvaFall17V2noIso WP90 as in calibration v4, FSR-subtracted isolation < 0.35,
SIP < 4, corrected pT above 3 / 5 GeV) with 60 < m_ll < 120 GeV after the data correction.  Uncertainties:
  * scale_shift: the solution covariance of the last iteration propagated to the average (gradient exp(u) times the
    model design per leg), in quadrature with the iteration and the larger closure component of the systematics
    configuration (the terms of the per-flavour scale uncertainty of the fits);
  * smear: the covariance propagated through r = sqrt(v) (gradient 1 / (2 r) per leg), in quadrature with the
    closure variance deviation / (2 r_avg);
  * sel_eff: the statistical parts of the bins in quadrature with their data-lepton weights, plus the fit-model
    parts added linearly (coherent within a flavour).
Writes production_v3/results/<select>/<label>/average.json (never overwritten).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import uproot

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_sf  # noqa: E402

BRANCHES = ["flavour", "mass", "weight"] + [f"l{i}_{v}" for i in (1, 2) for v in ("pt", "eta", "eta_sc", "iso_fsr", "sip", "flags")]
FLAVOURS = {"muon": (13, "mm"), "electron": (11, "ee")}


def muon_an_tight(flags):
    flags = flags.astype(np.int64)
    glob, tracker, pf, stations = flags & 1, flags >> 1 & 1, flags >> 3 & 1, flags >> 10 & 1
    return (pf == 1) & ((glob == 1) | ((tracker == 1) & (stations == 1)))


def electron_wp90(flags):
    return (flags.astype(np.int64) >> 1 & 1) == 1


class Factorized:
    """The factorized payload of one flavour (analysis_v3/common/include/h4l/calibration.h)."""

    def __init__(self, node: dict):
        model = node["model"]
        self.eta = np.array(model["eta_edges"])
        self.regions = np.array(model["region_edges"])
        self.pt = np.array(model["pt_edges"])
        self.ref = model["reference_pt_bin"]
        self.eta_sc = model["eta_variable"] == "abs_eta_sc"
        self.a, self.c = np.array(node["a"]), np.array(node["c"])
        self.b, self.d = np.array(node["b"]), np.array(node["d"])
        self.nodes = np.array(node["node_pt"])
        self.n_eta, self.n_regions, self.n_pt = len(self.eta), len(self.regions), len(self.pt)
        self.b_index = {}
        for r in range(self.n_regions):
            for p in range(self.n_pt):
                if p != self.ref:
                    self.b_index[(r, p)] = self.n_eta + len(self.b_index)
        self.n_par = self.n_eta + len(self.b_index)

    def bins(self, pt, abs_eta):
        e = np.searchsorted(self.eta, abs_eta, side="right") - 1
        r = np.searchsorted(self.regions, abs_eta, side="right") - 1
        return e, r

    def interpolation(self, pt, r):
        """Per leg: the two pT nodes of its region and their weights (constant beyond the outermost nodes)."""
        i0 = np.zeros(len(pt), dtype=int)
        w1 = np.zeros(len(pt))
        for region in range(self.n_regions):
            sel = r == region
            nodes = self.nodes[region]
            x = np.clip(pt[sel], nodes[0], nodes[-1])
            i = np.clip(np.searchsorted(nodes, x, side="right") - 1, 0, len(nodes) - 2)
            i0[sel] = i
            w1[sel] = (x - nodes[i]) / (nodes[i + 1] - nodes[i])
        return i0, w1

    def evaluate(self, pt, abs_eta, eta_values, pt_values):
        e, r = self.bins(pt, abs_eta)
        i0, w1 = self.interpolation(pt, r)
        return eta_values[e] + (1 - w1) * pt_values[r, i0] + w1 * pt_values[r, i0 + 1]

    def design_mean(self, pt, abs_eta, factor):
        """Mean over legs of factor x d(u or v)/d(parameters) in the solution's parameter order."""
        e, r = self.bins(pt, abs_eta)
        i0, w1 = self.interpolation(pt, r)
        g = np.zeros(self.n_par)
        np.add.at(g, e, factor)
        for region in range(self.n_regions):
            for p in range(self.n_pt):
                if p == self.ref:
                    continue
                k = self.b_index[(region, p)]
                sel_lo = (r == region) & (i0 == p)
                sel_hi = (r == region) & (i0 + 1 == p)
                g[k] += np.sum(factor[sel_lo] * (1 - w1[sel_lo])) + np.sum(factor[sel_hi] * w1[sel_hi])
        return g / len(pt)


def read_legs(path: str) -> dict:
    with uproot.open(path) as f:
        a = f["CalibPairs"].arrays(BRANCHES, library="np")
    return a


def propagate(g: np.ndarray, solution: dict) -> float:
    cov = solution["covariance"]
    idx = np.array(cov["parameters"], dtype=int)
    matrix = np.array(cov["matrix"])
    gd = g[idx]
    return float(math.sqrt(max(gd @ matrix @ gd, 0.0)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--calibration", required=True, help="<version>/<run> under production_v3/calibration")
    parser.add_argument("--sf", type=Path, required=True)
    parser.add_argument("--systematics", type=Path, required=True)
    parser.add_argument("--select", default="v5")
    parser.add_argument("--label", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    sf_path = args.sf if args.sf.is_absolute() else REPO / args.sf
    syst_path = args.systematics if args.systematics.is_absolute() else REPO / args.systematics
    out = PRODUCTION / "results" / args.select / args.label
    if (out / "average.json").exists():
        raise SystemExit(f"{out / 'average.json'} exists")
    run_dir = PRODUCTION / "calibration" / args.calibration
    payload = json.loads((run_dir / "payload.json").read_text(encoding="utf-8"))
    solution = json.loads(Path(payload["solution"]).read_text(encoding="utf-8"))
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    config = json.loads(Path(run["config"]).read_text(encoding="utf-8"))
    scan = json.loads(Path(run["extract_scan"]).read_text(encoding="utf-8"))
    syst = json.loads(syst_path.read_text(encoding="utf-8"))
    sf = h4l_sf.LeptonSF(sf_path)
    files = sorted(o["path"] for o in scan["outputs"] if o["kind"] == "data")
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        parts = list(pool.map(read_legs, files))
    legs = {k: np.concatenate([p[k] for p in parts]) for k in BRANCHES}
    report = {"schema": "h4l_v3_data_weighted_average/1", "calibration": str(run_dir), "sf": str(sf_path),
              "systematics": str(syst_path), "data_files": files,
              "definition": "average over the data leptons of the Z -> ll calibration sample (both legs of the data control "
                            "pairs of the calibration extract passing the calibration leg selection, 60 < m_ll < 120 GeV after "
                            "the data correction)", "flavours": {}}
    for flavour, (pdg, key) in FLAVOURS.items():
        model = Factorized(payload["flavours"][flavour])
        sel_cfg = config["selection"][flavour]
        pair = legs["flavour"] == pdg
        leg_pt, leg_eta, leg_u, passing = [], [], [], None
        cols = {}
        for i in (1, 2):
            pt = legs[f"l{i}_pt"][pair].astype(float)
            abs_eta = np.abs(legs[f"l{i}_eta_sc" if model.eta_sc else f"l{i}_eta"][pair].astype(float))
            abs_eta_lepton = np.abs(legs[f"l{i}_eta"][pair].astype(float))
            e, r = model.bins(pt, abs_eta)
            inside = (e >= 0) & (r >= 0) & (pt >= model.pt[0])
            u = np.where(inside, model.evaluate(pt, np.where(inside, abs_eta, 0.0), model.a, model.b), 0.0)
            ident = muon_an_tight(legs[f"l{i}_flags"][pair]) if pdg == 13 else electron_wp90(legs[f"l{i}_flags"][pair])
            ok = (inside & ident & (legs[f"l{i}_iso_fsr"][pair] < sel_cfg["max_iso_fsr"]) & (legs[f"l{i}_sip"][pair] < sel_cfg["max_sip"])
                  & (pt * np.exp(-u) >= sel_cfg["min_pt"]) & (abs_eta_lepton < sel_cfg["max_abs_eta"]))
            cols[i] = (pt, abs_eta, u, legs[f"l{i}_eta"][pair].astype(float), legs[f"l{i}_eta_sc"][pair].astype(float))
            passing = ok if passing is None else passing & ok
        mass = legs["mass"][pair].astype(float) * np.exp(-0.5 * (cols[1][2] + cols[2][2]))
        passing &= (mass > 60.0) & (mass < 120.0)
        pt = np.concatenate([cols[i][0][passing] for i in (1, 2)])
        abs_eta = np.concatenate([cols[i][1][passing] for i in (1, 2)])
        u = np.concatenate([cols[i][2][passing] for i in (1, 2)])
        eta = np.concatenate([cols[i][3][passing] for i in (1, 2)])
        eta_sc = np.concatenate([cols[i][4][passing] for i in (1, 2)])
        n = len(pt)
        # scale shift s = exp(u) - 1
        s = np.exp(u) - 1.0
        scale_stat = propagate(model.design_mean(pt, abs_eta, np.exp(u)), solution[key]["scale"])
        comps = syst["details"]["per_flavour"]["mu" if pdg == 13 else "e"]
        closure = max(comps["scale"]["closure_bias_rms"].values())
        scale_unc = math.sqrt(scale_stat ** 2 + comps["scale"]["iteration"] ** 2 + closure ** 2)
        # smear r = sqrt(max(v, 0))
        v = model.evaluate(pt, abs_eta, model.c, model.d)
        r_leg = np.sqrt(np.clip(v, 0.0, None))
        r_avg = float(np.mean(r_leg))
        factor = np.where(r_leg > 0, 0.5 / np.where(r_leg > 0, r_leg, 1.0), 0.0)
        smear_stat = propagate(model.design_mean(pt, abs_eta, factor), solution[key]["smear"])
        smear_closure = comps["resolution"]["smear_variance_closure_max_dev"] / (2.0 * r_avg) if r_avg > 0 else 0.0
        smear_unc = math.hypot(smear_stat, smear_closure)
        # tag-and-probe SF of the full selection on the corrected pT (the tag-and-probe bins use calibrated momenta)
        pt_corr = pt * np.exp(-u)
        pdg_arr = np.full(n, pdg)
        value, _ = sf.lepton(pdg_arr, pt_corr, eta, eta_sc)
        k = sf.bins(pdg_arr, pt_corr, eta, eta_sc)
        stat, fit_model = sf.parts[pdg]
        counts = np.bincount(k[k >= 0], minlength=len(stat)).astype(float)
        wk = counts / counts.sum()
        sf_stat = float(math.sqrt(np.sum((wk * stat) ** 2)))
        sf_model = float(np.sum(wk * fit_model))
        sf_value = float(np.mean(value))
        report["flavours"][flavour] = {
            "n_legs": int(n),
            "scale_shift": {"value": float(np.mean(s)), "unc": scale_unc,
                            "components": {"stat": scale_stat, "iteration": comps["scale"]["iteration"], "closure": closure}},
            "smear": {"value": r_avg, "unc": smear_unc, "components": {"stat": smear_stat, "closure": smear_closure}},
            "sel_eff": {"value": sf_value, "unc": math.hypot(sf_stat, sf_model),
                        "components": {"stat": sf_stat, "fit_model": sf_model}},
        }
        print(f"[average] {flavour}: {n} legs; scale {np.mean(s):+.5f} +- {scale_unc:.5f} (stat {scale_stat:.5f}); smear {r_avg:.5f} +- "
              f"{smear_unc:.5f} (stat {smear_stat:.5f}); SF {sf_value:.5f} +- {math.hypot(sf_stat, sf_model):.5f}", flush=True)
    out.mkdir(parents=True, exist_ok=True)
    (out / "average.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[average] {out / 'average.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
