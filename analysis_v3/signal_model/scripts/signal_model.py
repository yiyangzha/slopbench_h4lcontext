"""Signal model of the H -> ZZ* -> 4l fits (stage 6).

    pixi run py -- analysis_v3/signal_model/scripts/signal_model.py --select v4 --yr4 v1 --label sm_v1 \
        [--categories analysis_v3/reconstruction/config/categories_ul16_v2.json]

Inputs: the final selection scan (signal-region rows of the ggH, VBF and VH samples at
m_H = 125 GeV), the YR4 table (yr4_xsbr.py) and the category configuration.
Per final state (4mu, 4e, 2e2mu), with the Z1-refitted mass m4l_refit and without
the refit (m4l), in the window [lo, hi] (default 105-140 GeV):
  * resonant shape: a double-sided Crystal Ball (DCB) fitted (weighted, unbinned) to
    the ggH + VBF events;
  * the per-event width model of the 3D fit: the DCB conditional on the relative
    per-event mass uncertainty D_mass = dm_i / m4l_i (the third observable; dm_i =
    m4l_refit_err or m4l_err), with width = s x D_mass x 125 GeV (the generated m_H), s
    fitted (it absorbs the MC calibration of the per-event errors); the likelihood uses
    s x D_mass x m_H, a density in m4l normalized at fixed D_mass;
  * VH: the resonant DCB of the final state plus a non-resonant Landau component
    (leptons not from the Higgs decay), fraction and Landau fitted to the VH events;
  * m_H dependence by rest-frame scaling of the 4l system (user decision): every
    signal event is scaled by k = m_H / 125 (m4l, the Z masses, the lepton pT, the
    per-event error), so the resonant DCB keeps its tails with mean(m_H) = k mean_125
    and width(m_H) = k width_125; the acceptance x efficiency ratio
    A(m_H)/A(125) per production mode, final state and category comes from
    re-applying the kinematic selection (lepton pT 5/7 GeV, 20/10 GeV leading,
    40 < m_Z1 < 120 GeV, z2_low < m_Z2 < 120 GeV, the window) to the scaled events;
    the non-resonant VH part does not move with m_H; the alternative "shift"
    morphing (m4l + (m_H - 125) at constant width) and its acceptance are recorded
    for the morphing systematic;
  * normalization: sigma_eff(mode, m_H) x L x A x eff(mode, fs, cat)(m_H), with
    sigma_eff from YR4 (tabulated 120-130 GeV; outside, log-quadratic extrapolation of
    the table, documented), A x eff(125) = the weighted signal-region yield / (sigma_eff
    L) from the MC (scale factors applied downstream by the model builder).
Writes production_v3/signal_model/<select>/<label>/signal_model.json and plots/.
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
import h4l_categories as cats  # noqa: E402
import h4l_select_io as io  # noqa: E402
import h4l_shapes as shapes  # noqa: E402
from iminuit import Minuit  # noqa: E402

MODES = {"GluGluToHToZZ_M125": "ggH", "VBF_HToZZ_M125": "VBF", "VHToZZ_M125": "VH"}
YR4_MODES = {"ggH": "GluGluToHToZZ", "VBF": "VBF_HToZZ", "VH": "VHToZZ"}
FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}
BRANCHES = ["final_state", "m4l", "m4l_err", "m4l_refit", "m4l_refit_err", "mz1", "mz2", "l_pdg", "l_pt", "d_bkg_kin"] + \
           cats.CATEGORY_BRANCHES
MH_GRID = [110.0, 115.0, 118.0, 120.0, 122.0, 124.0, 125.0, 126.0, 128.0, 130.0, 132.0, 135.0, 140.0]
MH_MC = 125.0


def sigma_eff_function(yr4: dict, mode: str):
    """sigma_eff(mode, m_H) from the YR4 table; log-quadratic extrapolation outside it."""
    rows = yr4["sigma_eff"][YR4_MODES[mode]]
    masses = np.array([r["mH"] for r in rows])
    values = np.array([r["sigma_eff_pb"] for r in rows])
    coefficients = np.polyfit(masses - 125.0, np.log(values), 2)

    def f(m):
        m = float(m)
        if masses[0] <= m <= masses[-1]:
            return float(np.interp(m, masses, values))
        return float(np.exp(np.polyval(coefficients, m - 125.0)))
    return f, {"table_range": [float(masses[0]), float(masses[-1])], "log_quadratic": coefficients.tolist()}


def passes_kinematics(rows: dict, k: float, z2_low: float, lo: float, hi: float, variable: str) -> np.ndarray:
    pt = rows["l_pt"] * k
    pdg = np.abs(rows["l_pdg"])
    lepton_ok = np.where(pdg == 13, pt > 5.0, pt > 7.0).all(axis=1)
    ordered = np.sort(pt, axis=1)
    mz1, mz2, mass = rows["mz1"] * k, rows["mz2"] * k, rows[variable] * k
    return lepton_ok & (ordered[:, 3] > 20.0) & (ordered[:, 2] > 10.0) & (mz1 > 40.0) & (mz1 < 120.0) & (mz2 > z2_low) & \
        (mz2 < 120.0) & (mass > lo) & (mass < hi)


def fit_vh(m4l, w, dcb: dict, lo: float, hi: float) -> dict:
    """VH = f_res DCB (shape fixed from ggH + VBF) + (1 - f_res) Landau."""
    sel = (m4l > lo) & (m4l < hi)
    x, ww = m4l[sel], w[sel]
    v = dcb["values"]
    res = shapes.dcb_pdf(x, v["mean"], v["width"], v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"], lo, hi)

    def nll(f_res, mpv, width):
        dens = f_res * res + (1 - f_res) * shapes.landau_pdf(x, mpv, width, lo, hi)
        return -np.sum(ww * np.log(np.maximum(dens, 1e-300))) / np.sum(ww) * len(x)

    m = Minuit(nll, f_res=0.8, mpv=120.0, width=15.0)
    m.errordef = Minuit.LIKELIHOOD
    m.limits["f_res"] = (0.0, 1.0)
    m.limits["mpv"] = (80.0, 200.0)
    m.limits["width"] = (1.0, 100.0)
    m.migrad()
    m.hesse()
    return {"f_res": float(m.values["f_res"]), "f_res_error": float(m.errors["f_res"]), "landau_mpv": float(m.values["mpv"]),
            "landau_width": float(m.values["width"]), "valid": bool(m.valid), "n_events": int(len(x))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--categories", type=Path, default=REPO / "analysis_v3/reconstruction/config/categories_ul16_v2.json")
    parser.add_argument("--window", type=float, nargs=2, default=[105.0, 140.0])
    args = parser.parse_args()
    out_dir = PRODUCTION / "signal_model" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    lo, hi = args.window
    scan = io.load_scan(args.select)
    lumi = io.lumi_fb(scan)
    plan = json.loads(Path(scan["plan"]).read_text(encoding="utf-8"))
    regions = {r["name"]: r["z2_low"] for r in plan["tasks"][0]["selection_config"]["candidate"]["regions"]}
    z2_low = regions["SR"]
    categories = json.loads((args.categories if args.categories.is_absolute() else REPO / args.categories).read_text(encoding="utf-8"))
    yr4 = json.loads((PRODUCTION / "signal_model" / "yr4" / args.yr4 / "yr4.json").read_text(encoding="utf-8"))
    sigma_eff = {}
    extrapolation = {}
    for mode in MODES.values():
        sigma_eff[mode], extrapolation[mode] = sigma_eff_function(yr4, mode)
    rows = {}
    for sample, mode in MODES.items():
        r = io.read(scan, [sample], "SR", BRANCHES)
        r["category"] = cats.category_index(r, categories)
        rows[mode] = r
        print(f"[signal] {mode}: {len(r['w'])} SR rows, weighted {np.sum(r['w']):.3f} events (all m4l)", flush=True)
    model = {"schema": "h4l_v3_signal_model/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
             "select_plan_sha256": scan["plan_sha256"], "yr4": str(PRODUCTION / "signal_model" / "yr4" / args.yr4 / "yr4.json"),
             "categories": categories, "window": [lo, hi], "luminosity_fb": lumi, "z2_low": z2_low,
             "sigma_eff_extrapolation": extrapolation, "final_states": {}}
    plots = out_dir / "plots"
    for fs, fs_name in FINAL_STATES.items():
        entry = {}
        for variable, error in (("m4l_refit", "m4l_refit_err"), ("m4l", "m4l_err")):
            gv = {k: np.concatenate([rows[m][k][rows[m]["final_state"] == fs] for m in ("ggH", "VBF")]) for k in ("w", variable, error)}
            dcb = shapes.fit_dcb(gv[variable], gv["w"], lo, hi)
            dcb_event = shapes.fit_dcb(gv[variable], gv["w"], lo, hi, sigma_events=gv[error] / gv[variable] * MH_MC,
                                       start={k: dcb["values"][k] for k in ("mean", "alpha_l", "n_l", "alpha_r", "n_r")})
            vh = rows["VH"]
            in_fs = vh["final_state"] == fs
            vh_fit = fit_vh(vh[variable][in_fs], vh["w"][in_fs], dcb, lo, hi)
            entry[variable] = {"dcb": dcb, "dcb_per_event_width": dcb_event, "vh_nonresonant": vh_fit, "error_variable": error}
        # Yields and acceptance ratios per mode and category (refitted mass defines the window).
        yields = {}
        for mode, r in rows.items():
            in_fs = r["final_state"] == fs
            per_cat = {}
            for c, name in enumerate(categories["order"]):
                sel = in_fs & (r["category"] == c)
                base = passes_kinematics({k: r[k][sel] for k in ("l_pt", "l_pdg", "mz1", "mz2", "m4l_refit")}, 1.0, z2_low, lo, hi,
                                         "m4l_refit")
                y125 = float(np.sum(r["w"][sel][base]))
                ratio_scale, ratio_shift = [], []
                for mh in MH_GRID:
                    k = mh / 125.0
                    sub = {kk: r[kk][sel] for kk in ("l_pt", "l_pdg", "mz1", "mz2", "m4l_refit")}
                    acc = passes_kinematics(sub, k, z2_low, lo, hi, "m4l_refit")
                    ratio_scale.append(float(np.sum(r["w"][sel][acc]) / y125) if y125 > 0 else 1.0)
                    shifted = dict(sub, m4l_refit=sub["m4l_refit"] + (mh - 125.0))
                    acc_shift = passes_kinematics(shifted, 1.0, z2_low, lo, hi, "m4l_refit")
                    ratio_shift.append(float(np.sum(r["w"][sel][acc_shift]) / y125) if y125 > 0 else 1.0)
                n_raw = int(np.sum(sel & ((r["m4l_refit"] > lo) & (r["m4l_refit"] < hi))))
                per_cat[name] = {"yield_125": y125, "raw_events": n_raw,
                                 "acceptance_efficiency_125": y125 / (sigma_eff[mode](125.0) * 1000.0 * lumi),
                                 "acceptance_ratio_scale": ratio_scale, "acceptance_ratio_shift": ratio_shift}
            yields[mode] = per_cat
        entry["yields"] = yields
        model["final_states"][fs_name] = entry
        print(f"[signal] {fs_name}: DCB(refit) mean {entry['m4l_refit']['dcb']['values']['mean']:.3f} width "
              f"{entry['m4l_refit']['dcb']['values']['width']:.3f}; per-event scale {entry['m4l_refit']['dcb_per_event_width']['values']['width']:.3f}; "
              f"VH f_res {entry['m4l_refit']['vh_nonresonant']['f_res']:.3f}; yield ggH/VBF/VH "
              f"{sum(v['yield_125'] for v in yields['ggH'].values()):.2f}/{sum(v['yield_125'] for v in yields['VBF'].values()):.2f}/"
              f"{sum(v['yield_125'] for v in yields['VH'].values()):.2f}", flush=True)
    model["mh_grid"] = MH_GRID
    model["sigma_eff_pb"] = {mode: {str(m): sigma_eff[mode](m) for m in MH_GRID} for mode in MODES.values()}
    model["morphing"] = ("rest-frame scaling by k = m_H/125: resonant DCB mean and width scale by k, tails fixed; acceptance "
                         "ratios from the scaled events; non-resonant VH fixed; 'shift' morphing recorded for the systematic")
    out_dir.mkdir(parents=True)
    plots.mkdir()
    for fs, fs_name in FINAL_STATES.items():
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        for ax, variable in zip(axes, ("m4l_refit", "m4l")):
            gv = {k: np.concatenate([rows[m][k][rows[m]["final_state"] == fs] for m in ("ggH", "VBF")]) for k in ("w", variable)}
            edges = np.arange(lo, hi + 1e-9, 0.5)
            h, _ = np.histogram(gv[variable], bins=edges, weights=gv["w"])
            centres = 0.5 * (edges[1:] + edges[:-1])
            v = model["final_states"][fs_name][variable]["dcb"]["values"]
            dens = shapes.dcb_pdf(centres, v["mean"], v["width"], v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"], lo, hi)
            total = np.sum(gv["w"][(gv[variable] > lo) & (gv[variable] < hi)])
            ax.errorbar(centres, h, np.sqrt(np.histogram(gv[variable], bins=edges, weights=gv["w"] ** 2)[0]), fmt="o", ms=2, color="black",
                        label="ggH + VBF MC")
            ax.plot(centres, dens * total * 0.5, color="#3f90da", label="DCB fit")
            ax.set_yscale("log")
            ax.set_ylim(max(1e-4, 0.5 * np.min(h[h > 0])) if np.any(h > 0) else 1e-4, 3 * np.max(h))
            ax.set_xlabel(f"{variable} [GeV]")
            ax.set_title(f"{fs_name}: mean {v['mean']:.2f}, width {v['width']:.2f}, aL {v['alpha_l']:.2f}, nL {v['n_l']:.1f}, "
                         f"aR {v['alpha_r']:.2f}, nR {v['n_r']:.1f}", fontsize=7)
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(plots / f"signal_dcb_{fs_name}.png", dpi=120)
        plt.close(fig)
    (out_dir / "signal_model.json").write_text(json.dumps(model, indent=1) + "\n", encoding="utf-8")
    print(f"[signal] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
