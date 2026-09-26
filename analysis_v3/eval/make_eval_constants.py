"""The method constants of the evaluation submission (my_analysis/h4l_eval/data/constants.json).

    pixi run py -- analysis_v3/eval/make_eval_constants.py

The submission re-measures every dataset-dependent quantity (lepton calibration, tag-and-probe efficiencies, fake
rates, Z+X yields, MC normalization, templates) on each dataset it runs on.  What it takes from the main analysis are
physics or method constants that do not depend on the dataset (user rule 2026-09-25: e.g. shapes that barely change
between datasets need no extensive search):
  * the true Z1 line shape of the Z1 kinematic refit and the FSR-photon pT resolution (signal MC generator level,
    production_v3/h4l_reco/refit/v3/refit_inputs.json);
  * the per-lepton momentum-error scale factors lambda of the MC (production_v3/calibration/v4/nominal/lambda_v2_r2),
    applied to data and MC alike (the data/MC resolution difference is carried by the measured smear);
  * the YR4 13 TeV cross sections (ggF N3LO, VBF, WH + ZH) and BR(H -> 4l) against m_H over 105-145 GeV: the SM
    sheets on their fine 120-130 GeV grid (the tables of production_v3/signal_model/yr4/v1), extended outside with the
    BSM sheet cross sections and BR(H -> ZZ) from the partial widths of the width sheet, each scaled to join the SM
    values at 120 and 130 GeV (read from the YR4 spreadsheet, the one authorized file under results/);
  * the Z+X m4l shapes (Landau + exponential per final state) and the relative MC-closure systematic of the Z+X
    method (production_v3/backgrounds/v5/zx_v2/zx.json), with the fake-rate binning;
  * the method systematics of the lepton calibration from its closures (analysis_v3/inference/config/
    systematics_ul16_v3.json: scale closure bias and iteration, smear variance closure deviation, per flavour).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
OUT = REPO / "my_analysis/h4l_eval/data/constants.json"
sys.path.insert(0, str(REPO / "analysis_v3/signal_model/scripts"))
import yr4_xsbr  # noqa: E402


def rounded(values, digits=7):
    return [float(f"{v:.{digits}g}") for v in values]


def yr4_extended() -> dict:
    """sigma(mode, m) [pb] for ggH, VBF, VH and BR4l(m) on a grid covering 105-145 GeV."""
    yr4 = json.loads((PRODUCTION / "signal_model/yr4/v1/yr4.json").read_text(encoding="utf-8"))
    tables = {k: np.array(v) for k, v in yr4["tables"].items()}
    rows = yr4_xsbr.read_sheets(yr4_xsbr.SPREADSHEET, {"YR4 BSM 13TeV", "YR4 BSM Width"})
    wide = {"ggF": yr4_xsbr.column_pairs(rows["YR4 BSM 13TeV"], 1, 2), "VBF": yr4_xsbr.column_pairs(rows["YR4 BSM 13TeV"], 17, 18),
            "WH": yr4_xsbr.column_pairs(rows["YR4 BSM 13TeV"], 26, 27), "ZH": yr4_xsbr.column_pairs(rows["YR4 BSM 13TeV"], 37, 38)}
    width = rows["YR4 BSM Width"]
    # Total width = the sum of the partial widths (fermions: mass column 1; bosons: mass column 21).
    fermions = [2, 5, 8, 11, 14, 17]
    bosons = [22, 25, 28, 31, 34]
    total, zz = {}, {}
    for cols, mcol in ((fermions, 1), (bosons, 21)):
        for c in cols:
            for m, v in yr4_xsbr.column_pairs(width, mcol, c):
                total[m] = total.get(m, 0.0) + v
                if c == 34:
                    zz[m] = v
    br_zz_wide = sorted((m, zz[m] / total[m]) for m in zz if m in total and total[m] > 0)
    grid = np.round(np.arange(105.0, 145.0001, 0.5), 3)
    fine = sorted(set(np.round(tables["ggF"][:, 0], 3)) | set(grid))
    grid = np.array([m for m in fine if 105.0 <= m <= 145.0])

    def joined(sm: np.ndarray, wide_pairs: list) -> np.ndarray:
        """SM values on [120, 130]; outside, the wide table (log-linear in m) scaled to join at 120 and 130 GeV."""
        wm = np.array([p[0] for p in wide_pairs])
        wv = np.log(np.array([p[1] for p in wide_pairs]))
        out = np.empty(len(grid))
        lo_scale = np.interp(120.0, sm[:, 0], sm[:, 1]) / np.exp(np.interp(120.0, wm, wv))
        hi_scale = np.interp(130.0, sm[:, 0], sm[:, 1]) / np.exp(np.interp(130.0, wm, wv))
        for i, m in enumerate(grid):
            if 120.0 <= m <= 130.0:
                out[i] = np.interp(m, sm[:, 0], sm[:, 1])
            else:
                out[i] = np.exp(np.interp(m, wm, wv)) * (lo_scale if m < 120.0 else hi_scale)
        return out

    sigma = {"ggH": joined(tables["ggF"], wide["ggF"]), "VBF": joined(tables["VBF"], wide["VBF"])}
    wh, zh = joined(tables["WH"], wide["WH"]), joined(tables["ZH"], wide["ZH"])
    sigma["VH"] = wh + zh
    br4l = joined(tables["BR4l"], br_zz_wide)
    return {"source": "YR4 13 TeV (Higgs_XSBR_YR4_update.xlsx): SM sheets 120-130 GeV, BSM sheet cross sections and "
                      "BR(H->ZZ) = Gamma_ZZ / sum of partial widths outside, joined at 120 and 130 GeV; BR4l is "
                      "BR(H -> 4l, l = e, mu)",
            "mass": rounded(grid), "sigma_pb": {k: rounded(v) for k, v in sigma.items()}, "br4l": rounded(br4l)}


def main() -> int:
    refit = json.loads((PRODUCTION / "h4l_reco/refit/v3/refit_inputs.json").read_text(encoding="utf-8"))
    z1 = refit["z1_lineshape"]
    edges = np.array(z1["edges"])
    density = np.array(z1["density"])
    log_density = np.log(np.clip(density, 1e-300, None))
    lam = json.loads((PRODUCTION / "calibration/v4/nominal/lambda_v2_r2/lambda.json").read_text(encoding="utf-8"))
    zx = json.loads((PRODUCTION / "backgrounds/v5/zx_v2/zx.json").read_text(encoding="utf-8"))
    syst = json.loads((REPO / "analysis_v3/inference/config/systematics_ul16_v3.json").read_text(encoding="utf-8"))
    lambda_mc = {}
    for flavour, regions in lam["flavours"].items():
        lambda_mc[flavour] = [{"abs_eta": r["region"]["abs_eta"], "rel_err": r["region"].get("rel_err", [0.0, None]),
                               "lambda": r["mc"]["lambda"]} for r in regions["regions"]]
    calib_syst = {}
    for flavour, key in (("muon", "mu"), ("electron", "e")):
        d = syst["details"]["per_flavour"][key]
        calib_syst[flavour] = {"scale_closure": max(d["scale"]["closure_bias_rms"].values()), "scale_iteration": d["scale"]["iteration"],
                               "smear_variance_closure": d["resolution"]["smear_variance_closure_max_dev"]}
    constants = {
        "schema": "h4l_eval_constants/1",
        "generated_by": "analysis_v3/eval/make_eval_constants.py",
        "z1_lineshape": {"lo": float(edges[0]), "step": float(edges[1] - edges[0]), "log_density": rounded(log_density, 8),
                         "definition": z1["definition"]},
        "fsr_photon_resolution": {"form": refit["fsr_photon_resolution"]["form"], "a": refit["fsr_photon_resolution"]["a"],
                                  "b": refit["fsr_photon_resolution"]["b"]},
        "lambda_mc": lambda_mc,
        "yr4": yr4_extended(),
        "zx": {"shapes": {fs: {k: zx["shapes"][fs][k] for k in ("form", "range", "mpv", "width", "exp_fraction", "exp_slope")}
                          for fs in ("4mu", "4e", "2e2mu")},
               "relative_systematic": zx["relative_systematic"],
               "fake_rate_pt_edges": {k: [x if x is not None else None for x in v] for k, v in zx["fake_rate_bins"].items()},
               "eta_split": zx["eta_split"],
               "source": "production_v3/backgrounds/v5/zx_v2/zx.json (main analysis, 20 fb^-1 pseudo-data and MC closure)"},
        "calibration_systematics": calib_syst,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(constants, indent=1) + "\n", encoding="utf-8")
    print(f"[constants] {OUT} ({OUT.stat().st_size / 1024:.1f} kB)")
    y = constants["yr4"]
    for m in (110.0, 120.0, 125.0, 130.0, 140.0):
        i = y["mass"].index(m)
        print(f"  m_H {m}: ggH {y['sigma_pb']['ggH'][i]:.3f} VBF {y['sigma_pb']['VBF'][i]:.4f} VH {y['sigma_pb']['VH'][i]:.4f} pb, "
              f"BR4l {y['br4l'][i]:.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
