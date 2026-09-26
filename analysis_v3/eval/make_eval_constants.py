"""The method constants of the evaluation submission (my_analysis/h4l_eval/data/constants.json).

    pixi run py -- analysis_v3/eval/make_eval_constants.py

The submission re-measures every dataset-dependent quantity (lepton calibration, tag-and-probe efficiencies, fake
rates, Z+X yields, MC normalization, templates) on each dataset it runs on.  What it takes from the main analysis are
physics or method constants that do not depend on the dataset (user rule 2026-09-25: e.g. shapes that barely change
between datasets need no extensive search):
  * the true Z1 line shape of the Z1 kinematic refit and the FSR-photon pT resolution (signal MC generator level,
    production_v3/h4l_reco/refit/v3/refit_inputs.json);
  * the YR4 effective cross sections sigma_eff(mode, m_H) of the main analysis's likelihood (production_v3/signal_model/
    yr4/v1/yr4.json, tabulated 120-130 GeV), as ratios to 125 GeV with the log-quadratic extrapolation outside the table
    (h4l_likelihood.SigmaEff);
  * the method systematics of the lepton calibration from its closures (analysis_v3/inference/config/
    systematics_ul16_v3.json: scale closure bias and iteration, smear variance closure deviation, per flavour);
  * fallbacks used only when a dataset cannot provide the in-situ quantity (logged): the lambda of the main analysis's
    final selection (MC role), the Z+X shapes and the relative MC-closure systematic (production_v3/backgrounds/v5/zx_v2).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
OUT = REPO / "my_analysis/h4l_eval/data/constants.json"


def rounded(values, digits=7):
    return [float(f"{v:.{digits}g}") for v in values]


def main() -> int:
    refit = json.loads((PRODUCTION / "h4l_reco/refit/v3/refit_inputs.json").read_text(encoding="utf-8"))
    z1 = refit["z1_lineshape"]
    edges = np.array(z1["edges"])
    density = np.array(z1["density"])
    log_density = np.log(np.clip(density, 1e-300, None))
    plan = json.loads((PRODUCTION / "h4l_select/v5/plan.json").read_text(encoding="utf-8"))
    lam_payload = plan["tasks"][0]["lambda_payload"]["flavours"]
    zx = json.loads((PRODUCTION / "backgrounds/v5/zx_v2/zx.json").read_text(encoding="utf-8"))
    syst = json.loads((REPO / "analysis_v3/inference/config/systematics_ul16_v3.json").read_text(encoding="utf-8"))
    lambda_fallback = {"13": [float(r["mc"]["lambda"]) for r in lam_payload["muon"]["regions"]],
                       "11": [float(r["mc"]["lambda"]) for r in lam_payload["electron"]["regions"]]}
    yr4 = json.loads((PRODUCTION / "signal_model/yr4/v1/yr4.json").read_text(encoding="utf-8"))
    sigma_eff = {}
    for mode, key in (("ggH", "GluGluToHToZZ"), ("VBF", "VBF_HToZZ"), ("VH", "VHToZZ")):
        rows = yr4["sigma_eff"][key]
        masses = np.array([r["mH"] for r in rows])
        values = np.array([r["sigma_eff_pb"] for r in rows])
        ref = float(np.interp(125.0, masses, values))
        sigma_eff[mode] = {"mass": masses.tolist(), "ratio": (values / ref).tolist(),
                           "log_quadratic": np.polyfit(masses - 125.0, np.log(values / ref), 2).tolist()}
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
        "sigma_eff": sigma_eff,
        "calibration_systematics": calib_syst,
        "lambda_fallback": lambda_fallback,
        "zx_fallback": {"shapes": {fs: {k: zx["shapes"][fs][k] for k in ("form", "range", "mpv", "width", "exp_fraction", "exp_slope")}
                                   for fs in ("4mu", "4e", "2e2mu")},
                        "relative_systematic": zx["relative_systematic"],
                        "source": "production_v3/backgrounds/v5/zx_v2/zx.json (main analysis, 20 fb^-1 pseudo-data and MC closure)"},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(constants, indent=1) + "\n", encoding="utf-8")
    print(f"[constants] {OUT} ({OUT.stat().st_size / 1024:.1f} kB)")
    for mode, t in sigma_eff.items():
        print(f"  sigma_eff {mode}: table {t['mass'][0]}-{t['mass'][-1]} GeV, ratio(120) {np.interp(120.0, t['mass'], t['ratio']):.4f}, "
              f"ratio(130) {np.interp(130.0, t['mass'], t['ratio']):.4f}, log-quadratic {t['log_quadratic']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
