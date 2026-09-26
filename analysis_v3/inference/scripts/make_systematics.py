"""The experimental nuisance magnitudes of the likelihood from the calibration products (stage 7).

    pixi run py -- analysis_v3/inference/scripts/make_systematics.py --calibration v4/nominal \
        --closures v4/closure_uniform_halves v4/closure_sloped_halves --signal v5/sm_v1 --select v5 \
        --base analysis_v3/inference/config/systematics_ul16_v1.json --out analysis_v3/inference/config/systematics_ul16_v2.json

Lepton momentum scale (relative, per flavour), the nuisance scale_<f> moving every lepton of the flavour:
    sqrt(stat^2 + iteration^2 + bias^2),
stat and iteration the uncertainties of the calibration at the report point (45 GeV, |eta| 1.2), bias the
largest rms over the closure profiles (|eta| and pT, points below 100 GeV where the signal leptons are) of
measured - injected scale among the independent-halves closures (the method's accuracy, template noise
included).
Lepton resolution (relative m4l width, per flavour), the nuisance res_<f>: the uncertainty of the smear
variance (the nominal statistical one and the largest closure deviation at the report point, in quadrature)
propagated to the width, 1/2 dv / sigma_l^2 with sigma_l = 2 sigma(m4l)/m4l of the four-lepton final state of
that flavour (a lepton carries about a quarter of the m4l^2 derivative), in quadrature with the residual
data/MC width difference of the dilepton peak after the calibration (rms over the charge x eta and charge x phi
slices of |IQR ratio - 1|, amplified by IQR^2 / IQR_res^2 to remove the Z natural width; IQR_res = 1.35 sigma_Z
with sigma_Z the resolution part of the dilepton peak, sigma_l m_Z / sqrt 2).
Flavour fractions per final state: scale, the mean m4l response to the flavour's scale in the signal MC
(dm4l_scale_<f> / m4l); resolution, the variance shares of 2e2mu from the 4mu and 4e widths (sigma_4mu^2 and
sigma_4e^2 halves), 1 for the pure final states.  Everything else is taken from the base configuration.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_select_io as io  # noqa: E402

Z_MASS = 91.1876
SIGNAL = ["GluGluToHToZZ_M125", "VBF_HToZZ_M125", "VHToZZ_M125"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--calibration", required=True)
    parser.add_argument("--closures", nargs="+", required=True)
    parser.add_argument("--signal", required=True, help="<select>/<signal-model label>")
    parser.add_argument("--select", required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out if args.out.is_absolute() else REPO / args.out
    if out.exists():
        raise SystemExit(f"{out} exists")
    base = json.loads((args.base if args.base.is_absolute() else REPO / args.base).read_text(encoding="utf-8"))
    cal = PRODUCTION / "calibration" / args.calibration
    report = json.loads((cal / "plots" / "report.json").read_text(encoding="utf-8"))
    diag = json.loads((cal / "diag_final" / "diagnostics.json").read_text(encoding="utf-8"))
    closures = {c: json.loads((PRODUCTION / "calibration" / c / "closure" / "closure.json").read_text(encoding="utf-8"))
                for c in args.closures}
    sm = json.loads((PRODUCTION / "signal_model" / args.signal / "signal_model.json").read_text(encoding="utf-8"))
    rel_width = {fs: sm["final_states"][fs]["m4l_refit"]["dcb"]["values"]["width"] / sm["final_states"][fs]["m4l_refit"]["dcb"]["values"]["mean"]
                 for fs in ("4mu", "4e", "2e2mu")}
    details = {}
    scale, resolution = {}, {}
    for key, flavour, fs, tag in (("mu", "muon", "4mu", "mm"), ("e", "electron", "4e", "ee")):
        f = report["flavours"][flavour]
        stat, iteration = f["scale_shift"]["stat"], f["scale_shift"]["iteration"]
        biases = {}
        dv_closure = 0.0
        for name, c in closures.items():
            # Profile points above 100 GeV (a few per cent of the signal leptons) are left out of the rms.
            prof = [q for q in c["flavours"][flavour].get("profile", []) if q["pt"] < 100.0]
            if prof:
                biases[name] = float(np.sqrt(np.mean([(p["s"] - p["s_true"]) ** 2 for p in prof])))
            sv = c["flavours"][flavour]["smear_variance"]
            dv_closure = max(dv_closure, abs(sv["measured"] - sv["injected"]))
        bias = max(biases.values()) if biases else 0.0
        scale[key] = math.sqrt(stat ** 2 + iteration ** 2 + bias ** 2)
        sigma_l = 2.0 * rel_width[fs]
        dv = math.hypot(f["smear_variance"]["stat"], dv_closure)
        from_smear = 0.5 * dv / sigma_l ** 2
        ratios = np.concatenate([np.abs(np.array(diag["slices"][tag][s]["iqr_ratio"]) - 1.0) for s in diag["slices"][tag]])
        sigma_z = sigma_l * Z_MASS / math.sqrt(2.0)
        iqr_res = 1.35 * sigma_z
        iqr_total = math.hypot(2.4952, iqr_res)
        from_width = float(np.sqrt(np.mean(ratios ** 2))) * iqr_total ** 2 / iqr_res ** 2
        resolution[key] = math.hypot(from_smear, from_width)
        details[key] = {"scale": {"stat": stat, "iteration": iteration, "closure_bias_rms": biases, "total": scale[key]},
                        "resolution": {"sigma_lepton_rel": sigma_l, "smear_variance_stat": f["smear_variance"]["stat"],
                                       "smear_variance_closure_max_dev": dv_closure, "from_smear": from_smear,
                                       "dilepton_iqr_rms_dev": float(np.sqrt(np.mean(ratios ** 2))), "iqr_amplification": iqr_total ** 2 / iqr_res ** 2,
                                       "from_width": from_width, "total": resolution[key]}}
    # Flavour fractions from the signal MC.
    scan = io.load_scan(args.select)
    rows = io.read(scan, SIGNAL, "SR", ["final_state", "m4l", "m4l_refit", "dm4l_scale_mu", "dm4l_scale_e"],
                   cut=lambda a: (a["m4l_refit"] > 105) & (a["m4l_refit"] < 140))
    fractions = {}
    for code, fs in ((0, "4mu"), (1, "4e"), (2, "2e2mu")):
        sel = rows["final_state"] == code
        w = rows["w"][sel]
        f_mu = float(np.sum(w * rows["dm4l_scale_mu"][sel] / rows["m4l"][sel]) / np.sum(w))
        f_e = float(np.sum(w * rows["dm4l_scale_e"][sel] / rows["m4l"][sel]) / np.sum(w))
        fractions[fs] = [f_mu, f_e]
    var_mu, var_e = rel_width["4mu"] ** 2, rel_width["4e"] ** 2
    res_fractions = {"4mu": [1.0, 0.0], "4e": [0.0, 1.0], "2e2mu": [var_mu / (var_mu + var_e), var_e / (var_mu + var_e)]}
    config = dict(base)
    config["version"] = "syst_v2"
    config["comment"] = (base.get("comment", "") + "  v2: the lepton scale and resolution magnitudes and the flavour fractions derived by "
                         "make_systematics.py from the calibration, its closures and the signal MC (details block).")
    config["scale_uncertainty"] = scale
    config["resolution_uncertainty"] = resolution
    config["flavour_fraction"] = fractions
    config["resolution_flavour_fraction"] = res_fractions
    config["details"] = {"calibration": str(cal), "closures": args.closures, "signal_model": args.signal, "select": args.select,
                         "per_flavour": details}
    out.write_text(json.dumps(config, indent=1) + "\n", encoding="utf-8")
    print(f"[syst] scale {scale}, resolution {resolution}\n[syst] scale fractions {fractions}\n[syst] resolution fractions "
          f"{res_fractions}\n[syst] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
