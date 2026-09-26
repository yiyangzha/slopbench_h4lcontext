"""MC closure of the fit on real simulated events (final review 2026-09-25): the model fitted to the MC itself.

    pixi run py -- analysis_v3/inference/scripts/validate_mc_closure.py --select v5 --model model_cat_r1 --yr4 v1 \
        --label closure_mc_v1 [--dimensions 1D 2Dmass 3D]

The pseudo-data are the simulated events themselves, weighted to the expected yields: the ggH, VBF and VH MC at
m_H = 125 GeV and the qqZZ and ggZZ MC (weights x scale factors, the model's channels and window), plus the Z+X
component as weighted points of its model density (as in the Asimov grid).  Unlike the Asimov dataset and the
toys, which are drawn from the model, these events carry the true per-event mass resolution and its correlation
with m4l (e.g. the lower m4l of poorly measured electrons), so the fit tests the signal model itself: the fitted mu
and m_H against the truth (1 and 125 GeV) per dimension, with the statistical precision of the MC.  Writes
production_v3/inference/<select>/<model>/<label>/closure.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/inference"))
sys.path.insert(0, str(REPO / "analysis_v3/inference/scripts"))
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import fit_model as fm  # noqa: E402
import h4l_categories as cats  # noqa: E402
import h4l_likelihood as lk  # noqa: E402
import h4l_select_io as io  # noqa: E402
import h4l_sf  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

SAMPLES = {"GluGluToHToZZ_M125": "signal", "VBF_HToZZ_M125": "signal", "VHToZZ_M125": "signal", "ZZTo4L": "qqZZ",
           "GGZZ4Mu": "ggZZ", "GGZZ4E": "ggZZ", "GGZZ2E2Mu": "ggZZ"}
FS_CODES = {"4mu": 0, "4e": 1, "2e2mu": 2}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--dimensions", nargs="+", default=["1D", "2Dmass", "3D"])
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.select / args.model
    out = model_dir / args.label
    if out.exists():
        raise SystemExit(f"{out} exists")
    model, _, signal_model, yr4 = fm.load(model_dir, args.yr4)
    syst = json.loads(args.systematics.read_text(encoding="utf-8"))
    categories = signal_model["categories"]
    sf = h4l_sf.LeptonSF(model["sf_payload"])
    scan = io.load_scan(args.select)
    lo, hi = model["window"]
    branches = ["final_state", "m4l", "m4l_err", "m4l_refit", "m4l_refit_err", "d_bkg_kin", "l_pdg", "l_pt", "l_eta", "l_eta_sc",
                "pt4l", "jet_pt", "jet_eta"] + cats.CATEGORY_BRANCHES
    rows = {}
    for sample, kind in SAMPLES.items():
        r = io.read(scan, [sample], "SR", branches, cut=lambda a: (a["m4l_refit"] > lo - 5) & (a["m4l_refit"] < hi + 5))
        r["category"] = cats.category_index(r, categories)
        r["sf"] = sf.event(r)["sf"]
        rows[sample] = r
    report = {"schema": "h4l_v3_mc_closure/1", "model": str(model_dir), "truth": {"mu": 1.0, "mH": 125.0}, "fits": {}}
    for dim in args.dimensions:
        for refit in (True, False):
            lik = lk.Model(model, signal_model, yr4, dim, refit, "b", syst, "inclusive")
            var, err = lik.variable, lik.error
            events = {}
            p0 = lik.default_parameters(1.0, 125.0)
            for ch_name, ch in model["channels"].items():
                fs, cat = ch["final_state"], categories["order"].index(ch["category"])
                parts = {"m": [], "dm": [], "d": [], "w": []}
                for sample, r in rows.items():
                    sel = (r["final_state"] == FS_CODES[fs]) & (r["category"] == cat) & (r[var] > lo) & (r[var] < hi)
                    parts["m"].append(r[var][sel])
                    parts["dm"].append(r[err][sel])
                    parts["d"].append(r["d_bkg_kin"][sel])
                    parts["w"].append((r["w"] * r["sf"])[sel])
                # Z+X: weighted points of its model density on the Asimov grid (its yield in the channel).
                grid = fm.asimov_events(lik, p0)[ch_name]
                fse = model["final_states"][fs]
                zx_dens = shapes.zx_pdf(grid["m"], fse["zx"]["shape"], lik.lo, lik.hi)
                if lik.use_d:
                    zx_dens = zx_dens * lik.d_density("zx", fs, grid["m"], grid["d"])
                if lik.use_e:
                    zx_dens = zx_dens * lik.e_density("zx", fs, grid["dm"] / grid["m"])
                step = 0.25
                n_e = fm.E_SUB if lik.use_e else 1.0
                parts["m"].append(grid["m"])
                parts["dm"].append(grid["dm"])
                parts["d"].append(grid["d"])
                parts["w"].append(lik.yields(p0, ch_name)["zx"] * zx_dens * step / n_e)
                events[ch_name] = {k: np.concatenate(v) for k, v in parts.items()}
            fitter = fm.Fitter(lik, events)
            best = fitter.fit()
            key = f"{dim}_{'refit' if refit else 'norefit'}"
            report["fits"][key] = {"mu": float(best.values["mu"]), "mH": float(best.values["mH"]), "mu_error": float(best.errors["mu"]),
                                   "mH_error": float(best.errors["mH"]), "valid": bool(best.valid),
                                   "mH_bias": float(best.values["mH"]) - 125.0, "mu_bias": float(best.values["mu"]) - 1.0}
            print(f"[closure] {key}: mu {best.values['mu']:.4f}, mH {best.values['mH']:.4f} (+- {best.errors['mH']:.3f} data-like), "
                  f"valid {best.valid}", flush=True)
    out.mkdir(parents=True)
    (out / "closure.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[closure] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
