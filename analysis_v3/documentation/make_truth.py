"""The truth of the pseudo-data for the comparison after unblinding (user, 2026-09-25/26: the injection record
slopbench_code_fork/injection/records/20260923T183922Z_917f9326efa44ea9.json; its profiles interpreted with the
generator's catalog and configuration, benchmark/configs/h4l_pfnano_seed_plan_v3.json, and pfnano_response.cpp).

    pixi run py -- analysis_v3/documentation/make_truth.py --label truth_v1

Truth: m_H = 125 GeV, mu = 1 (every mode 1), 20 fb^-1, response profile "medium", efficiency profile "effA".
Lepton scale factor k(pT, eta) = 1 - 0.018 - 0.004 S(pT / 100 GeV) - 0.008 S(|eta| / 2.5), S(x) = x^2 / (1 + x^2)
(pfnano_response.cpp scaleFactor, saturating_downward); the truth scale_shift per flavour is its average over the same
data legs as the RESULT.json average (the Z -> ll calibration sample, corrected pT as the true pT); lepton smear: relative
pT sigma 0.01 (plus eta and phi smears of 0.01); efficiency effA thins Muon_loose/medium/tight/softId, the electron
MVA WP90 / WP80 and cut-based flags and the single-lepton trigger decisions: the flags of this analysis's selection (PF
muon, mvaFall17V2noIso WPL electron, isolation, SIP) are not thinned, so the truth sel_eff of the full lepton selection is
1 for both flavours.  Writes production_v3/results/v5/<label>/truth.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
import data_weighted_average as dwa  # noqa: E402

RECORD = Path("/eos/home-y/yiyangz/codex/slopbench_code_fork/injection/records/20260923T183922Z_917f9326efa44ea9.json")


def saturating_square(x):
    x = np.abs(x)
    return x * x / (1 + x * x)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    out = REPO / "production_v3/results/v5" / args.label
    if (out / "truth.json").exists():
        raise SystemExit("exists")
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    cfg = record["configuration"]
    run = json.loads((REPO / "production_v3/calibration/v4/nominal/run.json").read_text(encoding="utf-8"))
    payload = json.loads((REPO / "production_v3/calibration/v4/nominal/payload.json").read_text(encoding="utf-8"))
    config = json.loads(Path(run["config"]).read_text(encoding="utf-8"))
    scan = json.loads(Path(run["extract_scan"]).read_text(encoding="utf-8"))
    files = sorted(o["path"] for o in scan["outputs"] if o["kind"] == "data")
    with ProcessPoolExecutor(max_workers=9) as pool:
        parts = list(pool.map(dwa.read_legs, files))
    legs = {k: np.concatenate([p[k] for p in parts]) for k in dwa.BRANCHES}
    compare = {"mu": float(cfg["mu"]), "mH": float(cfg["higgs_mass_gev"])}
    detail = {}
    for flavour, (pdg, _) in dwa.FLAVOURS.items():
        model = dwa.Factorized(payload["flavours"][flavour])
        sel_cfg = config["selection"][flavour]
        pair = legs["flavour"] == pdg
        passing, cols = None, {}
        for i in (1, 2):
            pt = legs[f"l{i}_pt"][pair].astype(float)
            abs_eta_model = np.abs(legs[f"l{i}_eta_sc" if model.eta_sc else f"l{i}_eta"][pair].astype(float))
            eta = legs[f"l{i}_eta"][pair].astype(float)
            e, r = model.bins(pt, abs_eta_model)
            inside = (e >= 0) & (r >= 0) & (pt >= model.pt[0])
            u = np.where(inside, model.evaluate(pt, np.where(inside, abs_eta_model, 0.0), model.a, model.b), 0.0)
            ident = dwa.muon_an_tight(legs[f"l{i}_flags"][pair]) if pdg == 13 else dwa.electron_wp90(legs[f"l{i}_flags"][pair])
            ok = (inside & ident & (legs[f"l{i}_iso_fsr"][pair] < sel_cfg["max_iso_fsr"]) & (legs[f"l{i}_sip"][pair] < sel_cfg["max_sip"])
                  & (pt * np.exp(-u) >= sel_cfg["min_pt"]) & (np.abs(eta) < sel_cfg["max_abs_eta"]))
            cols[i] = (pt * np.exp(-u), eta, u)
            passing = ok if passing is None else passing & ok
        mass = legs["mass"][pair].astype(float) * np.exp(-0.5 * (cols[1][2] + cols[2][2]))
        passing &= (mass > 60.0) & (mass < 120.0)
        pt_true = np.concatenate([cols[i][0][passing] for i in (1, 2)])
        eta = np.concatenate([cols[i][1][passing] for i in (1, 2)])
        s_true = -(0.018 + 0.004 * saturating_square(pt_true / 100.0) + 0.008 * saturating_square(eta / 2.5))
        compare[f"scale_shift {flavour}"] = float(np.mean(s_true))
        compare[f"smear {flavour}"] = 0.01
        compare[f"sel_eff {flavour}"] = 1.0
        detail[flavour] = {"n_legs": int(len(pt_true)), "scale_shift_at_45GeV_eta1p2": float(-(0.018 + 0.004 * saturating_square(0.45)
                                                                                               + 0.008 * saturating_square(1.2 / 2.5)))}
    truth = {"schema": "h4l_v3_truth/1", "record": str(RECORD), "configuration": cfg, "signal_injected": record["signal"]["total_injected"],
             "compare": compare, "detail": detail,
             "source": ("Truth from the injection record 20260923T183922Z_917f9326efa44ea9 (m_H = 125 GeV, mu = 1 for every mode, "
                        "20 fb^-1, response profile medium, efficiency profile effA) with the profile definitions of the generator "
                        "(h4l_pfnano_seed_plan_v3.json, pfnano_response.cpp)."),
             "notes": ["Lepton scale: k(pT, eta) = 1 - 0.018 - 0.004 S(pT/100 GeV) - 0.008 S(|eta|/2.5), S(x) = x^2/(1+x^2); the truth "
                       "scale_shift per flavour is its average over the data legs of the Z -> ll calibration sample (as RESULT.json); at "
                       "pT = 45 GeV, |eta| = 1.2 it is -0.0202.",
                       "Lepton smear: relative pT sigma 0.01 plus Gaussian eta and phi smears of 0.01; the pT-only smear model of the "
                       "analysis absorbs the angular smears, so the measured smear is expected somewhat above 0.01.",
                       "Efficiency effA thins the muon loose/medium/tight/soft IDs, the electron MVA WP90/WP80 and cut-based flags and "
                       "the single-lepton trigger decisions, none of the flags of this selection (PF muon, mvaFall17V2noIso WPL "
                       "electron, isolation, SIP): the truth data/MC efficiency ratio of the lepton selection is 1; the trigger "
                       "thinning acts on the tag selection and the event trigger (taken from the MC)."]}
    out.mkdir(parents=True, exist_ok=True)
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(compare, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
