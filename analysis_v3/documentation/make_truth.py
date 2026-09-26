"""The truth of the pseudo-data for the comparison after unblinding (user, 2026-09-25/26: the injection record
slopbench_code_fork/injection/records/20260923T183922Z_917f9326efa44ea9.json; its profiles interpreted with the
generator's catalog and configuration, benchmark/configs/h4l_pfnano_seed_plan_v3.json, and pfnano_response.cpp).

    pixi run py -- analysis_v3/documentation/make_truth.py --label truth_v1

Truth: m_H = 125 GeV, mu = 1 (every mode 1), 20 fb^-1, response profile "medium", efficiency profile "effA".
Lepton scale factor k(pT, eta) = 1 - 0.018 - 0.004 S(pT / 100 GeV) - 0.008 S(|eta| / 2.5), S(x) = x^2 / (1 + x^2)
(pfnano_response.cpp scaleFactor, saturating_downward); the truth scale_shift per flavour is its average over the same
data legs as the RESULT.json average (the Z -> ll calibration sample, corrected pT as the true pT); lepton smear: relative
pT sigma 0.01 (plus eta and phi smears of 0.01).  Efficiency effA (pfnano_response.cpp applyMuonEfficiency /
applyElectronEfficiency), one uniform draw per lepton and family, efficiency ratio r_B(pT, eta) = (1 - offset)
(1 - pt_coefficient / (1 + (pT/20)^2)) (1 - eta_coefficient S(|eta|/1.5)):
  * muons: a killed Muon_looseId (offset 0.0025, pT 0.01, eta 0.005) also clears Muon_isPFcand (the generator keeps
    looseId = PF and (global or tracker) consistent) and may remove the muon from the collection, so the PF muon of this
    selection survives with r_loose(pT, eta): truth sel_eff = r_loose (0.9939 at 45 GeV, |eta| 1.2);
  * electrons: mvaFall17V2noIso WPL is not configured and is cleared only where its constant low-pT cut lies above the
    failed WP90 / WP80 curve (pT < 10 GeV MVA categories: in practice endcap electrons below ~8.6 GeV), so the truth
    sel_eff is 1 except there (r_WP90 or r_WP80);
  * the trigger families kill trigger objects (the event trigger efficiency, taken from the MC in this analysis).
The per-flavour truth values are averages over the same data legs as RESULT.json.  Writes
production_v3/results/v5/<label>/truth.json.
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


def efficiency_ratio(pt, eta, offset, pt_coefficient, eta_coefficient):
    """pfnano_response.cpp efficiencyRatio with the effA common shapes (lowpt_turnon 20 GeV, saturating_square 1.5)."""
    return (1 - offset) * (1 - pt_coefficient / (1 + (pt / 20.0) ** 2)) * (1 - eta_coefficient * saturating_square(eta / 1.5))


# Fall17 V2 noIso thresholds of the generator (C - A exp(-pT / tau), raw MVA), categories (pT < 10 | >= 10) x (EB1, EB2, EE).
NOISO_WPL = [(0.89441115863, 0.0, 1.0), (0.79196646463, 0.0, 1.0), (1.4710485717, 0.0, 1.0),
             (-0.29396295866, 0.0, 1.0), (-0.25042475858, 0.0, 1.0), (-0.13098517903, 0.0, 1.0)]
NOISO_WP90 = [(2.7707238734, 8.1630486018, 3.8150091215), (1.8560231781, 11.856893682, 2.1869765494),
              (1.7348930781, 17.013880078, 2.0163211971), (5.9175992258, 9.3196623269, 13.480729454),
              (5.0159883725, 8.7941819377, 13.12804515), (4.1692134321, 9.0072091321, 13.201722462)]
NOISO_WP80 = [(3.2644962047, 8.8466978357, 3.3265714922), (2.835578385, 11.097801657, 2.1515048765),
              (2.9199494518, 24.024807824, 1.6987547752), (7.1336238874, 8.2253122239, 16.56052688),
              (6.1863827578, 7.4976456532, 15.269463428), (5.4317586574, 7.5689969229, 15.429007595)]


def electron_wpl_ratio(pt, eta, eta_sc):
    """Survival ratio of the WPL flag under effA: WPL is cleared with a killed WP80 where its cut is above the WP80 curve,
    else with a killed WP90 where its cut is above the WP90 curve, else never."""
    abs_sc = np.abs(eta_sc)
    cat = np.where(pt < 10.0, 0, 3) + np.where(abs_sc < 0.8, 0, np.where(abs_sc < 1.479, 1, 2))
    out = np.ones(len(pt))
    r90 = efficiency_ratio(pt, eta, 0.005, 0.02, 0.01)
    r80 = efficiency_ratio(pt, eta, 0.0075, 0.03, 0.015)
    for k in range(6):
        sel = cat == k
        wpl = NOISO_WPL[k][0]
        t90 = NOISO_WP90[k][0] - NOISO_WP90[k][1] * np.exp(-pt[sel] / NOISO_WP90[k][2])
        t80 = NOISO_WP80[k][0] - NOISO_WP80[k][1] * np.exp(-pt[sel] / NOISO_WP80[k][2])
        out[sel] = np.where(wpl >= t80, r80[sel], np.where(wpl >= t90, r90[sel], 1.0))
    return out


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
            cols[i] = (pt * np.exp(-u), eta, u, legs[f"l{i}_eta_sc"][pair].astype(float))
            passing = ok if passing is None else passing & ok
        mass = legs["mass"][pair].astype(float) * np.exp(-0.5 * (cols[1][2] + cols[2][2]))
        passing &= (mass > 60.0) & (mass < 120.0)
        pt_true = np.concatenate([cols[i][0][passing] for i in (1, 2)])
        eta = np.concatenate([cols[i][1][passing] for i in (1, 2)])
        eta_sc = np.concatenate([cols[i][3][passing] for i in (1, 2)])
        s_true = -(0.018 + 0.004 * saturating_square(pt_true / 100.0) + 0.008 * saturating_square(eta / 2.5))
        compare[f"scale_shift {flavour}"] = float(np.mean(s_true))
        compare[f"smear {flavour}"] = 0.01
        if pdg == 13:
            eff = efficiency_ratio(pt_true, eta, 0.0025, 0.01, 0.005)
            eff45 = float(efficiency_ratio(np.array([45.0]), np.array([1.2]), 0.0025, 0.01, 0.005)[0])
        else:
            eff = electron_wpl_ratio(pt_true, eta, eta_sc)
            eff45 = float(electron_wpl_ratio(np.array([45.0]), np.array([1.2]), np.array([1.2]))[0])
        compare[f"sel_eff {flavour}"] = float(np.mean(eff))
        detail[flavour] = {"n_legs": int(len(pt_true)), "scale_shift_at_45GeV_eta1p2": float(-(0.018 + 0.004 * saturating_square(0.45)
                                                                                               + 0.008 * saturating_square(1.2 / 2.5))),
                           "sel_eff_at_45GeV_eta1p2": eff45}
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
                       "Efficiency effA: a killed Muon_looseId (ratio (1-0.0025)(1-0.01/(1+(pT/20)^2))(1-0.005 S(|eta|/1.5))) also "
                       "clears Muon_isPFcand and may remove the muon, so the PF-muon selection has the truth ratio r_loose (0.9939 at "
                       "45 GeV, |eta| 1.2); the electron mvaFall17V2noIso WPL flag is cleared only where its low-pT cut lies above a "
                       "killed WP90/WP80 curve (endcap electrons below about 8.6 GeV), so the electron truth ratio is 1 for the "
                       "Z legs; the trigger families remove trigger objects (event trigger efficiency, taken from the MC)."]}
    out.mkdir(parents=True, exist_ok=True)
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(compare, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
