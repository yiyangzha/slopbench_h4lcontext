"""The result files of stage 8: RESULT.json, results_full.json, summary.md, the judge-readable tables and the
paper-comparison figure.

    pixi run py -- analysis_v3/inference/scripts/make_results.py --select v5 --tag r1 --model-cat model_cat_v1 \
        --model-incl model_incl_v1 --model-pt4l model_pt4l_v1 --model-njets model_njets_v1 --model-ptj1 model_ptj1_v1 \
        --fiducial production_v3/signal_model/v5/gen_v1/fiducial_xsec.json --yr4 v1 --calibration v4/nominal \
        --sf production_v3/tnp/v2/v7/<label>/sf.json --nm1 production_v3/optimization/v5/nm1_v1/nm1.json \
        --average production_v3/results/v5/average_v1/average.json [--pyhf production_v3/results/v5/r1/MODEL.json]

Inputs: the outputs of run_results.py (the same arguments name them), the data-weighted per-flavour averages of
the calibration and the tag-and-probe SF (data_weighted_average.py; RESULT.json), the calibration report (scale
shift and smear at pT = 45 GeV, |eta| = 1.2) and the tag-and-probe SF payload (the full single-lepton SF at the
same point) for the AN, the systematics configuration of the fits, the N-1 table and the paper values
(analysis_v3/inference/config/paper_jhep11_2017_047.json).  Conventions (AGENTS.md): RESULT.json follows
EVAL_CONTRACT v0.4 with result set (b) and the per-flavour calibration block (user 2026-09-25: {"scale_shift",
"smear", "sel_eff"} x {"muon", "electron"}, no flavour-combined value); the declared cuts and the systematic groups
(set (a)) go to results_full.json; uncertainties are symmetrized (half the MINOS interval), syst from the
fixed-nuisance impacts.
Writes production_v3/results/<select>/<tag>/.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/inference/scripts"))
import run_results as rr  # noqa: E402

PAPER = REPO / "analysis_v3/inference/config/paper_jhep11_2017_047.json"
DISCOVERY_Z = 5.0
FINAL_STATES = ("4mu", "4e", "2e2mu")


def sym(interval) -> float:
    return 0.5 * (abs(interval[0]) + abs(interval[1]))


def poi(fit: dict, name: str) -> dict:
    r = fit["fit"]["pois"][name]
    out = {"value": r["value"], "total": r["total"], "stat": r["stat"], "total_sym": sym(r["total"]), "stat_sym": sym(r["stat"])}
    if "syst" in r:
        out["syst"] = r["syst"]
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("--select", "--tag", "--model-cat", "--model-incl", "--model-pt4l", "--model-njets", "--model-ptj1", "--yr4",
                 "--calibration"):
        parser.add_argument(name, required=True)
    parser.add_argument("--fiducial", type=Path, required=True)
    parser.add_argument("--sf", type=Path, required=True)
    parser.add_argument("--nm1", type=Path, required=True)
    parser.add_argument("--pyhf", type=Path, default=None)
    parser.add_argument("--average", type=Path, required=True,
                        help="average.json of data_weighted_average.py: the per-flavour calibration of RESULT.json")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--scan-workers", type=int, default=1)
    parser.add_argument("--gof-toys", type=int, default=500)
    parser.add_argument("--only", default=None)
    parser.add_argument("--out-tag", default=None, help="output directory tag (default --tag, the tag of the run_results fits)")
    parser.add_argument("--toys", nargs="*", default=[], help="run_toys labels of the categories model (validation)")
    parser.add_argument("--paired", default=None, help="run_toys --paired label")
    parser.add_argument("--mc-closure", default=None, help="validate_mc_closure.py label")
    parser.add_argument("--tnp-closure", nargs="*", default=[], help="production_v3 paths of the T&P efficiency closures")
    args = parser.parse_args()
    for key in ("fiducial", "sf", "nm1", "pyhf", "average"):
        value = getattr(args, key)
        if value is not None and not value.is_absolute():
            setattr(args, key, REPO / value)
    out_dir = PRODUCTION / "results" / args.select / (args.out_tag or args.tag)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("RESULT.json", "results_full.json", "summary.md"):
        if (out_dir / name).exists():
            raise SystemExit(f"{out_dir / name} exists; results are never overwritten (use a new --tag)")
    outputs = {j["name"]: j["output"] for j in rr.jobs(args)}
    missing = [n for n, p in outputs.items() if not p.exists()]
    if missing:
        raise SystemExit(f"missing run_results outputs: {missing}")
    load = lambda name: json.loads(outputs[name].read_text(encoding="utf-8"))  # noqa: E731
    paper = json.loads(PAPER.read_text(encoding="utf-8"))
    model = json.loads((PRODUCTION / "inference" / args.select / args.model_cat / "model.json").read_text(encoding="utf-8"))
    fiducial = json.loads(args.fiducial.read_text(encoding="utf-8"))
    b_data, b_asimov, b_post = load("b3d_data"), load("b3d_asimov"), load("b3d_asimov_postfit")
    a_data, a_asimov = load("a3d_data"), load("a3d_asimov")
    syst = b_data["systematics"]

    # ---------------------------------------------------------------- calibration and efficiency
    report = json.loads((PRODUCTION / "calibration" / args.calibration / "plots" / "report.json").read_text(encoding="utf-8"))
    calib = {}
    for flavour, key in (("muon", "mu"), ("electron", "e")):
        f = report["flavours"][flavour]
        details = syst["details"]["per_flavour"][key]["resolution"]
        smear = f["smear"]["value"]
        dv = math.hypot(details["smear_variance_stat"], details["smear_variance_closure_max_dev"])
        calib[flavour] = {
            "scale_shift": {"value": f["scale_shift"]["value"], "unc": syst["scale_uncertainty"][key],
                            "unc_components": syst["details"]["per_flavour"][key]["scale"]},
            "smear": {"value": smear, "unc": dv / (2.0 * smear) if smear > 0 else None,
                      "unc_components": {"variance_stat": details["smear_variance_stat"],
                                         "variance_closure_max_dev": details["smear_variance_closure_max_dev"]}}}
    sf = json.loads(args.sf.read_text(encoding="utf-8"))
    effscale = {}
    for flavour in ("muon", "electron"):
        rp = sf["flavours"][flavour]["report_point"]["sf_full"]
        effscale[flavour] = {"value": rp.get("value"), "unc": rp.get("error"),
                             "definition": "tag-and-probe SF of the full single-lepton selection (chain id | loose, sip | id, "
                                           "iso | id and sip) at pT = 45 GeV, |eta| = 1.2 (bilinear between bin centres)"}
    # Cross-check (user decision 2026-09-25: reconstruction efficiency from the MC, the yield ratio not propagated): the
    # data/MC ratio of the Z -> ll yields of the calibration pairs (both legs passing the full selection, no SFs on the
    # MC) against the square of the inclusive tag-and-probe SF(full).
    yield_ratio = {}
    for flavour in ("muon", "electron"):
        ratio = report["flavours"][flavour]["zpeak_totals"]["before"]["mc_normalization"]
        inclusive = next(r for r in sf["flavours"][flavour]["product"] if r["bin"] == "all")
        yield_ratio[flavour] = {"z_yield_ratio_data_over_mc": ratio, "tnp_sf_full_inclusive": inclusive["value"],
                                "tnp_sf_full_squared": inclusive["value"] ** 2,
                                "remaining_per_lepton": math.sqrt(ratio / inclusive["value"] ** 2)}
    working_points = {"muon": "AN loose muon (global or tracker, |dxy| < 0.5 cm, |dz| < 1 cm) + PF ID (tracker high-pT ID above "
                              "200 GeV) + SIP3D < 4 + FSR-subtracted relative isolation (dR 0.3) < 0.35",
                      "electron": "AN loose electron (|dxy| < 0.5 cm, |dz| < 1 cm) + mvaFall17V2noIso WPL + SIP3D < 4 + "
                                  "FSR-subtracted relative isolation (dR 0.3) < 0.35"}

    # ---------------------------------------------------------------- RESULT.json
    mu_b, mh_b = poi(b_data, "mu"), poi(b_data, "mH")
    z_obs = b_data["fit"]["significance"]["Z"]
    gof = load("gof")

    def contract_poi(name, p):
        return {"name": name, "value": p["value"], "stat": p["stat_sym"], "syst": p["syst"], "total": p["total_sym"]}

    # Per-flavour calibration (EVAL_CONTRACT v0.4 with the user's per-flavour format, 2026-09-25): the registered names
    # scale_shift, smear and sel_eff, each {"muon": {value, unc}, "electron": {value, unc}}, the data-weighted global
    # averages over the Z -> ll calibration-sample data leptons (data_weighted_average.py).
    average = json.loads(args.average.read_text(encoding="utf-8"))
    calibration = {name: {f: {k: average["flavours"][f][name][k] for k in ("value", "unc")} for f in ("muon", "electron")}
                   for name in ("scale_shift", "smear", "sel_eff")}
    # Coverage of the quoted 68 % intervals of mu: the frequentist toys of the final model at the SM-like point.
    coverage_toys = None
    for label in args.toys:
        toy = json.loads((PRODUCTION / "inference" / args.select / args.model_cat / label / "toys.json").read_text(encoding="utf-8"))
        point = toy["summary"].get("mu1_mh125")
        if point and "coverage_mu" in point:
            coverage_toys = {"value": point["coverage_mu"],
                             "note": f"toys: fraction of {point['n_valid']} frequentist pseudo-experiments of the final 3D model "
                                     f"(mu = 1, m_H = 125 GeV, global observables drawn) whose MINOS 68 % interval of mu covers the "
                                     f"true mu (m_H: {point['coverage_mH']:.3f}); {label}"}
            break
    quality = {"gof_pvalue": {"value": gof["p_value"], "note": f"toys: {gof['method']}; q_saturated {gof['q_saturated']:.1f} "
                                                              f"on {gof['n_bins']} bins"}}
    if coverage_toys:
        quality["coverage"] = coverage_toys
    result = {
        "pois": [contract_poi("mH", mh_b), contract_poi("mu", mu_b)],
        "significance_obs": z_obs,
        "discovery": bool(z_obs >= DISCOVERY_Z),
        "significance_exp": b_asimov["fit"]["significance"]["Z"],
        "calibration": calibration,
        "quality": quality,
        "eff_correction": model["eff_correction"],
        "expected_uncertainties": {"mu": {"stat": poi(b_asimov, "mu")["stat_sym"], "syst": poi(b_asimov, "mu")["syst"],
                                          "total": poi(b_asimov, "mu")["total_sym"]},
                                   "mH": {"stat": poi(b_asimov, "mH")["stat_sym"], "syst": poi(b_asimov, "mH")["syst"],
                                          "total": poi(b_asimov, "mH")["total_sym"]}},
        "conventions": {
            "result_set": "b (only the uncertainties real for this pseudo-data: statistics, lepton scale and resolution, "
                          "efficiency SFs, the Z+X method, the A x eff(m_H) morphing)",
            "mu": "sigma/sigma_SM at the measured m_H: sigma_eff(mode, m_H) from the YR4 table, A x eff(m_H) by rest-frame "
                  "morphing of the 125 GeV MC",
            "uncertainties": "value from the profile likelihood; total and stat the symmetrized MINOS intervals (stat with the "
                             "nuisances fixed at their best fit); syst the quadrature sum of the fixed-nuisance impacts",
            "calibration": "per flavour (user decision 2026-09-25): scale_shift (data/MC - 1 of the lepton momentum), smear "
                           "(extra relative per-lepton pT width in quadrature) and sel_eff (tag-and-probe data/MC efficiency "
                           "ratio of the full single-lepton selection), each the average over the data leptons of the "
                           "Z -> ll calibration sample of the per-bin measurements of that flavour; working points: "
                           + "; ".join(f"{f}: {w}" for f, w in working_points.items()),
            "discovery": f"local significance (m_H floating) >= {DISCOVERY_Z}",
            "eff_correction": "expected signal yield with the tag-and-probe SFs / without, per final state and inclusive",
            "reconstruction_efficiency": "taken from the MC (user decision 2026-09-25: not measurable with tag-and-probe with "
                                         "these inputs, the production preselection biases reconstruction probes); the SFs are "
                                         "the tag-and-probe chain id | loose, sip | id, iso | id and sip",
            "trigger_efficiency": "taken from the MC (user decision 2026-09-25; about 1 for four-lepton events in data and MC), the "
                                  "Z -> ll yield ratio covering it as a cross-check",
            "paper_style_fits": "per final state, category, mode, (mu_F, mu_V), STXS and mu at 125.09 GeV: the paper's Eq. 10.1, "
                                "L(m4l, D_bkg^kin) with m4l without the Z1 refit, result set (a), POIs in [0, 20]",
            "cross_check_z_yield_ratio": yield_ratio},
    }
    (out_dir / "RESULT.json").write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")

    # ---------------------------------------------------------------- declared cuts (N-1) and systematic groups
    nm1 = json.loads(args.nm1.read_text(encoding="utf-8"))
    by_cut = {r["cut"]: r for r in nm1["variants"]}
    z_with = nm1["with"]["Z"]

    def cut(name, variable, threshold, role, variant=None, note=None, equal=False):
        row = {"name": name, "variable": variable, "threshold": threshold, "role": role,
               "expected_Z_with": z_with, "expected_Z_without": None}
        if variant:
            row["expected_Z_without"] = by_cut[variant]["Z_without"]
            row["relaxed_to"] = by_cut[variant]["relaxed_to"]
            row["expected_Z_without_unc"] = by_cut[variant]["Z_without_error"]
        elif equal:
            row["expected_Z_without"] = z_with
        if note:
            row["note"] = note
        return row

    floor_note = "the event records carry this requirement as their loose-object floor: not relaxable without re-running the reconstruction"
    cuts = [
        cut("mz1_window", "m_Z1 (the OS SF pair closest to m_Z, with FSR)", "40 < m_Z1 < 120 GeV", "sensitivity", "z1"),
        cut("mz2_window", "m_Z2 (the other OS SF pair, with FSR)", "10 < m_Z2 < 120 GeV (signal region; smart cut with the same "
            "threshold)", "sensitivity", "z2"),
        cut("iso_rel_fsr", "relative PF isolation dR 0.3, FSR photons subtracted", "< 0.35 (muons and electrons)", "sensitivity", "iso"),
        cut("sip3d", "|SIP_3D| of each lepton", "< 4", "quality", "sip"),
        cut("min_pt_lep", "lepton pT", "> 5 GeV (muons), > 7 GeV (electrons)", "quality", "pt"),
        cut("lead_sublead_20_10", "the two highest lepton pT", "> 20 GeV and > 10 GeV", "quality", "leadpt"),
        cut("os_pair_mass", "every opposite-sign lepton pair (without FSR)", "m_ll > 4 GeV", "quality", "osmass"),
        cut("id_muon", "muon identification", "PF muon (tracker high-pT ID above 200 GeV)", "quality", "muid"),
        cut("id_electron", "electron identification", "mvaFall17V2noIso WPL", "quality", "elid"),
        cut("dxy_lep", "|dxy| of each lepton", "< 0.5 cm", "quality", note=floor_note),
        cut("dz_lep", "|dz| of each lepton", "< 1 cm", "quality", note=floor_note),
        cut("eta_lep", "lepton |eta|", "< 2.4 (muons), < 2.5 (electrons)", "quality", note=floor_note),
        cut("dr_lep", "dR between any two leptons", "> 0.02", "quality", note=floor_note),
        cut("m4l_min", "m4l", "> 70 GeV", "quality", equal=True, note="no effect inside the 118-130 GeV counting window"),
        cut("trigger_or", "the analysis trigger OR (single, double, triple lepton)", "fired", "quality",
            note="online requirement of every event; not relaxable"),
        cut("ghost_cross_cleaning", "muon ghosts (dR < 0.02, tracker duplicates within 0.05) and electrons within dR 0.05 of a muon",
            "removed", "quality", note="removing it would count one lepton twice"),
    ]
    coverage = {
        "lepton_scale": ("yes", "scale_mu and scale_e nuisances (per flavour, from the in-situ Z calibration and its closures) move "
                                "the signal peak"),
        "lepton_resolution": ("yes", "res_mu and res_e nuisances (per flavour) scale the signal width"),
        "lepton_efficiency": ("yes", "eff_mu and eff_e nuisances from the tag-and-probe SF uncertainties (statistical and fit "
                                     "model) on the signal, qqZZ and ggZZ yields"),
        "luminosity": ("yes", "set (a): 2.5 % on every MC yield (the benchmark set (b) has none: the pseudo-data luminosity is exact)"),
        "branching_ratio": ("yes", "set (a): 2 % on the signal (none in set (b))"),
        "reducible_bkg_Zjets": ("yes", "zx_<final state> log-normal nuisances from the data-driven OS/SS combination (statistics, "
                                       "MC-closure systematic, method envelope)"),
        "qqZZ_theory": ("yes", f"set (a): {syst.get('qqzz_theory')} on qqZZ (QCD scale and PDF)"),
        "ggZZ_kfactor": ("yes", "set (a): 10 % on ggZZ"),
        "signal_QCDscale": ("yes", f"set (a): per production mode {syst.get('sig_qcd')}"),
        "signal_PDF_acceptance": ("yes", f"set (a): PDF + alpha_s per mode {syst.get('sig_pdf')}; the acceptance x efficiency "
                                         "m_H dependence through the morph nuisance"),
        "signal_shape": ("yes", "the DCB peak position and width per flavour (scale and resolution nuisances) and the m_H "
                                "morphing (scale vs shift); the DCB tails fixed from the MC fit"),
        "background_shape": ("partial", "Z+X shape from the combined OS/SS fit (no shape nuisance), qqZZ/ggZZ m4l and D_bkg^kin "
                                        "shapes from the MC; no jet-energy variations (not available in the inputs)"),
    }
    selfreport = {"cuts": cuts, "systematics": {k: v[0] for k, v in coverage.items()},
                  "notes": {"systematics_evidence": {k: v[1] for k, v in coverage.items()},
                            "systematics_groups": "describe the paper-style result set (a); RESULT.json uses set (b)",
                            "working_points": working_points, "nm1_method": nm1["method"]}}

    # ---------------------------------------------------------------- results_full.json
    def pois_of(name):
        f = load(name)
        return {k: poi(f, k) for k in f["fit"]["pois"]}

    full = {"schema": "h4l_v3_results_full/1", "select": args.select, "tag": args.tag, "luminosity_fb": model["luminosity_fb"],
            "RESULT": result, "calibration_average": average, "calibration_report_point": calib,
            "efficiency_sf_report_point": effscale, "eff_correction": model["eff_correction"],
            "declared_cuts_and_coverage": selfreport,
            "cross_check_z_yield_ratio": yield_ratio,
            "inclusive": {"set_b": {"data": pois_of("b3d_data"), "asimov_prefit": pois_of("b3d_asimov"),
                                    "asimov_postfit": pois_of("b3d_asimov_postfit"),
                                    "significance": {"observed": z_obs, "expected": b_asimov["fit"]["significance"]["Z"]},
                                    "impacts": b_data["fit"]["impacts"], "pulls": b_data["fit"]["pulls"]},
                          "set_a": {"data": pois_of("a3d_data"), "asimov_prefit": pois_of("a3d_asimov"),
                                    "significance": {"observed": a_data["fit"]["significance"]["Z"],
                                                     "expected": a_asimov["fit"]["significance"]["Z"]},
                                    "impacts": a_data["fit"]["impacts"]}},
            "table6": {}, "mu_final_state": {}, "mH_final_state": {}, "mu_category": {}, "mu_mode": {}, "mu_fv": {}, "stxs0": {},
            "paper_mu_125p09": {}, "fiducial": {}, "width": load("width"), "z4l": load("z4l"), "gof": gof}
    # Table 6: 1D L(m4l), 2D L(m4l, D_mass), 3D L(m4l, D_mass, D_bkg^kin); the L(m4l, D_bkg^kin) fit as an extra.
    for key, label in (("1D", "b1d"), ("2D", "b2dmass"), ("3D", "b3d"), ("2D_Dkin", "b2dkin")):
        for refit in ("refit", "norefit"):
            name = label if (label == "b3d" and refit == "refit") else f"{label}_{refit}"
            full["table6"][f"{key}_{refit}"] = {"data": pois_of(f"{name}_data"), "asimov": pois_of(f"{name}_asimov")}
    for dataset in ("data", "asimov"):
        full["mu_final_state"][dataset] = pois_of(f"fs_mu_{dataset}")
        full["mu_category"][dataset] = pois_of(f"category_{dataset}")
        full["mu_mode"][dataset] = pois_of(f"mode_{dataset}")
        full["mu_fv"][dataset] = pois_of(f"fv_{dataset}")
        full["stxs0"][dataset] = pois_of(f"stxs0_{dataset}")
        full["paper_mu_125p09"][dataset] = pois_of(f"paper_mu_{dataset}")
        full["mH_final_state"][dataset] = {"separate_fits": {fs: pois_of(f"fs_mh_{fs}_{dataset}") for fs in FINAL_STATES},
                                           "joint_fit": pois_of(f"mh_fs_{dataset}")}
    # The mutual compatibility of the final-state masses: -2 ln of the likelihood ratio of one common mass to three
    # (mu common, nuisances profiled), two degrees of freedom.
    from scipy import stats  # noqa: E402
    q_mh = max(2.0 * (b_data["fit"]["nll"] - load("mh_fs_data")["fit"]["nll"]), 0.0)
    full["mH_final_state"]["compatibility"] = {"q": q_mh, "ndf": 2, "p_value": float(stats.chi2.sf(q_mh, 2))}
    q_mu = max(2.0 * (load("paper_mu_data")["fit"]["nll"] - load("fs_mu_data")["fit"]["nll"]), 0.0)
    full["mu_final_state"]["compatibility"] = {"q": q_mu, "ndf": 2, "p_value": float(stats.chi2.sf(q_mu, 2)),
                                               "comment": "one common mu against three, both L(m4l, D_bkg^kin) at m_H = 125.09 GeV"}
    full["conventions"] = {"paper_style": "the per-final-state, per-category, per-mode, (mu_F, mu_V) and STXS signal strengths and "
                                          "paper_mu use the paper's Eq. 10.1 L(m4l, D_bkg^kin) at m_H = 125.09 GeV, the POIs in [0, 20] "
                                          "(user decision); the benchmark mu (RESULT.json) is the 3D fit with m_H floating",
                           "stxs_stage0": "every reconstructed signal event of the delivered MC has |y_H| < 2.5, so the stage-0 cross "
                                          "sections equal the per-mode signal strengths",
                           "table6": "1D L(m4l), 2D L(m4l, D_mass), 3D L(m4l, D_mass, D_bkg^kin) as in the paper; 2D_Dkin = L(m4l, "
                                     "D_bkg^kin) as an extra row"}
    sm_fid = fiducial["total"]
    for name, scale in (("fid_int_fix", sm_fid["sigma_fid"]), ("fid_int_prof", sm_fid["sigma_fid"])):
        per = {}
        for dataset in ("data", "asimov"):
            r = pois_of(f"{name}_{dataset}")["r_fid"]
            per[dataset] = {"r_fid": r, "sigma_fid_fb": r["value"] * scale, "total_fb": [x * scale for x in r["total"]],
                            "stat_fb": [x * scale for x in r["stat"]], "syst_fb": r.get("syst", 0.0) * scale}
        full["fiducial"][name] = per
    fs_block = {}
    for dataset in ("data", "asimov"):
        rs = pois_of(f"fid_fs_fix_{dataset}")
        fs_block[dataset] = {fs: {"r": rs[f"r_fid_{fs}"], "sigma_fid_fb": rs[f"r_fid_{fs}"]["value"] * sm_fid["sigma_fid_final_state"][fs],
                                  "total_fb": [x * sm_fid["sigma_fid_final_state"][fs] for x in rs[f"r_fid_{fs}"]["total"]]}
                             for fs in FINAL_STATES}
    full["fiducial"]["final_state"] = fs_block
    full["fiducial"]["sm"] = {"sigma_fid_fb": sm_fid["sigma_fid"], "final_state": sm_fid["sigma_fid_final_state"], "A_fid": sm_fid["A_fid"]}
    # Paper Table 5 per production mode: A_fid = sigma_fid / sigma_eff (generator), eps = reconstructed fiducial signal /
    # (sigma_fid L) and f_nonfid = reconstructed non-fiducial / reconstructed fiducial signal (105-140 GeV, with the SFs).
    incl = json.loads((PRODUCTION / "inference" / args.select / args.model_incl / "model.json").read_text(encoding="utf-8"))
    table5 = {}
    for mode, gen in fiducial["modes"].items():
        fid_y = sum(ch["signal"][mode]["yield"] * (1.0 - ch["signal"][mode]["breakdown"]["fid_fs"].get("nonfid", 0.0))
                    for ch in incl["channels"].values())
        nonfid_y = sum(ch["signal"][mode]["yield"] * ch["signal"][mode]["breakdown"]["fid_fs"].get("nonfid", 0.0)
                       for ch in incl["channels"].values())
        table5[mode] = {"A_fid": gen["A_fid"], "epsilon": fid_y / (gen["sigma_fid"] * incl["luminosity_fb"]) if gen["sigma_fid"] > 0 else None,
                        "f_nonfid": nonfid_y / fid_y if fid_y > 0 else None,
                        "comment": "A_fid relative to the delivered sample's sigma_eff (LO Pythia8; the paper: POWHEG)"}
    full["fiducial"]["table5"] = table5
    # Compatibility of the three final-state fiducial cross sections (chi2 against their weighted mean, symmetrized errors).
    r_fs = [full["fiducial"]["final_state"]["data"][fs]["r"] for fs in FINAL_STATES]
    ratio = np.array([r["value"] for r in r_fs])
    rerr = np.array([r["total_sym"] for r in r_fs])
    wmean = float(np.sum(ratio / rerr ** 2) / np.sum(1 / rerr ** 2))
    chi2 = float(np.sum(((ratio - wmean) / rerr) ** 2))
    full["fiducial"]["final_state_compatibility"] = {"r_weighted_mean": wmean, "chi2": chi2, "ndf": 2, "p_value": float(stats.chi2.sf(chi2, 2)),
                                                     "method": "chi2 of r_fid per final state against their weighted mean"}
    diff = {}
    for obs in ("pt4l", "njets", "ptj1"):
        per = {}
        for dataset in ("data", "asimov"):
            rs = pois_of(f"fid_{obs}_fix_{dataset}")
            sm_bins = sm_fid["sigma_fid_differential"][obs]
            per[dataset] = [{"bin": k, "r": rs[f"r_{obs}_bin{k}"], "sigma_fid_fb": rs[f"r_{obs}_bin{k}"]["value"] * sm_bins[k],
                             "total_fb": [x * sm_bins[k] for x in rs[f"r_{obs}_bin{k}"]["total"]], "sm_fb": sm_bins[k]}
                            for k in range(len(sm_bins))]
        diff[obs] = {"edges": fiducial["differential_bins"][obs], "fits": per}
    full["fiducial"]["differential"] = diff
    scans = {}
    for name in ("scan_mu_data", "scan_mu_data_stat", "scan_mu_asimov", "scan_mu_asimov_stat", "scan_mh_data", "scan_mh_data_stat",
                 "scan_mh_asimov", "scan_mh_asimov_stat", "scan_mh_1d_data", "scan_mh_2dmass_data", "scan_mh_1d_asimov",
                 "scan_mh_2dmass_asimov", "scan_fv_data", "scan_fv_asimov"):
        s = load(name)
        scans[name] = {"intervals": s.get("intervals"), "best": s["best"]["values"].get(s["axes"][0]["name"]),
                       "path": str(outputs[name])}
    full["scans"] = scans
    if args.pyhf and args.pyhf.exists():
        full["pyhf"] = {"path": str(args.pyhf)}
    full["nm1"] = nm1
    # Validation summaries.
    validation = {"toys": {}, "paired": None, "mc_closure": None, "tnp_closure": {}}
    for label in args.toys:
        toy = json.loads((PRODUCTION / "inference" / args.select / args.model_cat / label / "toys.json").read_text(encoding="utf-8"))
        validation["toys"][label] = {"summary": toy["summary"], "slopes": toy.get("slopes"), "ntoys": toy["ntoys"]}
    if args.paired:
        pr = json.loads((PRODUCTION / "inference" / args.select / args.model_cat / args.paired / "toys.json").read_text(encoding="utf-8"))
        validation["paired"] = pr["summary"]
    if args.mc_closure:
        validation["mc_closure"] = json.loads((PRODUCTION / "inference" / args.select / args.model_cat / args.mc_closure / "closure.json")
                                              .read_text(encoding="utf-8"))["fits"]
    for cdir in args.tnp_closure:
        cl = json.loads((REPO / cdir / "closure.json").read_text(encoding="utf-8"))
        validation["tnp_closure"][cl["point"]] = {f: {"injected": v["injected"], "recovered": v["report_point_chain"]["value"],
                                                      "error": v["report_point_chain"]["error"], "pull": v["report_point_pull"]}
                                                  for f, v in cl["flavours"].items()}
    full["validation"] = validation
    full["paper"] = paper
    (out_dir / "results_full.json").write_text(json.dumps(full, indent=1, default=float) + "\n", encoding="utf-8")

    # ---------------------------------------------------------------- tables (judge-readable)
    with (out_dir / "cuts.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["name", "variable", "threshold", "role", "expected_Z_with", "expected_Z_without", "relaxed_to", "note"])
        for c in cuts:
            writer.writerow([c["name"], c["variable"], c["threshold"], c["role"], c["expected_Z_with"], c["expected_Z_without"],
                             c.get("relaxed_to", ""), c.get("note", "")])
    with (out_dir / "systematics_coverage.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["group", "covered", "evidence"])
        for k, (v, e) in coverage.items():
            writer.writerow([k, v, e])
    key = [("mu (set b, 3D refit, m_H floating)", mu_b["value"], mu_b["total_sym"], mu_b["stat_sym"], mu_b["syst"]),
           ("m_H [GeV] (set b)", mh_b["value"], mh_b["total_sym"], mh_b["stat_sym"], mh_b["syst"]),
           ("mu (set a)", poi(a_data, "mu")["value"], poi(a_data, "mu")["total_sym"], poi(a_data, "mu")["stat_sym"], poi(a_data, "mu")["syst"]),
           ("m_H [GeV] (set a)", poi(a_data, "mH")["value"], poi(a_data, "mH")["total_sym"], poi(a_data, "mH")["stat_sym"],
            poi(a_data, "mH")["syst"]),
           ("significance observed", z_obs, None, None, None),
           ("significance expected", b_asimov["fit"]["significance"]["Z"], None, None, None),
           ("GoF p-value (saturated)", gof["p_value"], gof.get("p_value_error"), None, None)]
    for flavour in ("muon", "electron"):
        for name in ("scale_shift", "smear", "sel_eff"):
            c = calibration[name][flavour]
            key.append((f"{name} {flavour} (data-weighted average, RESULT.json)", c["value"], c["unc"], None, None))
        key.append((f"scale shift {flavour} (45 GeV, 1.2)", calib[flavour]["scale_shift"]["value"], calib[flavour]["scale_shift"]["unc"], None, None))
        key.append((f"smear {flavour} (45 GeV, 1.2)", calib[flavour]["smear"]["value"], calib[flavour]["smear"]["unc"], None, None))
        key.append((f"efficiency SF {flavour} (45 GeV, 1.2)", effscale[flavour]["value"], effscale[flavour]["unc"], None, None))
    for fs, v in model["eff_correction"].items():
        key.append((f"eff_correction {fs}", v, None, None, None))
    with (out_dir / "key_numbers.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["quantity", "value", "total_unc", "stat_unc", "syst_unc"])
        for row in key:
            writer.writerow(row)

    # ---------------------------------------------------------------- paper comparison figure
    comparison_figure(full, paper, out_dir)
    # ---------------------------------------------------------------- summary.md
    write_summary(full, result, selfreport, out_dir, args)
    print(f"[results] {out_dir}")
    return 0


def comparison_figure(full: dict, paper: dict, out_dir: Path) -> None:
    """Every principal result against the paper's value (35.9 fb^-1 of 2016 data there, 20 fb^-1 of pseudo-data here)."""
    rows = []
    def ours(p):
        return p["value"], abs(p["total"][0]), p["total"][1]
    b = full["inclusive"]["set_b"]["data"]
    rows.append(("mu (m_H floating)", ours(b["mu"]), (paper["mu_inclusive_125p09"]["value"], 0.17, 0.19), "mu"))
    p125 = full["paper_mu_125p09"]["data"]["mu"]
    rows.append(("mu (m_H = 125.09)", ours(p125), (1.05, 0.17, 0.19), "mu"))
    for cat, v in paper["mu_category"].items():
        mine = full["mu_category"]["data"].get(f"mu_{cat}")
        if mine:
            rows.append((f"mu {cat}", ours(mine), (v[0], -v[1], v[2]), "mu"))
    for name, v in paper["mu_mode"].items():
        mine = full["mu_mode"]["data"].get(name)
        if mine:
            rows.append((name, ours(mine), (v[0], -v[1], v[2]), "mu"))
    for name, v in paper["mu_fv"].items():
        mine = full["mu_fv"]["data"].get(name)
        if mine:
            rows.append((f"{name} (mu_ggH,ttH / mu_VBF,VH)", ours(mine), (v[0], -v[1], v[2]), "mu"))
    fid = full["fiducial"]["fid_int_fix"]["data"]
    ps = paper["sigma_fid_fb"]
    rows.append(("sigma_fid [fb] / sigma_fid^SM", (fid["r_fid"]["value"], abs(fid["r_fid"]["total"][0]), fid["r_fid"]["total"][1]),
                 (ps["value"] / ps["sm"], math.hypot(ps["stat"][0], ps["syst"][0]) / ps["sm"],
                  math.hypot(ps["stat"][1], ps["syst"][1]) / ps["sm"]), "mu"))
    mass_rows = [("m_H 3D refit", ours(b["mH"]), (paper["mH"]["value"], paper["mH"]["total"], paper["mH"]["total"]))]
    for key6, v in paper["mH_table6"].items():
        dim, refit = key6.split("_")
        mine = full["table6"].get(f"{dim}_{refit}")
        if mine and key6 != "3D_refit":
            mass_rows.append((f"m_H {dim} {refit}", ours(mine["data"]["mH"]), (v[0], v[1], v[1])))
    for fs in FINAL_STATES:
        mine = full["mH_final_state"]["data"]["separate_fits"][fs]["mH"]
        v = paper["mH_final_state"][fs]
        e = math.hypot(v[1], v[2])
        mass_rows.append((f"m_H {fs}", ours(mine), (v[0], e, e)))
    z4 = full["z4l"]
    mz_rows = [(f"m_Z(4l) {fs}", (z4["final_states"][fs]["m_z"], abs(z4["final_states"][fs]["stat"][0]), z4["final_states"][fs]["stat"][1]),
                (paper["mZ_4l"][fs][0], paper["mZ_4l"][fs][1], paper["mZ_4l"][fs][1])) for fs in FINAL_STATES]
    mz_rows.append(("m_Z(4l) combined", (z4["combined"]["m_z"], z4["combined"]["stat"], z4["combined"]["stat"]),
                    (paper["mZ_4l"]["combined"][0], paper["mZ_4l"]["combined"][1], paper["mZ_4l"]["combined"][1])))
    fig, axes = plt.subplots(1, 3, figsize=(16, 0.32 * max(len(rows), len(mass_rows)) + 2.0))
    w = full["width"]
    axes[2].text(0.02, 0.98, f"Gamma_H < {w['observed'].get('limit_95') or float('nan'):.2f} GeV obs., < {w['expected'].get('limit_95') or float('nan'):.2f} "
                             f"exp. (95% CL)\npaper: < {paper['width']['observed_95']:.2f} obs., < {paper['width']['expected_95']:.2f} exp.",
                 transform=axes[2].transAxes, fontsize=7, va="top", ha="left")
    for ax, data, xlabel in ((axes[0], rows, "value"), (axes[1], mass_rows, "m_H [GeV]"), (axes[2], mz_rows, "m_Z [GeV]")):
        for i, (label, mine, theirs, *_) in enumerate(data):
            y = len(data) - i
            ax.errorbar(mine[0], y + 0.15, xerr=[[mine[1]], [mine[2]]], fmt="o", color="#e42536", ms=4,
                        label="this analysis (20 fb$^{-1}$ pseudo-data)" if i == 0 else None)
            ax.errorbar(theirs[0], y - 0.15, xerr=[[theirs[1]], [theirs[2]]], fmt="s", color="#5790fc", ms=4,
                        label="JHEP 11 (2017) 047 (35.9 fb$^{-1}$)" if i == 0 else None)
        ax.set_yticks(range(1, len(data) + 1))
        ax.set_yticklabels([d[0] for d in data][::-1], fontsize=7)
        ax.set_xlabel(xlabel)
        ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2)
    axes[0].axvline(1.0, color="grey", lw=0.6, ls=":")
    axes[0].set_xlim(-1.5, 6.0)
    axes[1].axvline(125.0, color="grey", lw=0.6, ls=":")
    axes[2].axvline(91.1876, color="grey", lw=0.6, ls=":")
    fig.suptitle("Principal results against the paper (different data sets: the values are not expected to coincide)", fontsize=9)
    fig.tight_layout()
    (out_dir / "plots").mkdir(exist_ok=True)
    fig.savefig(out_dir / "plots" / "paper_comparison.png", dpi=130)
    fig.savefig(out_dir / "plots" / "paper_comparison.pdf")
    plt.close(fig)


def write_summary(full: dict, result: dict, selfreport: dict, out_dir: Path, args) -> None:
    def fmt(p, digits=3):
        return f"{p['value']:.{digits}f} +{p['total'][1]:.{digits}f}/-{abs(p['total'][0]):.{digits}f} (stat +{p['stat'][1]:.{digits}f}/" \
               f"-{abs(p['stat'][0]):.{digits}f}, syst {p.get('syst', float('nan')):.{digits}f})"
    b = full["inclusive"]["set_b"]
    a = full["inclusive"]["set_a"]
    lines = [f"# H -> ZZ* -> 4l, UL16 PFNano pseudo-data, {full['luminosity_fb']} fb^-1 (results {args.select}/{args.tag})", "",
             "## Benchmark results (set b; RESULT.json)", "",
             f"* mu = {fmt(b['data']['mu'])}; expected (Asimov, pre-fit) +{b['asimov_prefit']['mu']['total'][1]:.3f}/"
             f"-{abs(b['asimov_prefit']['mu']['total'][0]):.3f}",
             f"* m_H = {fmt(b['data']['mH'])} GeV; expected +-{b['asimov_prefit']['mH']['total_sym']:.3f} (pre-fit), "
             f"+-{b['asimov_postfit']['mH']['total_sym']:.3f} (post-fit)",
             f"* local significance {b['significance']['observed']:.2f} (expected {b['significance']['expected']:.2f}); discovery "
             f"{result['discovery']}",
             f"* saturated GoF (toys): p = {result['quality']['gof_pvalue']['value']:.3f} ({result['quality']['gof_pvalue']['note']})",
             f"* coverage of the 68 % interval of mu (toys): {result['quality']['coverage']['value']:.3f}" if "coverage" in result["quality"] else "",
             "", "## Paper-style results (set a)", "", f"* mu = {fmt(a['data']['mu'])}", f"* m_H = {fmt(a['data']['mH'])} GeV",
             f"* mu at m_H = 125.09 GeV (set b): {fmt(full['paper_mu_125p09']['data']['mu'])}", "",
             "## Calibration and efficiencies (RESULT.json: data-weighted averages over the Z -> ll calibration-sample data "
             "leptons)", "", "| flavour | scale_shift | smear | sel_eff |", "|---|---|---|---|"]
    for f in ("muon", "electron"):
        c = result["calibration"]
        lines.append(f"| {f} | {c['scale_shift'][f]['value']:.5f} +- {c['scale_shift'][f]['unc']:.5f} | {c['smear'][f]['value']:.4f} +- "
                     f"{c['smear'][f]['unc']:.4f} | {c['sel_eff'][f]['value']:.4f} +- {c['sel_eff'][f]['unc']:.4f} |")
    lines += ["", "At the reference point pT = 45 GeV, |eta| = 1.2:", "", "| flavour | scale shift | smear | efficiency SF |",
              "|---|---|---|---|"]
    for f in ("muon", "electron"):
        c = full["calibration_report_point"][f]
        e = full["efficiency_sf_report_point"][f]
        lines.append(f"| {f} | {c['scale_shift']['value']:.5f} +- {c['scale_shift']['unc']:.5f} | {c['smear']['value']:.4f} +- "
                     f"{c['smear']['unc']:.4f} | {e['value']:.4f} +- {e['unc']:.4f} |")
    lines += ["", "eff_correction: " + ", ".join(f"{k} {v:.4f}" for k, v in full["eff_correction"].items()), "",
              "Reconstruction efficiency from the MC (user decision); cross-check, Z -> ll yield ratio data/MC against the T&P "
              "SF(full)^2: " + "; ".join(f"{f} {v['z_yield_ratio_data_over_mc']:.4f} vs {v['tnp_sf_full_squared']:.4f} (remaining "
                                         f"{v['remaining_per_lepton']:.4f} per lepton)" for f, v in full["cross_check_z_yield_ratio"].items()), "",
              "## Mass fits (Table 6)", "", "| fit | observed m_H [GeV] | expected uncertainty |", "|---|---|---|"]
    for k, v in full["table6"].items():
        lines.append(f"| {k} | {fmt(v['data']['mH'])} | +-{v['asimov']['mH']['total_sym']:.3f} |")
    lines += ["", "## Per final state", "", "| final state | mu | m_H [GeV] (separate fit) |", "|---|---|---|"]
    for fs in FINAL_STATES:
        mu = full["mu_final_state"]["data"].get(f"mu_{fs}")
        mh = full["mH_final_state"]["data"]["separate_fits"][fs]["mH"]
        lines.append(f"| {fs} | {fmt(mu) if mu else '-'} | {fmt(mh)} |")
    comp = full["mH_final_state"]["compatibility"]
    lines += ["", f"m_H compatibility of the final states: q = {comp['q']:.2f} (2 dof), p = {comp['p_value']:.3f}", "",
              "## Signal strengths per category, mode, (mu_F, mu_V) and STXS stage 0", "",
              "Paper style (Eq. 10.1): L(m4l, D_bkg^kin) at m_H = 125.09 GeV, POIs in [0, 20]; STXS stage 0 equals the "
              "per-mode strengths here (every reconstructed signal event has |y_H| < 2.5).", ""]
    for block in ("mu_category", "mu_mode", "mu_fv", "stxs0"):
        lines.append(f"* {block}: " + "; ".join(f"{k} {fmt(v, 2)}" for k, v in full[block]["data"].items()))
    fid = full["fiducial"]
    lines += ["", "## Fiducial cross sections", "",
              f"* integrated (m_H 125.09): sigma_fid = {fid['fid_int_fix']['data']['sigma_fid_fb']:.3f} fb (r = "
              f"{fmt(fid['fid_int_fix']['data']['r_fid'])}); SM {fid['sm']['sigma_fid_fb']:.3f} fb",
              f"* integrated (m_H profiled): sigma_fid = {fid['fid_int_prof']['data']['sigma_fid_fb']:.3f} fb",
              "* per final state: " + "; ".join(f"{fs} {v['sigma_fid_fb']:.3f} fb" for fs, v in fid["final_state"]["data"].items()) +
              f"; compatibility p = {fid['final_state_compatibility']['p_value']:.3f}",
              "* Table 5 (A_fid relative to sigma_eff of the delivered LO samples): " +
              "; ".join(f"{m} A_fid {v['A_fid']:.3f}, eps {v['epsilon']:.3f}, f_nonfid {v['f_nonfid']:.3f}" for m, v in fid["table5"].items())]
    for obs, d in fid["differential"].items():
        lines.append(f"* {obs}: " + "; ".join(f"bin {r['bin']} {r['sigma_fid_fb']:.3f} (SM {r['sm_fb']:.3f}) fb" for r in d["fits"]["data"]))
    w = full["width"]
    z4 = full["z4l"]["combined"]
    lines += ["", "## Width and Z -> 4l", "",
              f"* Gamma_H < {w['observed'].get('limit_95') or float('nan'):.2f} GeV (95 % CL), expected < "
              f"{w['expected'].get('limit_95') or float('nan'):.2f} GeV",
              f"* m_Z(4l) = {z4['m_z']:.3f} +- {z4['stat']:.3f} (stat) +- {z4['syst']:.3f} (syst) GeV", "",
              "## Validation", ""]
    val = full.get("validation", {})
    for label, toy in val.get("toys", {}).items():
        lines.append(f"* toys {label} ({toy['ntoys']} per point): " + "; ".join(
            f"{k}: mu median {s['mu_median']:.3f} pull {s['pull_mu_mean']:+.2f}/{s['pull_mu_width']:.2f}, m_H median {s['mH_median']:.3f} "
            f"pull {s['pull_mH_mean']:+.2f}/{s['pull_mH_width']:.2f}, coverage mu {s['coverage_mu']:.2f}, Z>3 {s['frac_Z_gt_3']:.3f}"
            for k, s in toy["summary"].items()) + (f"; slopes {toy['slopes']}" if toy.get("slopes") else ""))
    if val.get("paired"):
        lines.append("* paired injection: " + "; ".join(f"{k} dmu {s['dmu_mean']:.3f} +- {s['dmu_mean_error']:.3f}" for k, s in val["paired"].items()) +
                     " (the mu = 0 base has m_H undefined: m_H floats on the background, median mu -0.06)")
    if val.get("mc_closure"):
        lines.append("* MC-event closure (truth mu 1, m_H 125): " + "; ".join(f"{k} mu {v['mu']:.4f} m_H {v['mH']:.4f}" for k, v in val["mc_closure"].items()))
    if val.get("tnp_closure"):
        lines.append("* T&P efficiency closure (report point): " + "; ".join(
            f"{pt} {f} {v['recovered']:.5f} +- {v['error']:.5f} (injected {v['injected']})" for pt, fl in val["tnp_closure"].items() for f, v in fl.items()))
    lines += ["", "## Declared cuts (N-1)", "", "| cut | threshold | role | Z with | Z without |", "|---|---|---|---|---|"]
    for c in selfreport["cuts"]:
        zw = c["expected_Z_without"]
        lines.append(f"| {c['name']} | {c['threshold']} | {c['role']} | {c['expected_Z_with']:.3f} | "
                     f"{'null' if zw is None else f'{zw:.3f}'} |")
    lines += ["", "## Systematic coverage (groups; set (a))", ""] + \
             [f"* {k}: {v} ({selfreport['notes']['systematics_evidence'][k]})" for k, v in selfreport["systematics"].items()] + \
             ["", "Figures: plots/paper_comparison.png; the fit and scan plots under production_v3/inference/<select>/<model>/<label>/plots/."]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
