"""N-1 expected significances of the declared selection cuts (eval_selfreport.json "cuts").

    pixi run py -- analysis_v3/reconstruction/scripts/nm1_significance.py --reference v5_nm1d_ref --nominal v5 \
        --prefix v5_nm1d --variants iso sip z1 z2 pt leadpt osmass muid elid --zx-label zx_nm1 --label nm1_v1

User decisions (2026-09-25): "without" a cut = relaxed to the floor of the event records (only the
isolation and the identifications are removed outright; make_nm1_configs.py), the reducible background
of every selection its own data-driven Z+X (zx_estimate.py on that selection, fake-rate numerator and
denominator following it), signal and ZZ from the MC; one frozen program for the reference and every
N-1 selection, the reference checked row by row against the nominal selection (v5).
For the reference and every N-1 selection, in the counting window 118 < m4l < 130 GeV (m4l without the
Z1 refit, as in the paper's Table 2):
    S = ggH + VBF + VH MC (m_H = 125 GeV),  B = qqZZ + ggZZ MC + Z+X (combined OS/SS yield of each final
    state x its fraction in the window, from the rows),
no scale factors, and the expected significance of the counting experiment
    Z = sqrt(2 [(S + B) ln(1 + S/B) - S]).
The uncertainty of Z propagates the MC statistics of S and B and the total Z+X uncertainty (recorded).
Every N-1 selection must come from the reference's frozen program, event records, payloads, refit inputs and
MELA block and its configuration differs only in the relaxed requirement (checked).  Writes
production_v3/optimization/<nominal>/<label>/nm1.json, nm1.csv and nm1.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_select_io as io  # noqa: E402

SIGNAL = ["GluGluToHToZZ_M125", "VBF_HToZZ_M125", "VHToZZ_M125"]
MC_BACKGROUNDS = {"qqZZ": ["ZZTo4L"], "ggZZ": ["GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"]}
FINAL_STATES = ("4mu", "4e", "2e2mu")
NOTES = {"iso": "FSR-subtracted isolation < 0.35 removed",
         "sip": "SIP < 4 relaxed to SIP < 8 (event-record floor)",
         "z1": "40 < m_Z1 < 120 GeV relaxed to 30 < m_Z1 < 130 GeV (event-record floor)",
         "z2": "m_Z2 > 10 GeV (and the smart cut's threshold), m_Z2 < 120 GeV relaxed to 3-130 GeV (event-record floor; "
               "effectively 4 GeV through the kept m_ll(OS) > 4 GeV)",
         "pt": "lepton pT > 5 (mu) / 7 (e) GeV relaxed to 3 / 5 GeV (event-record floor)",
         "leadpt": "leading / subleading lepton pT > 20 / 10 GeV relaxed to 15 / 7 GeV (event-record floor)",
         "osmass": "m_ll(OS) > 4 GeV relaxed to 2 GeV (event-record floor)",
         "muid": "muon identification (PF, or tracker high-pT) removed",
         "elid": "electron identification (mvaFall17V2noIso WPL) removed"}
# The configuration keys each variant may change (make_nm1_configs.py), besides version and description.
CHANGED = {"iso": [("leptons", "max_iso_fsr")], "sip": [("leptons", "max_sip")], "z1": [("candidate", "z1_mass")],
           "z2": [("candidate", "z2_high"), ("candidate", "regions")], "pt": [("leptons", "muon"), ("leptons", "electron")],
           "leadpt": [("candidate", "leading_pt"), ("candidate", "subleading_pt")], "osmass": [("candidate", "min_os_mass")],
           "muid": [("leptons", "muon")], "elid": [("leptons", "electron")]}
COMPARED = ["final_state", "m4l", "m4l_refit", "mz1", "mz2", "weight", "d_bkg_kin"]


def plan_of(scan: dict) -> dict:
    return json.loads(Path(scan["plan"]).read_text(encoding="utf-8"))


def check_variant(reference: str, variant: str, cut: str) -> None:
    """The N-1 selection differs from the reference only in its relaxed requirement."""
    scans = [io.load_scan(v) for v in (reference, variant)]
    plans = [plan_of(s) for s in scans]
    if scans[0]["program_sha256"] != scans[1]["program_sha256"]:
        raise SystemExit(f"{variant}: another selection program than {reference}")
    if scans[0]["manifests"] != scans[1]["manifests"]:
        raise SystemExit(f"{variant}: other manifests than {reference}")
    for key in ("reco_scan_sha256", "mela"):
        if json.dumps(plans[0][key], sort_keys=True) != json.dumps(plans[1][key], sort_keys=True):
            raise SystemExit(f"{variant}: the plan's {key} differs from {reference}")
    for key in ("calibration_payload", "lambda_payload", "refit_inputs"):
        if plans[0][key]["sha256"] != plans[1][key]["sha256"]:
            raise SystemExit(f"{variant}: the {key} differs from {reference}")
    for name, info in scans[0]["samples"].items():
        other = scans[1]["samples"].get(name)
        if other is None or any(info.get(k) != other.get(k) for k in ("kind", "files", "genEventSumw", "sigma_eff_pb", "lumi_fb")):
            raise SystemExit(f"{variant}: the sample {name} differs from {reference}")
    configs = [plan["tasks"][0]["selection_config"] for plan in plans]
    allowed = {("version",), ("description",)} | set(CHANGED[cut])
    differences = set()
    for section in set(configs[0]) | set(configs[1]):
        a, b = configs[0].get(section), configs[1].get(section)
        if isinstance(a, dict) and isinstance(b, dict):
            differences |= {(section, k) for k in set(a) | set(b) if a.get(k) != b.get(k)}
        elif a != b:
            differences.add((section,))
    if not differences <= allowed or not differences & set(CHANGED[cut]):
        raise SystemExit(f"{variant}: configuration differences {sorted(differences)} are not those of the {cut} variant")
    if {t["task_id"] for t in plans[0]["tasks"]} != {t["task_id"] for t in plans[1]["tasks"]}:
        raise SystemExit(f"{variant}: its tasks are not the reference's")


def check_reference(reference: str, nominal: str) -> dict:
    """The reference (the nominal configuration with the program of the N-1 selections) reproduces the nominal
    selection: the same configuration, inputs and payloads, and identical signal-region rows of every sample."""
    scans = [io.load_scan(v) for v in (nominal, reference)]
    plans = [plan_of(s) for s in scans]
    if plans[0]["selection_config"]["sha256"] != plans[1]["selection_config"]["sha256"]:
        raise SystemExit(f"{reference}: another selection configuration than {nominal}")
    for key in ("calibration_payload", "lambda_payload", "refit_inputs"):
        if plans[0][key]["sha256"] != plans[1][key]["sha256"]:
            raise SystemExit(f"{reference}: the {key} differs from {nominal}")
    rows = {}
    for sample in sorted(scans[0]["samples"]):
        arrays = [io.read(s, [sample], "SR", COMPARED) for s in scans]
        if len(arrays[0]["weight"]) != len(arrays[1]["weight"]):
            raise SystemExit(f"{reference}: {sample} has {len(arrays[1]['weight'])} SR rows, {nominal} {len(arrays[0]['weight'])}")
        for name in COMPARED:
            if not np.array_equal(arrays[0][name], arrays[1][name]):
                raise SystemExit(f"{reference}: {sample} SR column {name} differs from {nominal}")
        rows[sample] = int(len(arrays[0]["weight"]))
    return {"nominal": nominal, "reference": reference, "identical_sr_rows": rows,
            "programs": [scans[0]["program_sha256"], scans[1]["program_sha256"]]}


def asimov_z(s: float, b: float) -> float:
    if s <= 0 or b <= 0:
        return 0.0
    return math.sqrt(max(2.0 * ((s + b) * math.log(1.0 + s / b) - s), 0.0))


def expectation(version: str, zx_label: str, window: tuple) -> dict:
    scan = io.load_scan(version)
    cut = lambda a: (a["m4l"] > window[0]) & (a["m4l"] < window[1])  # noqa: E731
    out = {"selection": version, "selection_version": scan["selection_version"], "program_sha256": scan["program_sha256"],
           "components": {}}
    r = io.read(scan, SIGNAL, "SR", ["m4l"], cut=cut)
    out["components"]["signal"] = {"yield": float(np.sum(r["w"])), "mc_stat": float(np.sqrt(np.sum(r["w"] ** 2))), "entries": int(len(r["w"]))}
    for name, samples in MC_BACKGROUNDS.items():
        r = io.read(scan, samples, "SR", ["m4l"], cut=cut)
        out["components"][name] = {"yield": float(np.sum(r["w"])), "mc_stat": float(np.sqrt(np.sum(r["w"] ** 2))),
                                   "entries": int(len(r["w"]))}
    zx = json.loads((PRODUCTION / "backgrounds" / version / zx_label / "zx.json").read_text(encoding="utf-8"))
    if zx["select_scan_plan_sha256"] != scan["plan_sha256"]:
        raise SystemExit(f"{version}/{zx_label}: the Z+X estimate belongs to another plan")
    zx_yield, zx_var = 0.0, 0.0
    per_fs = {}
    for fs in FINAL_STATES:
        c = zx["combined"][fs]
        if list(c["table_window"]) != list(window):
            raise SystemExit(f"{version}/{zx_label}: its table window {c['table_window']} is not {window}")
        value = c["value"] * c["table_window_fraction"]
        # The envelope half-width of the combination as the uncertainty (symmetrized), scaled to the window.
        error = 0.5 * (c["envelope"][1] - c["envelope"][0]) * c["table_window_fraction"]
        per_fs[fs] = {"yield": value, "error": error}
        zx_yield += value
        zx_var += error ** 2
    out["components"]["zx"] = {"yield": zx_yield, "error": math.sqrt(zx_var), "per_final_state": per_fs,
                               "estimate": str(PRODUCTION / "backgrounds" / version / zx_label / "zx.json")}
    c = out["components"]
    s = c["signal"]["yield"]
    b = c["qqZZ"]["yield"] + c["ggZZ"]["yield"] + c["zx"]["yield"]
    var_s = c["signal"]["mc_stat"] ** 2
    var_b = c["qqZZ"]["mc_stat"] ** 2 + c["ggZZ"]["mc_stat"] ** 2 + c["zx"]["error"] ** 2
    z = asimov_z(s, b)
    # Linear propagation: dZ/dS and dZ/dB numerically.
    eps_s, eps_b = 1e-4 * s, 1e-4 * b
    dzds = (asimov_z(s + eps_s, b) - asimov_z(s - eps_s, b)) / (2 * eps_s)
    dzdb = (asimov_z(s, b + eps_b) - asimov_z(s, b - eps_b)) / (2 * eps_b)
    out.update({"S": s, "B": b, "S_error": math.sqrt(var_s), "B_error": math.sqrt(var_b), "Z": z,
                "Z_error": math.sqrt(dzds ** 2 * var_s + dzdb ** 2 * var_b),
                "Z_irreducible_only": asimov_z(s, c["qqZZ"]["yield"] + c["ggZZ"]["yield"])})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="the nominal configuration run with the N-1 program")
    parser.add_argument("--nominal", required=True, help="the nominal selection the reference must reproduce")
    parser.add_argument("--prefix", required=True, help="the N-1 selections are <prefix>_<cut>")
    parser.add_argument("--variants", nargs="+", required=True)
    parser.add_argument("--zx-label", required=True, help="the Z+X label of every selection")
    parser.add_argument("--label", required=True)
    parser.add_argument("--window", type=float, nargs=2, default=[118.0, 130.0])
    args = parser.parse_args()
    out_dir = PRODUCTION / "optimization" / args.nominal / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    window = tuple(args.window)
    reference_check = check_reference(args.reference, args.nominal)
    with_cut = expectation(args.reference, args.zx_label, window)
    rows = []
    for cut in args.variants:
        version = f"{args.prefix}_{cut}"
        check_variant(args.reference, version, cut)
        v = expectation(version, args.zx_label, window)
        rows.append({"cut": cut, "relaxed_to": NOTES.get(cut, ""), "Z_with": with_cut["Z"], "Z_with_error": with_cut["Z_error"],
                     "Z_without": v["Z"], "Z_without_error": v["Z_error"], "delta_Z": with_cut["Z"] - v["Z"],
                     "S_without": v["S"], "B_without": v["B"], "B_without_error": v["B_error"],
                     "Z_irreducible_only_without": v["Z_irreducible_only"], "selection": v})
    report = {"schema": "h4l_v3_nm1/2", "method": "counting in 118 < m4l < 130 GeV: S = ggH + VBF + VH MC, B = qqZZ + ggZZ MC + "
              "the data-driven Z+X of the same selection; no SFs; Z = sqrt(2 [(S + B) ln(1 + S/B) - S]); N-1 = one "
              "requirement relaxed to the event-record floor (user decisions 2026-09-25)",
              "window": list(window), "reference_check": reference_check, "with": with_cut, "variants": rows}
    out_dir.mkdir(parents=True)
    (out_dir / "nm1.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    with (out_dir / "nm1.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cut", "relaxed_to", "Z_with", "Z_with_error", "Z_without", "Z_without_error", "delta_Z", "S_without",
                         "B_without", "B_without_error", "Z_irreducible_only_without"])
        for r in rows:
            writer.writerow([r["cut"], r["relaxed_to"], f"{r['Z_with']:.4f}", f"{r['Z_with_error']:.4f}", f"{r['Z_without']:.4f}",
                             f"{r['Z_without_error']:.4f}", f"{r['delta_Z']:+.4f}", f"{r['S_without']:.3f}", f"{r['B_without']:.3f}",
                             f"{r['B_without_error']:.3f}", f"{r['Z_irreducible_only_without']:.4f}"])
    c = with_cut["components"]
    lines = [f"# N-1 expected significances ({window[0]:g} < m4l < {window[1]:g} GeV)", "", report["method"], "",
             f"With every cut ({args.reference}, identical to {args.nominal}): S = {with_cut['S']:.2f}, B = {with_cut['B']:.2f} +- "
             f"{with_cut['B_error']:.2f} (qqZZ {c['qqZZ']['yield']:.2f}, ggZZ {c['ggZZ']['yield']:.2f}, Z+X {c['zx']['yield']:.2f} +- "
             f"{c['zx']['error']:.2f}), Z = {with_cut['Z']:.3f} +- {with_cut['Z_error']:.3f}", "",
             "| cut | relaxed to | Z with | Z without | delta Z | S without | B without |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['cut']} | {r['relaxed_to']} | {r['Z_with']:.3f} | {r['Z_without']:.3f} +- {r['Z_without_error']:.3f} | "
                     f"{r['delta_Z']:+.3f} | {r['S_without']:.2f} | {r['B_without']:.2f} +- {r['B_without_error']:.2f} |")
    (out_dir / "nm1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[4:]))
    print(f"[nm1] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
