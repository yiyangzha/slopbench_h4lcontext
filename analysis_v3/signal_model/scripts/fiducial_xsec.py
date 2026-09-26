"""Generator-level SM cross sections of the fiducial and stage-0 measurements (stage 6).

    pixi run py -- analysis_v3/signal_model/scripts/fiducial_xsec.py --select v5 --label gen_v1

For every signal sample (ggH, VBF, VH at m_H = 125 GeV) the skim GenTable (one row per generated
event of the sample's files, analysis_v3/skims/src/skim_v3.cpp) gives the generator-level
fiducial selection of JHEP 11 (2017) 047 Table 4 (dressed leptons) and the stage-0 bin.  The
cross sections use the normalization of the selection scan (the same sigma_eff and
sum(genEventSumw) over the same files as every reconstructed yield):
    sigma(class) = sigma_eff x sum(genWeight over the class) / sum(genEventSumw).
Written per production mode and in total: sigma_eff (the 4l cross section of the sample),
sigma_fid (inclusive, per generator final state, per differential bin of pT(H), N(jets) and the
leading-jet pT; h4l_fiducial.DIFFERENTIAL), the fiducial acceptance A_fid = sigma_fid / sigma_eff,
and the stage-0 cross sections per class (ggH, VBF, VH_had, VH_lep) with |y_H| < 2.5 and above.
The GenTable weight sum is compared with sum(genEventSumw) of the scan: the delivered MC carries the
production preselection (two muons or two electrons, any HLT bit), so the table misses the generated
events that failed it (about 10 %, essentially all outside the fiducial volume; documented
limitation: the fiducial events among them are not counted in sigma_fid).  Writes production_v3/signal_model/<select>/<label>/fiducial_xsec.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import uproot

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_fiducial as fid  # noqa: E402
import h4l_select_io as io  # noqa: E402

SIGNAL = {"GluGluToHToZZ_M125": "ggH", "VBF_HToZZ_M125": "VBF", "VHToZZ_M125": "VH"}
GEN_BRANCHES = ["weight", "h_y", "stage0", "vh_class", "fid_pass", "fid_final_state", "fid_pt4l", "fid_njets", "fid_jet1_pt"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--skims", default="v2")
    args = parser.parse_args()
    out_dir = PRODUCTION / "signal_model" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    scan = io.load_scan(args.select)
    plan = json.loads((PRODUCTION / "skims" / args.skims / "plan.json").read_text(encoding="utf-8"))
    report = {"schema": "h4l_v3_fiducial_xsec/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
              "skim_plan": str(PRODUCTION / "skims" / args.skims / "plan.json"), "differential_bins": {k: [e if np.isfinite(e) else None for e in v]
                                                                                                        for k, v in fid.DIFFERENTIAL.items()},
              "stxs_abs_y": fid.STXS_ABS_Y, "units": "fb", "modes": {}}
    totals = {}
    for sample, mode in SIGNAL.items():
        info = scan["samples"][sample]
        tasks = [t for t in plan["tasks"] if t["config"]["sample"] == sample]
        rows = {b: [] for b in GEN_BRANCHES}
        files = []
        for task in tasks:
            sidecar = json.loads(Path(task["outputs"]["json"]).read_text(encoding="utf-8"))
            files += [f["path"] if isinstance(f, dict) and "path" in f else f for f in sidecar.get("files", [])]
            with uproot.open(task["outputs"]["root"]) as f:
                arrays = f["GenTable"].arrays(GEN_BRANCHES, library="np")
            for b in GEN_BRANCHES:
                rows[b].append(arrays[b])
        rows = {b: np.concatenate(v) for b, v in rows.items()}
        w = rows["weight"].astype(np.float64)
        sumw_table = float(np.sum(w))
        sigma_eff_fb = info["sigma_eff_pb"] * 1000.0
        scale = sigma_eff_fb / info["genEventSumw"]
        entry = {"sigma_eff_fb": sigma_eff_fb, "genEventSumw_scan": info["genEventSumw"], "genWeight_sum_gentable": sumw_table,
                 "gentable_over_scan": sumw_table / info["genEventSumw"], "gentable_rows": int(len(w))}
        in_fid = rows["fid_pass"].astype(bool)
        entry["sigma_fid"] = float(np.sum(w[in_fid]) * scale)
        # The delivered MC carries the production preselection (two muons or two electrons, any HLT bit): the
        # GenTable misses the generated events that failed it (gentable_over_scan < 1), essentially all outside
        # the fiducial volume (four fiducial leptons are almost never left without two reconstructed leptons of
        # one flavour).  A_fid is relative to every generated event (sigma_eff).
        entry["A_fid"] = entry["sigma_fid"] / sigma_eff_fb
        entry["A_fid_of_preselected"] = entry["sigma_fid"] / (sumw_table * scale)
        entry["sigma_fid_final_state"] = {name: float(np.sum(w[in_fid & (rows["fid_final_state"] == fs)]) * scale)
                                          for fs, name in fid.FINAL_STATES.items()}
        entry["sigma_fid_differential"] = {}
        for obs, edges in fid.DIFFERENTIAL.items():
            idx = fid.bin_of(fid.gen_observable(obs, rows), edges)
            entry["sigma_fid_differential"][obs] = [float(np.sum(w[in_fid & (idx == k)]) * scale) for k in range(len(edges) - 1)]
        classes = fid.production_class(mode, rows["vh_class"])
        central = np.abs(rows["h_y"]) < fid.STXS_ABS_Y
        entry["stage0"] = {}
        for c in sorted(set(classes)):
            sel = classes == c
            entry["stage0"][c] = {"central": float(np.sum(w[sel & central]) * scale), "forward": float(np.sum(w[sel & ~central]) * scale)}
        entry["stage0_codes"] = {str(int(k)): float(np.sum(w[rows["stage0"] == k]) * scale) for k in np.unique(rows["stage0"])}
        report["modes"][mode] = entry
        for key in ("sigma_eff_fb", "sigma_fid"):
            totals[key] = totals.get(key, 0.0) + entry[key]
        print(f"[fid] {mode}: sigma_eff {sigma_eff_fb:.4f} fb, sigma_fid {entry['sigma_fid']:.4f} fb, A_fid {entry['A_fid']:.3f}; "
              f"GenTable/scan sum(genWeight) {entry['gentable_over_scan']:.6f}", flush=True)
    totals["sigma_fid_final_state"] = {name: sum(report["modes"][m]["sigma_fid_final_state"][name] for m in report["modes"])
                                       for name in fid.FINAL_STATES.values()}
    totals["sigma_fid_differential"] = {obs: [sum(report["modes"][m]["sigma_fid_differential"][obs][k] for m in report["modes"])
                                              for k in range(len(edges) - 1)] for obs, edges in fid.DIFFERENTIAL.items()}
    totals["A_fid"] = totals["sigma_fid"] / totals["sigma_eff_fb"]
    report["total"] = totals
    out_dir.mkdir(parents=True)
    (out_dir / "fiducial_xsec.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[fid] total sigma_fid {totals['sigma_fid']:.4f} fb (A_fid {totals['A_fid']:.3f}); per final state "
          f"{ {k: round(v, 4) for k, v in totals['sigma_fid_final_state'].items()} }\n[fid] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
