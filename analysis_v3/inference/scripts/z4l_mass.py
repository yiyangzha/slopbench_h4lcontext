"""The Z -> 4l mass cross-check of the lepton momentum scale (JHEP 11 (2017) 047 10.3).

    pixi run py -- analysis_v3/inference/scripts/z4l_mass.py --select v5 --label z4l_v1

Z -> 4l candidates (the SRZ4l region of h4l_select: m_Z2 > 4 GeV) with 80 < m4l < 100 GeV (no Z1
mass constraint).  Per final state a double-sided Crystal Ball is fitted (weighted unbinned
likelihood) to the ZZTo4L + gg -> ZZ MC (shape only: the scale factors do not enter); the data are
fitted with the MC shape (width and tails fixed) and a free mean, so that
    m_Z(data) = m_Z(PDG) + mean(data) - mean(MC),
the MC carrying the generated m_Z of 91.1876 GeV.  Statistical uncertainty from the profile
likelihood (MINOS); systematic: the lepton momentum-scale uncertainty (systematics
configuration) times m_Z per flavour fraction, and the MC-mean statistical uncertainty.  A common
shift fitted to the three final states gives the combined value; the compatibility of the three
channels and of the combination with the PDG value are chi-square p-values.  Writes
production_v3/inference/<select>/<label>/z4l_mass.json and plots/.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from iminuit import Minuit  # noqa: E402
from scipy import stats  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_select_io as io  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

Z_MASS = 91.1876
WINDOW = (80.0, 100.0)
FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}
FLAVOUR_FRACTION = {"4mu": (1.0, 0.0), "4e": (0.0, 1.0), "2e2mu": (0.5, 0.5)}
MC = ["ZZTo4L", "GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"]


def dcb_nll(m, w, q):
    mean, width, al, nl, ar, nr = q
    dens = shapes.dcb_pdf(m, mean, width, al, nl, ar, nr, *WINDOW)
    return -float(np.sum(w * np.log(np.maximum(dens, 1e-300))))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--systematics", type=Path, default=REPO / "analysis_v3/inference/config/systematics_ul16_v3.json")
    args = parser.parse_args()
    out_dir = PRODUCTION / "inference" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    syst = json.loads(args.systematics.read_text(encoding="utf-8"))
    d_scale = syst["scale_uncertainty"]
    scan = io.load_scan(args.select)
    inside = lambda a: (a["m4l"] > WINDOW[0]) & (a["m4l"] < WINDOW[1])  # noqa: E731
    data = io.read(scan, ["pseudo_data"], "SRZ4l", ["m4l", "final_state"], cut=inside)
    mc = io.read(scan, MC, "SRZ4l", ["m4l", "final_state"], cut=inside)
    results, shifts = {}, {}
    out_dir.mkdir(parents=True)
    (out_dir / "plots").mkdir()
    for fs, name in FINAL_STATES.items():
        mm, mw = mc["m4l"][mc["final_state"] == fs], mc["w"][mc["final_state"] == fs]
        norm = len(mm) / np.sum(mw)
        fit_mc = Minuit(lambda mean, width, al, nl, ar, nr: dcb_nll(mm, mw * norm, (mean, width, al, nl, ar, nr)),
                        mean=90.8, width=1.5, al=1.0, nl=3.0, ar=1.5, nr=5.0)
        fit_mc.errordef = Minuit.LIKELIHOOD
        fit_mc.limits = [(85, 95), (0.2, 6), (0.2, 10), (1.01, 80), (0.2, 10), (1.01, 80)]
        fit_mc.migrad()
        fit_mc.hesse()
        v = {k: float(fit_mc.values[k]) for k in ("mean", "width", "al", "nl", "ar", "nr")}
        # The MC-mean error of a weighted fit: HESSE of the normalized weights scaled by the effective statistics.
        n_eff = np.sum(mw) ** 2 / np.sum(mw ** 2)
        mc_mean_error = float(fit_mc.errors["mean"] * math.sqrt(len(mm) / n_eff))
        dm = data["m4l"][data["final_state"] == fs]
        fit_d = Minuit(lambda mean: dcb_nll(dm, np.ones(len(dm)), (mean, v["width"], v["al"], v["nl"], v["ar"], v["nr"])), mean=v["mean"])
        fit_d.errordef = Minuit.LIKELIHOOD
        fit_d.limits["mean"] = (85, 95)
        fit_d.migrad()
        fit_d.minos()
        me = fit_d.merrors["mean"]
        shift = float(fit_d.values["mean"]) - v["mean"]
        f_mu, f_e = FLAVOUR_FRACTION[name]
        # The muon and electron scale uncertainties are independent nuisances: their effects add in quadrature.
        syst_scale = Z_MASS * math.hypot(f_mu * d_scale["mu"], f_e * d_scale["e"])
        syst_total = math.hypot(syst_scale, mc_mean_error)
        results[name] = {"n_data": int(len(dm)), "n_mc_expected": float(np.sum(mw)), "mc_dcb": v, "mc_mean_error": mc_mean_error,
                         "data_mean": float(fit_d.values["mean"]), "stat": [float(me.lower), float(me.upper)],
                         "m_z": Z_MASS + shift, "syst_scale": syst_scale, "syst": syst_total}
        shifts[name] = (shift, 0.5 * (abs(me.lower) + abs(me.upper)), syst_total)
        fig, ax = plt.subplots(figsize=(6, 4.2))
        edges = np.arange(WINDOW[0], WINDOW[1] + 1e-9, 1.0)
        centres = 0.5 * (edges[1:] + edges[:-1])
        h, _ = np.histogram(dm, bins=edges)
        ax.errorbar(centres, h, yerr=np.sqrt(np.maximum(h, 1)), fmt="ko", ms=3, label=f"data ({len(dm)})")
        grid = np.linspace(*WINDOW, 400)
        dens = shapes.dcb_pdf(grid, fit_d.values["mean"], v["width"], v["al"], v["nl"], v["ar"], v["nr"], *WINDOW)
        ax.plot(grid, dens * len(dm) * 1.0, color="#e42536", label=f"fit: m_Z = {Z_MASS + shift:.2f} {me.lower:+.2f}/{me.upper:+.2f} GeV")
        mc_dens = shapes.dcb_pdf(grid, v["mean"], v["width"], v["al"], v["nl"], v["ar"], v["nr"], *WINDOW)
        ax.plot(grid, mc_dens * np.sum(mw), color="#5790fc", ls="--", label=f"MC ({np.sum(mw):.1f} expected)")
        ax.set_xlabel("m4l [GeV]")
        ax.set_ylabel("events / GeV")
        ax.set_title(f"Z -> 4l, {name} (no Z1 constraint)", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out_dir / "plots" / f"z4l_{name}.png", dpi=120)
        fig.savefig(out_dir / "plots" / f"z4l_{name}.pdf")
        plt.close(fig)
        print(f"[z4l] {name}: m_Z = {Z_MASS + shift:.3f} {me.lower:+.3f}/{me.upper:+.3f} (stat) +- {syst_total:.3f} (syst) GeV "
              f"({len(dm)} events, MC {np.sum(mw):.1f})", flush=True)
    # Combination: inverse-variance mean of the shifts (statistical), compatibility chi2.
    s = np.array([v[0] for v in shifts.values()])
    e = np.array([v[1] for v in shifts.values()])
    sy = np.array([v[2] for v in shifts.values()])
    w = 1 / e ** 2
    mean_shift = float(np.sum(w * s) / np.sum(w))
    chi2_mutual = float(np.sum(((s - mean_shift) / e) ** 2))
    comb_stat = float(1 / math.sqrt(np.sum(w)))
    comb_syst = float(np.sum(w * sy) / np.sum(w))
    combined = {"m_z": Z_MASS + mean_shift, "stat": comb_stat, "syst": comb_syst,
                "mutual_compatibility_p": float(stats.chi2.sf(chi2_mutual, len(s) - 1)),
                "pdg_compatibility_p": float(stats.chi2.sf((mean_shift / math.hypot(comb_stat, comb_syst)) ** 2, 1))}
    report = {"schema": "h4l_v3_z4l_mass/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
              "window": list(WINDOW), "method": "DCB of the ZZ MC per final state, data fitted with a free mean; m_Z = 91.1876 + shift",
              "final_states": results, "combined": combined, "systematics_config": str(args.systematics)}
    (out_dir / "z4l_mass.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[z4l] combined m_Z = {combined['m_z']:.3f} +- {comb_stat:.3f} (stat) +- {comb_syst:.3f} (syst) GeV; mutual p "
          f"{combined['mutual_compatibility_p']:.2f}, PDG p {combined['pdg_compatibility_p']:.2f}\n[z4l] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
