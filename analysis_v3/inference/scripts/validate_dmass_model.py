"""Validation of the per-event-width signal model of the 3D fit: m4l projections per D_mass bin (final review 2026-09-25).

    pixi run py -- analysis_v3/inference/scripts/validate_dmass_model.py --select v5 --model model_cat_r1 --label dmass_v1

For every final state and D_mass = dm/m4l template bin of the model: the ggH + VBF signal MC (m_H = 125 GeV, the
refitted mass, 105-140 GeV, weights with the scale factors) against the model's conditional density, the DCB of
width s x D_mass x 125 GeV averaged over the MC events of the bin (the density the likelihood uses for each event),
both normalized in the window; the chi2 per bin, the mean and the rms of the two are recorded.  Writes
production_v3/results/<select>/<label>/ (dmass.json, dmass_<fs>.png).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_select_io as io  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    out = PRODUCTION / "results" / args.select / args.label
    if out.exists():
        raise SystemExit(f"{out} exists")
    model = json.loads((PRODUCTION / "inference" / args.select / args.model / "model.json").read_text(encoding="utf-8"))
    signal_model = json.loads(Path(model["signal_model"]).read_text(encoding="utf-8"))
    lo, hi = model["window"]
    scan = io.load_scan(args.select)
    rows = io.read(scan, ["GluGluToHToZZ_M125", "VBF_HToZZ_M125"], "SR", ["m4l_refit", "m4l_refit_err", "final_state"],
                   cut=lambda a: (a["m4l_refit"] > lo) & (a["m4l_refit"] < hi))
    out.mkdir(parents=True)
    report = {"schema": "h4l_v3_dmass_validation/1", "model": args.model, "final_states": {}}
    edges = np.linspace(lo, hi, 71)
    centres = 0.5 * (edges[1:] + edges[:-1])
    grid = np.linspace(lo, hi, 701)
    for code, fs in FINAL_STATES.items():
        v = signal_model["final_states"][fs]["m4l_refit"]["dcb_per_event_width"]["values"]
        e_edges = np.array(model["final_states"][fs]["e_edges"]["m4l_refit"])
        sel = rows["final_state"] == code
        m, e, w = rows["m4l_refit"][sel], (rows["m4l_refit_err"] / rows["m4l_refit"])[sel], rows["w"][sel]
        fig, axes = plt.subplots(1, len(e_edges) - 1, figsize=(4 * (len(e_edges) - 1), 3.6))
        per_bin = []
        for k, ax in enumerate(axes):
            in_bin = (e >= e_edges[k]) & (e < e_edges[k + 1])
            if not np.any(in_bin):
                continue
            h, _ = np.histogram(m[in_bin], bins=edges, weights=w[in_bin])
            h2, _ = np.histogram(m[in_bin], bins=edges, weights=w[in_bin] ** 2)
            norm = float(np.sum(w[in_bin]))
            # The model density of the bin: the per-event conditional DCBs averaged over the bin's events (a subsample
            # of at most 2000 events for speed, weights kept).
            idx = np.arange(int(in_bin.sum()))
            if len(idx) > 2000:
                idx = np.random.default_rng(1).choice(len(idx), 2000, replace=False)
            eb, wb = e[in_bin][idx], w[in_bin][idx]
            dens = np.zeros(len(grid))
            for ei, wi in zip(eb, wb):
                dens += wi * shapes.dcb_pdf(grid, v["mean"], v["width"] * ei * 125.0, v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"], lo, hi)
            dens /= np.sum(wb)
            expected = np.interp(centres, grid, dens) * (edges[1] - edges[0]) * norm
            ok = h2 > 0
            chi2 = float(np.sum((h[ok] - expected[ok]) ** 2 / h2[ok]))
            mean_mc = float(np.sum(w[in_bin] * m[in_bin]) / norm)
            rms_mc = float(np.sqrt(np.sum(w[in_bin] * (m[in_bin] - mean_mc) ** 2) / norm))
            mean_model = float(np.trapezoid(grid * dens, grid))
            rms_model = float(np.sqrt(np.trapezoid((grid - mean_model) ** 2 * dens, grid)))
            per_bin.append({"bin": k, "range": [float(e_edges[k]), float(e_edges[k + 1])], "events": int(in_bin.sum()),
                            "chi2": chi2, "ndf": int(ok.sum()), "mean_mc": mean_mc, "mean_model": mean_model, "rms_mc": rms_mc,
                            "rms_model": rms_model})
            ax.errorbar(centres, h / norm, yerr=np.sqrt(h2) / norm, fmt="k.", ms=3, label="signal MC")
            ax.plot(centres, expected / norm, color="#e42536", label="model")
            ax.set_title(f"{fs}, D_mass [{e_edges[k]:.4f}, {min(e_edges[k + 1], 1):.4f})\nrms MC {rms_mc:.2f} / model {rms_model:.2f} GeV",
                         fontsize=8)
            ax.set_xlabel("m4l (refit) [GeV]")
            ax.set_yscale("log")
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out / f"dmass_{fs}.png", dpi=110)
        plt.close(fig)
        report["final_states"][fs] = per_bin
        print(f"[dmass] {fs}: " + "; ".join(f"bin {b['bin']} rms {b['rms_mc']:.3f}/{b['rms_model']:.3f} chi2/ndf {b['chi2']:.0f}/{b['ndf']}"
                                            for b in per_bin), flush=True)
    (out / "dmass.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(f"[dmass] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
