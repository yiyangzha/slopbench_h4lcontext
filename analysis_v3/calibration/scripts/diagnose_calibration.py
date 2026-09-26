"""Charge, signed-eta and phi diagnostics of a calibration payload (stage 3a).

    pixi run py -- analysis_v3/calibration/scripts/diagnose_calibration.py --version v2 --run nominal \\
        --payload-iteration 2 --label diag_iter03

The per-lepton calibration bins are in (pT, |eta|).  A momentum-scale effect
that depends on the charge, on the sign of eta or on phi (the Rochester-style
parameters of the AN muon corrections) averages out in those bins but leaves
pair-level inconsistencies.  This script
 1. re-runs the histogram jobs of the data and the MC roles with the payload of
    <run>/iter_<n>/payload.json and FullPairs trees (frozen zpeak_histograms);
 2. for each flavour and each leg charge, in slices of that leg's signed eta
    and of its phi (the other leg anywhere), fits the data mass distribution
    with the MC one scaled by k (event-level: every MC pair rescaled, Poisson
    likelihood scanned in ln k, parabola at the minimum) and reports
    ln k(data/MC) per slice, and the same for the relative width;
 3. plots them into <run>/<label>/plots/ and writes <run>/<label>/diagnostics.json.
A charge-antisymmetric curvature bias appears as opposite modulations for the
two charges; a signed-eta or phi dependence as a modulation common to both.
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
import uproot  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402
import run_calibration as rc  # noqa: E402

FLAVOURS = {"mm": ("muon", 13, r"\mu"), "ee": ("electron", 11, "e")}


def load_pairs(files: list, pdg: int) -> dict:
    columns = {"m": [], "w": [], "pt": [], "eta": [], "phi": [], "charge": []}
    for path in files:
        with uproot.open(path) as f:
            tree = f["FullPairs"]
            arrays = tree.arrays(["flavour", "m", "w", "pt", "eta", "phi", "charge"], library="np")
        keep = arrays["flavour"] == pdg
        for key in columns:
            columns[key].append(arrays[key][keep])
    return {key: np.concatenate(value) for key, value in columns.items()}


def shift_fit(data_m: np.ndarray, mc_m: np.ndarray, mc_w: np.ndarray, window=(82.0, 100.0), width=0.1) -> tuple:
    """ln k(data/MC) by an event-level Poisson likelihood scan; also the
    relative-width ratio from the interquartile ranges."""
    edges = np.arange(window[0], window[1] + 1e-9, width)
    d, _ = np.histogram(data_m, edges)
    grid = np.linspace(-0.006, 0.006, 61)
    nll = []
    for lnk in grid:
        t, _ = np.histogram(mc_m * math.exp(lnk), edges, weights=mc_w)
        t = np.maximum(t, 1e-9)
        t = t / t.sum() * d.sum()
        nll.append(float(np.sum(t - d * np.log(t))))
    nll = np.array(nll)
    best = int(np.argmin(nll))
    lo, hi = max(best - 4, 0), min(best + 5, len(grid))
    if hi - lo < 3:
        return math.nan, math.nan, math.nan
    a, b, _ = np.polyfit(grid[lo:hi], nll[lo:hi], 2)
    if a <= 0:
        return math.nan, math.nan, math.nan
    lnk = -b / (2 * a)
    error = 1 / math.sqrt(2 * a)
    def iqr(values, weights=None):
        order = np.argsort(values)
        v = values[order]
        cw = np.cumsum(weights[order] if weights is not None else np.ones_like(v))
        cw /= cw[-1]
        return np.interp(0.75, cw, v) - np.interp(0.25, cw, v)
    sel_d = (data_m > window[0]) & (data_m < window[1])
    sel_m = (mc_m > window[0]) & (mc_m < window[1])
    width_ratio = iqr(data_m[sel_d]) / iqr(mc_m[sel_m], mc_w[sel_m]) if sel_d.sum() > 50 and sel_m.sum() > 50 else math.nan
    return float(lnk), float(error), float(width_ratio)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--extract", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--payload-iteration", type=int, help="use <run>/iter_<n>/payload.json")
    parser.add_argument("--final", action="store_true", help="use the final <run>/payload.json")
    parser.add_argument("--label", required=True)
    parser.add_argument("--manifests", default="v2")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    base = rc.PRODUCTION / "calibration" / args.version / args.run
    out = base / args.label
    config = json.loads((base / "config.json").read_text(encoding="utf-8"))
    if args.final == (args.payload_iteration is not None):
        raise SystemExit("ERROR: give exactly one of --payload-iteration and --final")
    payload_path = base / "payload.json" if args.final else base / f"iter_{args.payload_iteration:02d}" / "payload.json"
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    if args.final:
        payload = payload["flavours"]
    scan, per_file, _ = rc.load_inputs(args.extract, args.manifests)
    frozen = rc.freeze(args.version)
    ctx = {"config": config, "scan": scan, "frozen": frozen, "workers": args.workers}
    mc_samples = [n for n, s in scan["samples"].items() if s["kind"] == "mc"]
    specs = {"data": {"role": "data", "samples": {"pseudo_data": 1.0}, "write_full_pairs": True},
             "mc": {"role": "mc", "samples": {s: rc.mc_scale(scan, per_file, s, -1) for s in mc_samples}, "write_full_pairs": True}}
    files = rc.histogram_jobs(ctx, out, payload, "hist", specs)
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    lumi = scan["samples"]["pseudo_data"]["lumi_fb"]
    results = {"payload": str(payload_path), "slices": {}}
    for tag, (name, pdg, symbol) in FLAVOURS.items():
        data = load_pairs(files["data"], pdg)
        mc = load_pairs(files["mc"], pdg)
        max_eta = 2.4 if tag == "mm" else 2.5
        variables = {"eta": np.linspace(-max_eta, max_eta, 13), "phi": np.linspace(-math.pi, math.pi, 13)}
        results["slices"][tag] = {}
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        for column, (variable, edges) in enumerate(variables.items()):
            for charge, colour in ((+1, "#e42536"), (-1, "#5790fc")):
                values, errors, widths = [], [], []
                for lo, hi in zip(edges[:-1], edges[1:]):
                    def select(sample):
                        mask = np.zeros(len(sample["m"]), dtype=bool)
                        for leg in (0, 1):
                            mask |= (sample["charge"][:, leg] == charge) & (sample[variable][:, leg] >= lo) & \
                                    (sample[variable][:, leg] < hi)
                        return mask
                    sd, sm = select(data), select(mc)
                    lnk, err, wr = shift_fit(data["m"][sd], mc["m"][sm], mc["w"][sm])
                    values.append(lnk)
                    errors.append(err)
                    widths.append(wr)
                centres = 0.5 * (edges[1:] + edges[:-1])
                sign = "+" if charge > 0 else "-"
                axes[0, column].errorbar(centres, 1e3 * np.array(values), 1e3 * np.array(errors), fmt="o-", color=colour,
                                         label=f"${symbol}^{sign}$ in slice")
                axes[1, column].plot(centres, widths, "o-", color=colour, label=f"${symbol}^{sign}$ in slice")
                results["slices"][tag][f"{variable}_{sign}"] = {"edges": edges.tolist(), "lnk": values, "error": errors,
                                                               "iqr_ratio": widths}
            axes[0, column].axhline(0, color="grey", lw=0.5)
            axes[0, column].set_ylabel(r"$\ln k$ (data/MC) $\times 10^3$")
            axes[1, column].axhline(1, color="grey", lw=0.5)
            axes[1, column].set_ylabel("IQR ratio data/MC")
            for row in (0, 1):
                axes[row, column].set_xlabel(f"{variable} of the {symbol} in the slice".replace("\\mu", "mu"))
                axes[row, column].legend(fontsize="small")
        fig.suptitle(f"{name} pairs, payload {'final' if args.final else 'of iteration ' + str(args.payload_iteration)}")
        h4l_style.save(fig, plots / f"charge_eta_phi_{tag}")
    (out / "diagnostics.json").write_text(json.dumps(results, indent=1) + "\n", encoding="utf-8")
    print(f"[diagnostics] {out / 'diagnostics.json'}; plots {plots}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
