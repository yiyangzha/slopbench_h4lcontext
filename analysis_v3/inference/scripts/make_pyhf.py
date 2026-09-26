"""MODEL.json: the binned pyhf (HistFactory) workspace of the H -> 4l likelihood (stage 7).

    pixi run py -- analysis_v3/inference/scripts/make_pyhf.py --model v5/model_v1 --fit fit_v1 --yr4 v1 --out MODEL.json

The unbinned likelihood (h4l_likelihood.py) is binned at the reported best-fit m_H: per final state
(categories merged), m4l (the Z1-refitted mass) in 2.5 GeV bins over the window x the D_bkg^kin
template bins.  Samples: signal (all production modes, the resonant DCB of the 2D model and the VH
non-resonant Landau, times the D templates; the POI mu is a normfactor), qqZZ, ggZZ, Z+X (their
shapes and D templates).  Modifiers (set (b) of the fit): the efficiency scale factors per flavour
(normsys from the flavour variations of the yields), the Z+X normalization per final state (normsys
with the asymmetric kappa), the lepton scale and resolution per flavour (histosys of the signal:
peak position and width moved by their uncertainties), the MC statistical uncertainty (staterror).
The observed counts are the binned data events.  The script fits the workspace with pyhf and
records mu_hat against the reported mu (the scoreboard requires |dmu| <= 0.25 sigma + 0.05).
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
import h4l_likelihood as lk  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

M_STEP = 2.5


def binned_shape(density, m_edges_fine, m_edges, d_table, m_template_edges):
    """Expected fraction per (m4l bin, D bin) of a density in m4l times P(D | m4l bin)."""
    out = np.zeros((len(m_edges) - 1, d_table.shape[1]))
    for i in range(len(m_edges) - 1):
        grid = np.linspace(m_edges[i], m_edges[i + 1], 26)
        centres = 0.5 * (grid[1:] + grid[:-1])
        mass = float(np.sum(density(centres)) * (grid[1] - grid[0]))
        ti = min(max(int(np.searchsorted(m_template_edges, 0.5 * (m_edges[i] + m_edges[i + 1]), side="right") - 1), 0), d_table.shape[0] - 1)
        out[i] = mass * d_table[ti]
    return out.ravel()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--fit", required=True, help="fit label whose best fit fixes m_H and the reported mu")
    parser.add_argument("--yr4", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    model_dir = PRODUCTION / "inference" / args.model
    model, events, signal_model, yr4 = fm.load(model_dir, args.yr4)
    fit = json.loads((model_dir / args.fit / "fit.json").read_text(encoding="utf-8"))
    syst = fit["systematics"]
    lik = lk.Model(model, signal_model, yr4, "2D", fit["refit"], "b", syst)
    best = fit["fit"]["values"]
    mh = best["mH"]
    p = np.zeros(len(lik.parameters))
    p[0], p[1] = 1.0, mh
    lo, hi = lik.lo, lik.hi
    m_edges = np.arange(lo, hi + 1e-9, M_STEP)
    channels_ws, observations = [], []
    d_scale = syst.get("scale_uncertainty", {"mu": 0.0005, "e": 0.0015})
    d_res = syst.get("resolution_uncertainty", {"mu": 0.10, "e": 0.10})
    for fs in model["final_states"]:
        fse = model["final_states"][fs]
        chans = [c for c in lik.channels if model["channels"][c]["final_state"] == fs]
        yields = {k: 0.0 for k in ("signal", "qqZZ", "ggZZ", "zx")}
        variations = {k: {"eff_mu": [0.0, 0.0], "eff_e": [0.0, 0.0]} for k in ("signal", "qqZZ", "ggZZ")}
        mcstat = {k: 0.0 for k in ("signal", "qqZZ", "ggZZ")}
        for c in chans:
            y = lik.yields(p, c)
            ys = sum(y[m] for m in lk.MODES)
            yields["signal"] += ys
            for name in ("qqZZ", "ggZZ", "zx"):
                yields[name] += y[name]
            ch = model["channels"][c]
            for name, blocks in (("signal", [ch["signal"][m] for m in lk.MODES]), ("qqZZ", [ch["backgrounds"]["qqZZ"]]),
                                 ("ggZZ", [ch["backgrounds"]["ggZZ"]])):
                for b in blocks:
                    scale_y = (ys / sum(ch["signal"][m]["yield"] for m in lk.MODES)) if name == "signal" and sum(ch["signal"][m]["yield"] for m in lk.MODES) > 0 else 1.0
                    for key, label in (("sf_muon", "eff_mu"), ("sf_electron", "eff_e")):
                        variations[name][label][0] += b["yield"] * scale_y * b[key][0]
                        variations[name][label][1] += b["yield"] * scale_y * b[key][1]
                    mcstat[name] += (b["mc_stat"] * scale_y) ** 2
        f_mu, f_e = lik.flavour_fraction[fs]
        g = dict(zip(("mu", "e"), lik.resolution_fraction[fs]))
        sm = fse["signal"][lik.variable]
        v = sm["dcb"]["values"]
        vh = sm["vh_nonresonant"]
        d_table_sig = np.array(fse["d_templates"]["signal"])
        m_template_edges = np.array(model["m_edges"])
        # Fraction of VH among the signal yield of the final state (non-resonant part).
        vh_share = sum(lik.yields(p, c)["VH"] for c in chans) / yields["signal"] if yields["signal"] > 0 else 0.0

        def signal_shape(mean_scale=1.0, width_scale=1.0):
            k = mh / 125.0

            def dens(m):
                res = shapes.dcb_pdf(m, v["mean"] * k * mean_scale, v["width"] * k * width_scale, v["alpha_l"], v["n_l"], v["alpha_r"],
                                     v["n_r"], lo, hi)
                nonres = shapes.landau_pdf(m, vh["landau_mpv"], vh["landau_width"], lo, hi)
                return (1 - vh_share * (1 - vh["f_res"])) * res + vh_share * (1 - vh["f_res"]) * nonres
            return binned_shape(dens, None, m_edges, d_table_sig, m_template_edges)

        samples = []
        nominal_sig = signal_shape() * yields["signal"]
        sig_mods = [{"name": "mu", "type": "normfactor", "data": None}]
        for label in ("eff_mu", "eff_e"):
            up, down = variations["signal"][label]
            if yields["signal"] > 0 and (abs(up) > 0 or abs(down) > 0):
                sig_mods.append({"name": label, "type": "normsys", "data": {"hi": 1 + up / yields["signal"], "lo": 1 + down / yields["signal"]}})
        for flav, frac in (("mu", f_mu), ("e", f_e)):
            if frac <= 0:
                continue
            up = signal_shape(mean_scale=1 + frac * d_scale[flav]) * yields["signal"]
            dn = signal_shape(mean_scale=1 - frac * d_scale[flav]) * yields["signal"]
            sig_mods.append({"name": f"scale_{flav}", "type": "histosys", "data": {"hi_data": up.tolist(), "lo_data": dn.tolist()}})
            up = signal_shape(width_scale=1 + g[flav] * d_res[flav]) * yields["signal"]
            dn = signal_shape(width_scale=max(1 - g[flav] * d_res[flav], 0.2)) * yields["signal"]
            sig_mods.append({"name": f"res_{flav}", "type": "histosys", "data": {"hi_data": up.tolist(), "lo_data": dn.tolist()}})
        samples.append({"name": "signal", "data": nominal_sig.tolist(), "modifiers": sig_mods})
        for name in ("qqZZ", "ggZZ"):
            coeff = fse["backgrounds"][name][f"bernstein_{lik.variable}"]
            shape = binned_shape(lambda m, c=coeff: shapes.bernstein_pdf(m, c, lo, hi), None, m_edges, np.array(fse["d_templates"][name]), m_template_edges)
            data = shape * yields[name]
            mods = []
            for label in ("eff_mu", "eff_e"):
                up, down = variations[name][label]
                if yields[name] > 0 and (abs(up) > 0 or abs(down) > 0):
                    mods.append({"name": label, "type": "normsys", "data": {"hi": 1 + up / yields[name], "lo": 1 + down / yields[name]}})
            rel = np.sqrt(mcstat[name]) / yields[name] if yields[name] > 0 else 0.0
            mods.append({"name": f"mcstat_{fs}", "type": "staterror", "data": (data * rel).tolist()})
            samples.append({"name": name, "data": data.tolist(), "modifiers": mods})
        zx = fse["zx"]
        shape = binned_shape(lambda m: shapes.zx_pdf(m, zx["shape"], lo, hi), None, m_edges,
                             np.array(fse["d_templates"]["zx"]), m_template_edges)
        samples.append({"name": "zx", "data": (shape * yields["zx"]).tolist(),
                        "modifiers": [{"name": f"zx_{fs}", "type": "normsys", "data": {"hi": 1 + zx["kappa_high"], "lo": 1 - zx["kappa_low"]}}]})
        channels_ws.append({"name": fs, "samples": samples})
        # Observed counts.
        counts = np.zeros((len(m_edges) - 1, len(model["d_edges"]) - 1))
        for c in chans:
            m = events[f"{c}__{lik.variable}"]
            d = events[f"{c}__d_bkg_kin"]
            h, _, _ = np.histogram2d(m, d, bins=[m_edges, np.array(model["d_edges"])])
            counts += h
        observations.append({"name": fs, "data": counts.ravel().tolist()})
    workspace = {"channels": channels_ws, "observations": observations,
                 "measurements": [{"name": "h4l_mu", "config": {"poi": "mu", "parameters": [{"name": "mu", "bounds": [[-20.0, 50.0]], "inits": [1.0]}]}}],
                 "version": "1.0.0"}
    import pyhf
    pyhf.set_backend("numpy", pyhf.optimize.minuit_optimizer(verbose=0))
    ws = pyhf.Workspace(workspace)
    pdf = ws.model()
    data = ws.data(pdf)
    bestfit = pyhf.infer.mle.fit(data, pdf)
    mu_hat = float(bestfit[pdf.config.poi_index])
    reported = fit["fit"]["pois"]["mu"]
    sigma = 0.5 * (abs(reported["total"][0]) + abs(reported["total"][1]))
    check = {"mu_pyhf": mu_hat, "mu_reported": reported["value"], "tolerance": 0.25 * sigma + 0.05,
             "consistent": abs(mu_hat - reported["value"]) <= 0.25 * sigma + 0.05, "mH_fixed": mh,
             "binning": f"m4l {M_STEP} GeV x {len(model['d_edges']) - 1} D bins per final state"}
    workspace["h4l_consistency"] = check
    out = args.out if args.out.is_absolute() else REPO / args.out
    out.write_text(json.dumps(workspace, indent=1) + "\n", encoding="utf-8")
    print(f"[pyhf] mu_hat {mu_hat:.4f} vs reported {reported['value']:.4f} (tolerance {check['tolerance']:.4f}): "
          f"{'consistent' if check['consistent'] else 'INCONSISTENT'}\n[pyhf] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
