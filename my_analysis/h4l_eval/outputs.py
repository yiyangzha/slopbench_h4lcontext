"""RESULT.json (EVAL_CONTRACT v0.4 with the per-flavour calibration block), MODEL.json (the fitted binned likelihood as a
pyhf workspace at the fitted m_H) and the diagnostics (results_full.json, summary.txt)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import config as C
from .model import NUISANCES

DISCOVERY_Z = 5.0


def sym(iv):
    return 0.5 * (abs(iv[0]) + abs(iv[1]))


def workspace(lik, mh: float) -> dict:
    """The fitted model as a pyhf workspace: at m_H = mh the signal templates are fixed; the lepton scale and resolution
    act as histosys, the efficiencies and Z+X as normsys, mu is the POI (normfactor)."""
    sig = lik.signal_templates(mh)
    channels, observations = [], []
    for ch, (nom, var), obs in zip(lik.channels, sig, lik.observed):
        mods = [{"name": "mu", "type": "normfactor", "data": None}]
        for name, v in var.items():
            if v is None:
                continue
            up, down = v
            mods.append({"name": name, "type": "histosys", "data": {"hi_data": up.tolist(), "lo_data": down.tolist()}})
        for fl in ("mu", "e"):
            k = ch["signal"]["kappa_eff"][fl]
            if k > 0:
                mods.append({"name": f"eff_{fl}", "type": "normsys", "data": {"hi": 1.0 + k, "lo": 1.0 - k}})
        samples = [{"name": "signal", "data": nom.tolist(), "modifiers": mods}]
        for role in ("qqZZ", "ggZZ"):
            b = ch[role]
            bm = [{"name": f"eff_{fl}", "type": "normsys", "data": {"hi": 1.0 + b["kappa_eff"][fl], "lo": 1.0 - b["kappa_eff"][fl]}}
                  for fl in ("mu", "e") if b["kappa_eff"][fl] > 0]
            samples.append({"name": role, "data": b["template"].tolist(), "modifiers": bm})
        z = ch["zx"]
        samples.append({"name": "ZX", "data": z["template"].tolist(),
                        "modifiers": [{"name": z["nuisance"], "type": "normsys", "data": {"hi": 1.0 + z["kappa"], "lo": 1.0 - z["kappa"]}}]})
        channels.append({"name": ch["name"], "samples": samples})
        observations.append({"name": ch["name"], "data": [float(x) for x in obs]})
    return {"channels": channels, "observations": observations,
            "measurements": [{"name": "h4l", "config": {"poi": "mu", "parameters": [{"name": "mu", "bounds": [[lik.mu_low, 50.0]], "inits": [1.0]}]}}],
            "version": "1.0.0"}


def write_outputs(out_dir: Path, ds, lik, inf, cal, cal_summary, tnp, sel_eff, zx_est, rates, meta, norm, timings, log):
    best = inf["best"]
    pois = []
    for name in ("mH", "mu"):
        total = sym(inf["total"][name])
        stat = sym(inf["stat"][name])
        pois.append({"name": name, "value": best[name], "stat": stat, "syst": inf["syst"][name], "total": total})
    calibration = {
        "scale_shift": {fl: {"value": cal_summary[fl]["scale_shift"]["value"], "unc": cal_summary[fl]["scale_shift"]["unc"]}
                        for fl in ("muon", "electron")},
        "smear": {fl: {"value": cal_summary[fl]["smear"]["value"], "unc": cal_summary[fl]["smear"]["unc"]} for fl in ("muon", "electron")},
        "sel_eff": {fl: {"value": sel_eff[fl]["value"], "unc": sel_eff[fl]["unc"]} for fl in ("muon", "electron")},
    }
    quality = {"gof_pvalue": {"value": inf["gof"]["p_value"],
                              "note": f"toys: saturated Poisson LR of the fitted binned model ({inf['gof']['n_bins']} bins, constraint "
                                      f"terms included), {inf['gof']['n_valid']} toys from the best fit with the global observables "
                                      f"drawn, each refitted with m_H floating"},
               "coverage": {"value": inf["coverage"]["value"],
                            "note": f"toys: fraction of {inf['coverage']['n']} pseudo-experiments from the best fit whose MINOS 68 % "
                                    f"interval of mu covers the generating mu"}}
    result = {"pois": pois, "significance_obs": inf["z_obs"], "discovery": bool(inf["z_obs"] >= DISCOVERY_Z),
              "significance_exp": inf["z_exp"], "calibration": calibration, "quality": quality,
              "conventions": {
                  "calibration": "per flavour (muon, electron) under the registered names; each value is the average of the "
                                 "per-bin measurements over the data leptons of the Z -> ll calibration sample; scale_shift = "
                                 "data/MC - 1 of the lepton momentum, smear = extra relative per-lepton pT resolution in "
                                 "quadrature, sel_eff = tag-and-probe data/MC efficiency ratio of the full single-lepton selection",
                  "uncertainties": "total and stat: symmetrized MINOS intervals (stat with the nuisances fixed at their best fit); "
                                   "syst: quadrature sum of the fixed-nuisance impacts",
                  "discovery": f"observed local significance >= {DISCOVERY_Z}"}}
    (out_dir / "RESULT.json").write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    ws = workspace(lik, best["mH"])
    (out_dir / "MODEL.json").write_text(json.dumps(ws) + "\n", encoding="utf-8")
    # The description of the workspace in a side file (the pyhf schema admits no extra top-level key).
    (out_dir / "MODEL_meta.json").write_text(json.dumps({
        "description": "H -> ZZ* -> 4l binned likelihood (m4l Z1-refitted, 0.5 GeV bins in 105-140 GeV; channels final state x "
                       "D_mass bin) at the fitted m_H", "mH_fixed_GeV": best["mH"],
        "interpolation": "pyhf defaults: histosys code 0, normsys code 1 (the fitted likelihood uses the same)"}, indent=1) + "\n",
        encoding="utf-8")
    full = {"schema": "h4l_eval_results/1", "dataset": ds["root"], "lumi_fb": ds["lumi_fb"], "RESULT": result,
            "fit": {k: v for k, v in inf.items()}, "calibration": {"summary": cal_summary, "report": cal["report"]},
            "tnp": {str(code): t for code, t in tnp.items()}, "sel_eff": sel_eff,
            "zx": {"os": zx_est, "fake_rates": {str(k): {kk: np.asarray(vv).tolist() for kk, vv in v.items()} for k, v in rates.items()}},
            "model": meta, "normalization": norm, "timings_s": timings}
    (out_dir / "results_full.json").write_text(json.dumps(full, indent=1, default=_default) + "\n", encoding="utf-8")
    lines = [f"mu = {best['mu']:.4f} +{inf['total']['mu'][1]:.4f}/{inf['total']['mu'][0]:.4f} (stat {sym(inf['stat']['mu']):.4f}, "
             f"syst {inf['syst']['mu']:.4f})",
             f"mH = {best['mH']:.3f} +{inf['total']['mH'][1]:.3f}/{inf['total']['mH'][0]:.3f} GeV (stat {sym(inf['stat']['mH']):.3f}, "
             f"syst {inf['syst']['mH']:.3f})",
             f"Z_obs = {inf['z_obs']:.2f}, Z_exp = {inf['z_exp']:.2f}; GoF p = {inf['gof']['p_value']:.3f}; coverage {inf['coverage']['value']:.3f}"]
    for fl in ("muon", "electron"):
        lines.append(f"{fl}: scale_shift {calibration['scale_shift'][fl]['value']:+.5f} +- {calibration['scale_shift'][fl]['unc']:.5f}, "
                     f"smear {calibration['smear'][fl]['value']:.5f} +- {calibration['smear'][fl]['unc']:.5f}, "
                     f"sel_eff {calibration['sel_eff'][fl]['value']:.5f} +- {calibration['sel_eff'][fl]['unc']:.5f}")
    lines.append("timings [s]: " + ", ".join(f"{k} {v:.0f}" for k, v in timings.items()))
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        log(line)


def _default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
