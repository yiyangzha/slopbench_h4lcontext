"""Closure tests of the Z+X method (stage 5 validation).

    pixi run py -- analysis_v3/backgrounds/scripts/zx_closure.py --select v4 --label closure_v1

Two tests with the functions of zx_estimate.py:
  1. MC closure: fake rates measured in the Z + 1 loose lepton rows of the DY and TTbar MC (no prompt
     subtraction: these samples have no prompt third lepton) are applied with the OS method (AN Eq. 19,
     no ZZ term) and with the SS method (AN Eq. 20, (OS/SS)_MC of the same samples, SS-space fake rates with
     the electron conversion correction as in the data) to the DY + TTbar
     control rows; the predictions are compared per final state with the DY + TTbar signal-region
     yield (m4l > 70 GeV and 105-140 GeV), all at the data luminosity.
  2. 3P1F from 2P2F: sum_2P2F [f3/(1-f3) + f4/(1-f4)] against the observed 3P1F count (for the data,
     minus the prompt ZZ 3P1F rows of the ZZTo4L and gg -> ZZ MC), in the data and in the MC.
Uncertainties are statistical: sqrt(sum w^2) of the MC rows and Poisson bootstrap replicas of all rows
for the predictions (fake rates and control rows together).
Writes production_v3/backgrounds/<select>/<label>/closure.json and plots/.
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
sys.path.insert(0, str(REPO / "analysis_v3/backgrounds/scripts"))
import h4l_select_io as io  # noqa: E402
import zx_estimate as zx  # noqa: E402

WINDOWS = {"m4l_gt_70": (70.0, np.inf), "m4l_105_140": (105.0, 140.0)}


empty_like, predict_3p1f, os_prediction, ss_prediction = zx.empty_like, zx.predict_3p1f, zx.os_prediction, zx.ss_prediction


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    out_dir = PRODUCTION / "backgrounds" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    scan = io.load_scan(args.select)
    zx.configure(json.loads(Path(scan["plan"]).read_text(encoding="utf-8"))["tasks"][0]["selection_config"])
    sel = lambda d, m: {k: v[m] for k, v in d.items()}  # noqa: E731

    # --- 1. MC closure -------------------------------------------------------------------------------
    zl_mc = io.read(scan, zx.FAKE_MC, "ZL", zx.ZL_BRANCHES)
    cr_mc_all = io.read(scan, zx.FAKE_MC, "CR", zx.CR_BRANCHES)
    sr_mc = io.read(scan, zx.FAKE_MC, "SR", ["m4l", "final_state"])
    cr_mc = sel(cr_mc_all, (cr_mc_all["region"] == 0) & (cr_mc_all["cr_type"] <= 1))
    ss_mc = sel(cr_mc_all, (cr_mc_all["region"] == 0) & (cr_mc_all["cr_type"] == 2))
    none = empty_like(zl_mc)
    rates_os = zx.fake_rates(zl_mc, none, zx.SS_WINDOWS["os_window"])
    # SS-method fake rates as in the data: SS phase space, electrons corrected for the conversion content of the
    # same-sign control rows with m4l > 100 GeV (their loose legs weighted with the MC weights).
    ss_high = sel(ss_mc, ss_mc["m4l"] > 100.0)
    rates_ss, ss_info = zx.ss_rates(zl_mc, none, dict(ss_high, boot=ss_high["w"]))
    os_ss = {name: r["ratio"] for name, r in zx.os_ss_ratio(cr_mc_all).items()}
    truth, truth_err = {}, {}
    for key, (lo, hi) in WINDOWS.items():
        inside = (sr_mc["m4l"] > lo) & (sr_mc["m4l"] < hi)
        truth[key] = {name: float(np.sum(sr_mc["w"][inside & (sr_mc["final_state"] == fs)])) for fs, name in zx.FINAL_STATES.items()}
        truth_err[key] = {name: float(np.sqrt(np.sum(sr_mc["w"][inside & (sr_mc["final_state"] == fs)] ** 2)))
                          for fs, name in zx.FINAL_STATES.items()}
    pred_os = {key: os_prediction(cr_mc, rates_os, window) for key, window in WINDOWS.items()}
    pred_ss = {key: ss_prediction(ss_mc, rates_ss, os_ss, window) for key, window in WINDOWS.items()}
    rng = np.random.default_rng(args.seed)
    boot_os = {key: {name: [] for name in zx.FINAL_STATES.values()} for key in WINDOWS}
    boot_ss = {key: {name: [] for name in zx.FINAL_STATES.values()} for key in WINDOWS}
    for _ in range(args.bootstrap):
        b_zl = rng.poisson(1.0, len(zl_mc["w"])).astype(float)
        b_cr = rng.poisson(1.0, len(cr_mc["w"])).astype(float)
        b_ss = rng.poisson(1.0, len(ss_mc["w"])).astype(float)
        r_os = zx.fake_rates(zl_mc, none, zx.SS_WINDOWS["os_window"], b_zl)
        r_ss, _ = zx.ss_rates(zl_mc, none, dict(ss_high, boot=ss_high["w"]), b_zl)
        for key, window in WINDOWS.items():
            o = os_prediction(cr_mc, r_os, window, b_cr)
            s = ss_prediction(ss_mc, r_ss, os_ss, window, b_ss)
            for name in zx.FINAL_STATES.values():
                boot_os[key][name].append(o[name]["yield"])
                boot_ss[key][name].append(s[name])
    mc_closure = {}
    for key in WINDOWS:
        mc_closure[key] = {}
        for name in zx.FINAL_STATES.values():
            t, te = truth[key][name], truth_err[key][name]
            po, pe = pred_os[key][name]["yield"], float(np.std(boot_os[key][name], ddof=1))
            so, se = pred_ss[key][name], float(np.std(boot_ss[key][name], ddof=1))
            mc_closure[key][name] = {"sr_mc": t, "sr_mc_stat": te,
                                     "os": {"value": po, "stat": pe, "terms": pred_os[key][name],
                                            "ratio_to_sr": po / t if t > 0 else None,
                                            "pull": (po - t) / np.hypot(pe, te) if np.hypot(pe, te) > 0 else None},
                                     "ss": {"value": so, "stat": se, "ratio_to_sr": so / t if t > 0 else None,
                                            "pull": (so - t) / np.hypot(se, te) if np.hypot(se, te) > 0 else None}}

    # --- 2. 3P1F from 2P2F, data and MC ------------------------------------------------------------------
    zl = io.read(scan, zx.DATA, "ZL", zx.ZL_BRANCHES)
    zl_prompt = io.read(scan, zx.PROMPT_MC, "ZL", zx.ZL_BRANCHES)
    cr_all = io.read(scan, zx.DATA, "CR", zx.CR_BRANCHES)
    cr_zz_all = io.read(scan, zx.PROMPT_MC, "CR", zx.CR_BRANCHES)
    cr = sel(cr_all, (cr_all["region"] == 0) & (cr_all["cr_type"] <= 1))
    cr_zz = sel(cr_zz_all, (cr_zz_all["region"] == 0) & (cr_zz_all["cr_type"] == 1))
    rates_data = zx.fake_rates(zl, zl_prompt, zx.SS_WINDOWS["os_window"])
    three = {}
    for label, rows, rates, prompt in (("data", cr, rates_data, cr_zz), ("mc_dy_ttbar", cr_mc, rates_os, None)):
        predicted = predict_3p1f(rows, rates)
        three[label] = {}
        for fs, name in zx.FINAL_STATES.items():
            in_fs = (rows["cr_type"] == 1) & (rows["final_state"] == fs)
            observed = float(np.sum(rows["w"][in_fs]))
            observed_err = float(np.sqrt(np.sum(rows["w"][in_fs] ** 2)))
            zz = float(np.sum(prompt["w"][prompt["final_state"] == fs])) if prompt is not None else 0.0
            three[label][name] = {"predicted_from_2p2f": predicted[name], "observed_3p1f": observed, "observed_stat": observed_err,
                                  "prompt_zz_3p1f": zz, "observed_minus_zz": observed - zz,
                                  "ratio": (observed - zz) / predicted[name] if predicted[name] > 0 else None}

    out_dir.mkdir(parents=True)
    plots = out_dir / "plots"
    plots.mkdir()
    lumi = io.lumi_fb(scan)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, key in zip(axes, WINDOWS):
        x = np.arange(3)
        names = list(zx.FINAL_STATES.values())
        c = mc_closure[key]
        ax.errorbar(x - 0.15, [c[n]["sr_mc"] for n in names], yerr=[c[n]["sr_mc_stat"] for n in names], fmt="ks", label="DY + TTbar MC, SR")
        ax.errorbar(x, [c[n]["os"]["value"] for n in names], yerr=[c[n]["os"]["stat"] for n in names], fmt="bo", label="OS method on MC")
        ax.errorbar(x + 0.15, [c[n]["ss"]["value"] for n in names], yerr=[c[n]["ss"]["stat"] for n in names], fmt="r^", label="SS method on MC")
        ax.set_xticks(x, names)
        ax.set_ylabel(f"events ({lumi} fb-1)")
        ax.set_title(f"Z+X MC closure, {key.replace('_', ' ')}", fontsize=9)
        ax.axhline(0, color="grey", lw=0.5)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plots / "zx_mc_closure.png", dpi=120)
    fig.savefig(plots / "zx_mc_closure.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    x = np.arange(3)
    names = list(zx.FINAL_STATES.values())
    for off, label, marker in ((-0.1, "data", "ko"), (0.1, "mc_dy_ttbar", "bs")):
        t = three[label]
        ratio = [t[n]["ratio"] if t[n]["ratio"] is not None else np.nan for n in names]
        err = [t[n]["observed_stat"] / t[n]["predicted_from_2p2f"] if t[n]["predicted_from_2p2f"] > 0 else np.nan for n in names]
        ax.errorbar(x + off, ratio, yerr=err, fmt=marker, label=label)
    ax.axhline(1.0, color="grey", ls="--")
    ax.set_xticks(x, names)
    ax.set_ylabel("observed 3P1F (- ZZ) / predicted from 2P2F")
    ax.set_title("3P1F prediction from 2P2F with the fake rates", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(plots / "zx_3p1f_from_2p2f.png", dpi=120)
    fig.savefig(plots / "zx_3p1f_from_2p2f.pdf")
    plt.close(fig)
    report = {"schema": "h4l_v3_zx_closure/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
              "select_scan_plan_sha256": scan["plan_sha256"], "luminosity_fb": lumi, "bootstrap": args.bootstrap, "seed": args.seed,
              "mc_samples": zx.FAKE_MC, "mc_fake_rates_os": {str(k): v.tolist() for k, v in rates_os.items()},
              "mc_fake_rates_ss": {str(k): v.tolist() for k, v in rates_ss.items()}, "mc_ss_conversion_correction": ss_info,
              "mc_os_ss_ratio": os_ss,
              "mc_closure": mc_closure, "three_from_two": three}
    (out_dir / "closure.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    for key in WINDOWS:
        for name in zx.FINAL_STATES.values():
            c = mc_closure[key][name]
            print(f"[closure] {key} {name}: SR MC {c['sr_mc']:.2f} +- {c['sr_mc_stat']:.2f}; OS {c['os']['value']:.2f} +- {c['os']['stat']:.2f}; "
                  f"SS {c['ss']['value']:.2f} +- {c['ss']['stat']:.2f}", flush=True)
    for label in three:
        for name, t in three[label].items():
            print(f"[closure] 3P1F {label} {name}: predicted {t['predicted_from_2p2f']:.1f}, observed {t['observed_3p1f']:.1f} "
                  f"(ZZ {t['prompt_zz_3p1f']:.1f}), ratio {t['ratio']}", flush=True)
    print(f"[closure] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
