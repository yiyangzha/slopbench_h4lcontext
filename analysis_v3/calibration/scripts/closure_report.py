"""Closure report of one calibration closure run (stage 3a).

    pixi run py -- analysis_v3/calibration/scripts/closure_report.py --version v3 --run closure_sloped_same

The run's data role is MC with the injection of closure_ul16_v1.json:
per lepton pT *= (1 + s(pT, |eta|)) (1 + r(pT, |eta|) N), s and r the
documented test functions, whose value at the report point (45 GeV, 1.2) is
their constant term.  This script compares, per flavour,
  * the extracted scale shift and smear at the report point (the converged
    payload plus the last residual, statistical uncertainty from the joint
    covariance) with the injected values;
  * the extracted ln(1 + s) and r^2 at every fine |eta| bin centre and pT node
    with the injected functions there;
and writes <run>/closure/closure.json and the profile plots
<run>/closure/closure_<flavour>.{pdf,png}.
"""

from __future__ import annotations

import argparse
import hashlib
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
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
import h4l_style  # noqa: E402
from plot_calibration import evaluate, full_covariance  # noqa: E402
from run_calibration import CLOSURE, Model  # noqa: E402

FLAVOURS = {"mm": "muon", "ee": "electron"}


def test_function(node: dict, pt: float, abs_eta: float) -> float:
    x = abs_eta - node.get("eta0", 1.2)
    return (node.get("a", 0.0) + node.get("b_eta", 0.0) * x + node.get("c_pt", 0.0) * (pt - node.get("pt0", 45.0)) / node.get("pt0", 45.0)
            + node.get("d_eta2", 0.0) * x * x)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--closure-config", type=Path, default=CLOSURE,
                        help="the closure configuration of the run (its sha256 must equal the run record)")
    args = parser.parse_args()
    base = PRODUCTION / "calibration" / args.version / args.run
    run = json.loads((base / "run.json").read_text(encoding="utf-8"))
    if not run.get("closure"):
        print("ERROR: not a closure run", file=sys.stderr)
        return 1
    closure_path = args.closure_config if args.closure_config.is_absolute() else REPO / args.closure_config
    recorded = run.get("closure_config_sha256")
    if recorded and hashlib.sha256(closure_path.read_bytes()).hexdigest() != recorded:
        print(f"ERROR: {closure_path} is not the closure configuration of the run", file=sys.stderr)
        return 1
    scenario = json.loads(closure_path.read_text(encoding="utf-8"))["scenarios"][run["closure"]]
    summary_path = base / "summary_final.json" if (base / "summary_final.json").exists() else base / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not (summary.get("converged") or summary.get("finalized_by_average")):
        print("ERROR: the closure run neither converged nor was finalized", file=sys.stderr)
        return 1
    config = json.loads((base / "config.json").read_text(encoding="utf-8"))
    final = json.loads((base / "payload.json").read_text(encoding="utf-8"))
    solution = json.loads((base / f"iter_{summary['converged_iteration']:02d}" / "solution.json").read_text(encoding="utf-8"))
    point = config["report_point"]
    out_dir = base / "closure"
    out_dir.mkdir(exist_ok=True)
    report = {"run": args.run, "scenario": run["closure"], "mode": scenario["mode"], "flavours": {}}
    for tag, name in FLAVOURS.items():
        model = Model(config["model"][name])
        payload = final["flavours"][name]
        inject = scenario["inject"][name]
        cov_u = full_covariance(solution[tag]["scale"], model.n_par)
        cov_v = full_covariance(solution[tag]["smear"], model.n_par)
        u, v, grad = evaluate(model, payload, point["pt"], point["abs_eta"])
        u += float(grad @ np.array(solution[tag]["scale"]["value"]))
        v += float(grad @ np.array(solution[tag]["smear"]["value"]))
        rms = final.get("iteration_rms", {}).get(tag)
        u_rms2 = float(np.sum((grad * np.array(rms["scale"])) ** 2)) if rms else 0.0
        v_rms2 = float(np.sum((grad * np.array(rms["smear"])) ** 2)) if rms else 0.0
        s_err = (1 + math.expm1(u)) * math.sqrt(float(grad @ cov_u @ grad) + u_rms2)
        v_err = math.sqrt(float(grad @ cov_v @ grad) + v_rms2)
        s_true = test_function(inject["scale"], point["pt"], point["abs_eta"])
        r_true = max(0.0, test_function(inject["smear"], point["pt"], point["abs_eta"]))
        r_value = math.sqrt(max(v, 0.0))
        entry = {"scale": {"measured": math.expm1(u), "stat": s_err, "injected": s_true,
                           "pull": (math.expm1(u) - s_true) / s_err if s_err > 0 else None},
                 "smear_variance": {"measured": v, "stat": v_err, "injected": r_true ** 2,
                                    "pull": (v - r_true ** 2) / v_err if v_err > 0 else None},
                 "smear": {"measured": r_value, "injected": r_true}}
        # Profiles at the fine |eta| bin centres and the pT nodes.
        profile = []
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        ref_node = payload["node_pt"][0][model.ref]
        etas, measured_s, true_s, errors_s, measured_v, true_v, errors_v = [], [], [], [], [], [], []
        for e in range(model.n_eta):
            upper = model.eta[e + 1] if e + 1 < model.n_eta else model.eta[e] + 0.2
            centre = 0.5 * (model.eta[e] + upper)
            ue, ve, ge = evaluate(model, payload, ref_node, centre)
            ue += float(ge @ np.array(solution[tag]["scale"]["value"]))
            ve += float(ge @ np.array(solution[tag]["smear"]["value"]))
            etas.append(centre)
            measured_s.append(math.expm1(ue))
            errors_s.append(math.sqrt(float(ge @ cov_u @ ge)))
            true_s.append(test_function(inject["scale"], ref_node, centre))
            measured_v.append(ve)
            errors_v.append(math.sqrt(float(ge @ cov_v @ ge)))
            true_v.append(max(0.0, test_function(inject["smear"], ref_node, centre)) ** 2)
            profile.append({"abs_eta": centre, "pt": ref_node, "s": measured_s[-1], "s_true": true_s[-1],
                            "v": ve, "v_true": true_v[-1]})
        axes[0, 0].errorbar(etas, 100 * np.array(measured_s), 100 * np.array(errors_s), fmt="o", label="extracted")
        axes[0, 0].plot(etas, 100 * np.array(true_s), "r-", label="injected")
        axes[0, 0].set_ylabel(f"s at {ref_node:.0f} GeV [%]")
        axes[1, 0].errorbar(etas, 1e5 * np.array(measured_v), 1e5 * np.array(errors_v), fmt="o", label="extracted")
        axes[1, 0].plot(etas, 1e5 * np.array(true_v), "r-", label="injected")
        axes[1, 0].set_ylabel(r"$r^2 \times 10^5$")
        for ax in axes[:, 0]:
            ax.set_xlabel("|eta| bin centre")
            ax.legend()
        for r in range(model.n_regions):
            centre = 0.5 * (model.regions[r] + (model.regions[r + 1] if r + 1 < model.n_regions else 2.4))
            nodes = payload["node_pt"][r]
            ms, ts, mv, tv = [], [], [], []
            for pt in nodes:
                ue, ve, ge = evaluate(model, payload, pt, centre)
                ue += float(ge @ np.array(solution[tag]["scale"]["value"]))
                ve += float(ge @ np.array(solution[tag]["smear"]["value"]))
                ms.append(math.expm1(ue))
                ts.append(test_function(inject["scale"], pt, centre))
                mv.append(ve)
                tv.append(max(0.0, test_function(inject["smear"], pt, centre)) ** 2)
                profile.append({"abs_eta": centre, "pt": pt, "s": ms[-1], "s_true": ts[-1], "v": ve, "v_true": tv[-1]})
            axes[0, 1].plot(nodes, 100 * np.array(ms), "o", label=f"extracted, |eta| {centre:.2f}")
            axes[0, 1].plot(nodes, 100 * np.array(ts), "-", label=f"injected, |eta| {centre:.2f}")
            axes[1, 1].plot(nodes, 1e5 * np.array(mv), "o", label=f"extracted, |eta| {centre:.2f}")
            axes[1, 1].plot(nodes, 1e5 * np.array(tv), "-", label=f"injected, |eta| {centre:.2f}")
        for ax in axes[:, 1]:
            ax.set_xscale("log")
            ax.set_xlabel("pT node [GeV]")
            ax.legend(fontsize="x-small")
        axes[0, 1].set_ylabel("s [%]")
        axes[1, 1].set_ylabel(r"$r^2 \times 10^5$")
        fig.suptitle(f"closure {run['closure']} ({scenario['mode']}), {name}")
        h4l_style.save(fig, out_dir / f"closure_{tag}")
        entry["profile"] = profile
        report["flavours"][name] = entry
    (out_dir / "closure.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({n: {k: v for k, v in f.items() if k != "profile"} for n, f in report["flavours"].items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
