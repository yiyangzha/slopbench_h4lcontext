"""Per-event mass-uncertainty calibration lambda (AN-16-442 5.3.1, stage 3b).

    pixi run py -- analysis_v3/calibration/scripts/run_lambda.py --version v3 --run nominal_r3 --label lambda_v2

Inputs: the FullPairs trees of the final calibration payload (data
corrected, MC smeared), written by diagnose_calibration.py --final into
production_v3/calibration/<version>/<run>/diag_final/hist/hist/, and the
configuration analysis_v3/calibration/config/lambda_ul16_v2.json.
 1. lambda_histograms per FullPairs file: the Z mass per pair class (the
    lambda regions a <= b of the two legs) and bin of the predicted relative
    mass uncertainty e, with the sums of w, w e, w e^2, w d_a^2, w d_b^2 of
    every bin in 80-100 GeV;
 2. fit_lambda per flavour: every configured fit (mode "same": both legs in
    the target region; "reference": one leg in a reference region with its
    lambda fixed) fits all e bins of the class simultaneously with
    BW (x) DCB + exponential, shared tails and sigma_k = (m_Z / 2)
    sqrt(lambda_a^2 <d_a^2>_k + lambda_b^2 <d_b^2>_k); MC first, data with the
    MC tails; plus per-bin closure fits with sigma free;
 3. lambda.json (lambda per region for data and MC, the payload used
    downstream: every lepton's relative momentum error is multiplied by the
    lambda of its region), the closure plots (measured against predicted
    sigma before and after the correction, AN Figure 18) and a gallery of
    every simultaneous fit.
Outputs in production_v3/calibration/<version>/<run>/<label>/: run.json,
hist/, fits/, lambda.json, plots/.  Nothing is overwritten; finished
identical jobs are reused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402
from run_calibration import PRODUCTION, canonical, fail, parallel, run_job, sha256_file, utc, write_new  # noqa: E402

LAMBDA_CONFIG = REPO / "analysis_v3/calibration/config/lambda_ul16_v2.json"
PROGRAMS = {name: REPO / "analysis_v3/calibration/bin" / name for name in ("lambda_histograms", "fit_lambda")}
FLAVOURS = {"mm": "muon", "ee": "electron"}


def freeze() -> dict:
    digests = {name: sha256_file(path) for name, path in PROGRAMS.items()}
    key = hashlib.sha256("".join(digests[n] for n in sorted(digests)).encode()).hexdigest()[:12]
    directory = PRODUCTION / "program_inputs" / "lambda_v2" / key
    directory.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for name, source in PROGRAMS.items():
        target = directory / name
        if not target.exists():
            temporary = target.with_name(name + f".partial.{os.getpid()}")
            shutil.copy2(source, temporary)
            if sha256_file(temporary) != digests[name]:
                fail(f"frozen copy of {name} differs")
            temporary.rename(target)
        elif sha256_file(target) != digests[name]:
            fail(f"frozen {target} differs from its content address")
        os.chmod(target, 0o755)
        frozen[name] = {"path": str(target), "sha256": digests[name], "source": str(source)}
    return {"key": key, "programs": frozen}


def closure_plot(report: dict, name: str, stem: Path) -> None:
    """Measured (per-bin closure fit) against predicted sigma, before and after lambda."""
    fits = report["fits"]
    fig, axes = plt.subplots(2, len(fits), figsize=(3.6 * len(fits), 7.2), squeeze=False)
    for column, fit in enumerate(fits):
        region = report["regions"][fit["target"]]["name"]
        for row, role in enumerate(("mc", "data")):
            ax = axes[row][column]
            result = fit.get(role, {})
            ax.set_title(f"{'MC' if role == 'mc' else 'data'}: {region}"
                         + (f" (ref {fit['reference']})" if fit["mode"] == "reference" else ""), fontsize=8)
            ax.set_xlabel("predicted sigma(m_ll) [GeV]", fontsize=8)
            ax.set_ylabel("fitted sigma(m_ll) [GeV]", fontsize=8)
            ax.tick_params(labelsize=7)
            if not result.get("categories"):
                ax.text(0.5, 0.5, result.get("skipped", "no fit"), ha="center", va="center", fontsize=8, transform=ax.transAxes)
                continue
            cats = [c for c in result["categories"] if c["closure"]["ok"]]
            before = np.array([c["predicted_sigma_before"] for c in cats])
            after = np.array([c["predicted_sigma_after"] for c in cats])
            measured = np.array([c["closure"]["sigma"] for c in cats])
            # A closure fit without a valid error (null in the report) is drawn without an error bar.
            error = np.array([c["closure"]["sigma_error"] if c["closure"]["sigma_error"] is not None else np.nan for c in cats],
                             dtype=float)
            lam = result["lambda"]["value"]
            ax.errorbar(before, measured, error, fmt="o", ms=3, color="#5790fc", label="uncorrected")
            ax.errorbar(after, measured, error, fmt="s", ms=3, color="#e42536",
                        label=f"corrected, lambda {lam:.3f} +- {result['lambda']['error']:.3f}")
            top = 1.2 * max(np.max(measured), np.max(after), np.max(before)) if len(cats) else 1
            grid = np.linspace(0, top, 10)
            ax.plot(grid, grid, color="grey", lw=0.8)
            ax.fill_between(grid, 0.8 * grid, 1.2 * grid, color="grey", alpha=0.15, label="+-20 %")
            ax.set_xlim(0, top)
            ax.set_ylim(0, top)
            ax.legend(fontsize=6, loc="upper left")
    fig.suptitle(f"per-event mass uncertainty closure, {name} pairs", fontsize=10)
    fig.tight_layout()
    h4l_style.save(fig, stem)
    plt.close(fig)


def gallery(report: dict, name: str, path: Path) -> int:
    """Every category of every simultaneous fit, data and MC, with pulls."""
    panels = []
    for fit in report["fits"]:
        for role in ("mc", "data"):
            for cat in fit.get(role, {}).get("categories", []):
                panels.append((fit, role, cat))
    with PdfPages(path) as pdf:
        for start in range(0, len(panels), 12):
            fig = plt.figure(figsize=(16, 12))
            grid = fig.add_gridspec(6, 4, height_ratios=[3, 1] * 3, hspace=0.08, wspace=0.25)
            for slot, (fit, role, cat) in enumerate(panels[start:start + 12]):
                r, c = divmod(slot, 4)
                top = fig.add_subplot(grid[2 * r, c])
                bottom = fig.add_subplot(grid[2 * r + 1, c], sharex=top)
                curve = cat["curve"]
                x, y, ey, model = (np.array(curve[k]) for k in ("x", "y", "ey", "model"))
                top.errorbar(x, y, ey, fmt="o", ms=1.2, color="black", lw=0.5, label=role)
                top.plot(x, model, color="#e42536", lw=1, label="BW (x) DCB + exp")
                top.plot(x, np.array(curve["background"]), color="#5790fc", lw=0.8, ls="--", label="background")
                top.set_yscale("log")
                top.set_ylim(max(0.5, 0.3 * np.min(model[model > 0])) if np.any(model > 0) else 0.5, 3 * np.max(y))
                result = fit[role]
                top.set_title(f"{report['regions'][fit['target']]['name']} {role} bin {cat['bin']}: sigma "
                              f"{cat['predicted_sigma_after']:.2f} (lambda {result['lambda']['value']:.3f})"
                              + ("" if result["ok"] else " NOT OK"), fontsize=7)
                top.tick_params(labelsize=6, labelbottom=False)
                pull = np.divide(y - model, ey, out=np.zeros_like(y), where=ey > 0)
                bottom.plot(x, pull, ".", ms=1.2, color="black")
                bottom.axhline(0, color="grey", lw=0.5)
                bottom.set_ylim(-5, 5)
                bottom.tick_params(labelsize=6)
                if slot == 0:
                    top.legend(fontsize=6)
            fig.suptitle(f"lambda fits, {name} (m_ll [GeV])", fontsize=11)
            pdf.savefig(fig)
            plt.close(fig)
    return len(panels)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--config", type=Path, default=LAMBDA_CONFIG, help="lambda configuration (default: %(default)s)")
    args = parser.parse_args()
    config_path = args.config if args.config.is_absolute() else REPO / args.config
    run_dir = PRODUCTION / "calibration" / args.version / args.run
    source = run_dir / "diag_final" / "hist" / "hist"
    base = run_dir / args.label
    lam = json.loads(config_path.read_text(encoding="utf-8"))
    inputs = {"data": [], "mc": []}
    for path in sorted(source.iterdir()):
        match = re.fullmatch(r"(data|mc)_calx_.+\.root", path.name)
        if not match:
            continue
        done = path.with_name(path.name + ".done.json")
        if not done.exists():
            fail(f"{path} has no .done.json record")
        record = json.loads(done.read_text(encoding="utf-8"))
        if record["output_sha256"] != sha256_file(path):
            fail(f"{path} differs from its .done.json record")
        with uproot.open(path) as f:
            entries = f["FullPairs"].num_entries
        inputs[match.group(1)].append({"path": str(path), "full_pairs": int(entries), "root_sha256": record["output_sha256"]})
    if not inputs["data"] or not inputs["mc"]:
        fail(f"no FullPairs files in {source}")
    frozen = freeze()
    record = {"source": str(source), "inputs": inputs, "lambda_config": str(config_path),
              "lambda_config_sha256": sha256_file(config_path), "frozen": frozen,
              "orchestrator_sha256": sha256_file(Path(__file__))}
    run_path = base / "run.json"
    if run_path.exists():
        existing = json.loads(run_path.read_text(encoding="utf-8"))
        if canonical({k: v for k, v in existing.items() if k != "started_utc"}) != canonical(record):
            fail(f"{run_path} records other inputs or programs; use a new --label")
    else:
        write_new(run_path, dict(record, started_utc=utc()))

    print(f"[{utc()}] histograms", flush=True)
    tasks, files = [], {"data": [], "mc": []}
    for role, items in inputs.items():
        for item in items:
            name = Path(item["path"]).stem
            job = {"lambda_config": lam, "weighted": role == "mc",
                   "inputs": [{"path": item["path"], "full_pairs": item["full_pairs"]}]}
            out = base / "hist" / f"{name}.root"
            tasks.append((run_job, frozen["programs"]["lambda_histograms"], job, base / "jobs" / f"{name}.json", out,
                          base / "logs" / f"{name}.log"))
            files[role].append(str(out))
    parallel(tasks, args.workers)

    print(f"[{utc()}] fits", flush=True)
    tasks = []
    for tag in FLAVOURS:
        job = {"lambda_config": lam, "flavour": tag, "data": files["data"], "mc": files["mc"], "curves": True}
        tasks.append((run_job, frozen["programs"]["fit_lambda"], job, base / "jobs" / f"fit_{tag}.json",
                      base / "fits" / f"fit_{tag}.json", base / "logs" / f"fit_{tag}.log"))
    records = parallel(tasks, 2)

    plots = base / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    payload = {"schema": "h4l_v3_lambda_payload/1", "version": lam["version"], "run": str(run_dir),
               "convention": "a lepton's relative momentum error d = dpT/pT is multiplied by lambda of its region "
                             "(first matching region; data and MC each with their own lambda)",
               "flavours": {}}
    lines = []
    for rec, (tag, name) in zip(records, FLAVOURS.items()):
        report = json.loads(Path(rec["output"]).read_text(encoding="utf-8"))
        closure_plot(report, name, plots / f"closure_{tag}")
        panels = gallery(report, name, plots / f"gallery_{tag}.pdf")
        regions = []
        for index, region in enumerate(report["regions"]):
            fit = next((f for f in report["fits"] if f["target"] == index), None)
            entry = {"region": region, "fit_mode": fit["mode"] if fit else None}
            for role in ("data", "mc"):
                result = fit.get(role) if fit else None
                if result and result.get("ok"):
                    entry[role] = {"lambda": result["lambda"]["value"], "error": result["lambda"]["error"],
                                   "chi2": result["chi2"], "ndf": result["ndf"], "categories": result["categories_used"]}
                else:
                    entry[role] = None
                    fail(f"{name} region {region['name']}: no usable {role} lambda fit "
                         f"({result.get('skipped', 'fit not ok') if result else 'no fit'}); revise the configuration")
            regions.append(entry)
            lines.append(f"{name} {region['name']}: data {entry['data']['lambda']:.4f} +- {entry['data']['error']:.4f}, "
                         f"MC {entry['mc']['lambda']:.4f} +- {entry['mc']['error']:.4f} ({entry['fit_mode']})")
        payload["flavours"][name] = {"regions": regions, "gallery_panels": panels}
    write_new(base / "lambda.json", payload)
    for line in lines:
        print(f"[lambda] {line}")
    print(f"[{utc()}] {base / 'lambda.json'}; plots {plots}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
