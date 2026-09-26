"""Plots and report values of one calibration run (stage 3a, factorized model).

    pixi run py -- analysis_v3/calibration/scripts/plot_calibration.py --version v3 --run nominal [--plot-dir plots]

Reads production_v3/calibration/<version>/<run>/ (iterations, solutions,
payloads, template fits, histograms) and writes into <run>/<plot dir>/:
  convergence                    worst residual and per-lepton chi2/ndf per
                                 iteration, and the scale and smear at the
                                 report point after each iteration;
  profile_eta_<flavour>          s and r along the fine |eta| bins at the
                                 reference pT (the a and c terms);
  profile_pt_<flavour>           the pT terms b_R and d_R per region;
  pulls_<flavour>                per-category residual pulls of the final
                                 joint least squares, per family;
  response_<flavour>             the per-category scale and smear responses;
  gallery_<family>_<flavour>_iter<n>.pdf   every fitted category of iterations
                                 0 and final: data, fitted template and
                                 untransformed template, with the residuals;
  zpeak_<flavour>_{before,after} the inclusive dilepton mass, data against MC,
                                 before (iteration 0) and after (final);
  kinematics_<family>_<flavour>.pdf  per category leg pT and pair pT, data
                                 against MC, final iteration;
and report.json: s and r at the report point (45 GeV, |eta| 1.2) with the
statistical uncertainty propagated from the joint covariance, and every term.
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
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
sys.path.insert(0, str(REPO / "analysis_v3/calibration/scripts"))
import h4l_style  # noqa: E402
from run_calibration import Model  # noqa: E402

FLAVOURS = {"mm": "muon", "ee": "electron"}
LABELS = {"mm": r"$\mu\mu$", "ee": r"$ee$"}
FAMILIES = ("A", "B")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_fits(directory: Path, family: str, tag: str) -> list:
    results = []
    for path in sorted(directory.glob(f"tfit_{family}_{tag}_*.json")):
        if path.name.endswith(".done.json"):
            continue
        results.extend(json.loads(path.read_text(encoding="utf-8"))["results"])
    return results


def lower_bin(edges, value):
    index = -1
    for k, edge in enumerate(edges):
        if value >= edge:
            index = k
    return index


def pt_weights(nodes, pt):
    """Linear interpolation weights along the pT nodes (constant beyond)."""
    w = np.zeros(len(nodes))
    if pt <= nodes[0]:
        w[0] = 1
    elif pt >= nodes[-1]:
        w[-1] = 1
    else:
        i = max(k for k in range(len(nodes) - 1) if nodes[k] <= pt)
        t = (pt - nodes[i]) / (nodes[i + 1] - nodes[i])
        w[i], w[i + 1] = 1 - t, t
    return w


def evaluate(model: Model, payload: dict, pt: float, abs_eta: float):
    """u, v and their gradients in the parameter vector of the joint least squares."""
    e = lower_bin(model.eta, abs_eta)
    r = model.region_of_eta[e]
    weights = pt_weights(payload["node_pt"][r], pt)
    grad = np.zeros(model.n_par)
    grad[e] = 1.0
    for p, weight in enumerate(weights):
        if p != model.ref and weight:
            grad[model.b_index[(r, p)]] += weight
    u = payload["a"][e] + float(weights @ np.array(payload["b"][r]))
    v = payload["c"][e] + float(weights @ np.array(payload["d"][r]))
    return u, v, grad


def full_covariance(solution: dict, n_par: int) -> np.ndarray:
    cov = np.zeros((n_par, n_par))
    parameters = solution["covariance"]["parameters"]
    matrix = np.array(solution["covariance"]["matrix"])
    for a, i in enumerate(parameters):
        for b, j in enumerate(parameters):
            cov[i, j] = matrix[a, b]
    return cov


def gallery(fits: list, title: str, path: Path) -> int:
    fitted = [e for e in fits if "fit" in e]
    with PdfPages(path) as pdf:
        for start in range(0, len(fitted), 12):
            fig = plt.figure(figsize=(16, 12))
            grid = fig.add_gridspec(6, 4, height_ratios=[3, 1] * 3, hspace=0.05, wspace=0.25)
            for slot, entry in enumerate(fitted[start:start + 12]):
                r, c = divmod(slot, 4)
                top = fig.add_subplot(grid[2 * r, c])
                bottom = fig.add_subplot(grid[2 * r + 1, c], sharex=top)
                curve = entry["fit"]["curve"]
                lo, hi = np.array(curve["x_lo"]), np.array(curve["x_hi"])
                width, x = hi - lo, 0.5 * (lo + hi)
                y, ey, model, before = (np.array(curve[k]) for k in ("y", "ey", "model", "template_before"))
                top.errorbar(x, y / width, ey / width, fmt="o", ms=1.5, color="black", lw=0.6, label="data")
                top.stairs(model / width, np.append(lo, hi[-1]), color="#e42536", label="fitted template")
                top.stairs(before / width, np.append(lo, hi[-1]), color="#5790fc", ls="--", label="MC, k=1, D=0")
                top.set_yscale("log")
                fit = entry["fit"]
                top.set_title(f"{entry['category']}: ln k {fit['lnk']['value']:.5f}, E {fit['E']['value']:.2e}"
                              + ("" if entry["ok"] else " NOT OK"), fontsize=8)
                top.tick_params(labelsize=7, labelbottom=False)
                pull = np.divide(y - model, ey, out=np.zeros_like(y), where=ey > 0)
                bottom.plot(x, pull, ".", ms=1.5, color="black")
                bottom.axhline(0, color="grey", lw=0.5)
                bottom.set_ylim(-5, 5)
                bottom.tick_params(labelsize=7)
                if slot == 0:
                    top.legend(fontsize=6)
            fig.suptitle(title, fontsize=11)
            pdf.savefig(fig)
            plt.close(fig)
    return len(fitted)


def summed(files: list, prefix: str, rebin: int = 5):
    total = variance = edges = None
    for path in files:
        with uproot.open(path) as f:
            for name in f.keys(recursive=False, cycle=False):
                if not name.startswith(prefix):
                    continue
                values, e = f[name].to_numpy()
                var = f[name].variances()
                if total is None:
                    total, variance, edges = values.copy(), var.copy(), e
                else:
                    total += values
                    variance += var
    if total is None:
        return None
    n = len(total) // rebin * rebin
    return total[:n].reshape(-1, rebin).sum(axis=1), variance[:n].reshape(-1, rebin).sum(axis=1), edges[: n + 1: rebin]


def zpeak_plot(data, mc, title: str, lumi: float, path: Path) -> dict:
    d, dv, edges = data
    m, _, _ = mc
    scale = d.sum() / m.sum()
    fig, (top, bottom) = h4l_style.figure(ratio=True)
    centres = 0.5 * (edges[1:] + edges[:-1])
    top.stairs(m * scale, edges, color="#5790fc", label=f"MC (normalized x{scale:.3f})", fill=True, alpha=0.5)
    top.errorbar(centres, d, np.sqrt(dv), fmt="o", ms=2.5, color="black", label="data")
    top.set_ylabel("pairs / 0.25 GeV")
    top.set_xlim(60, 120)
    top.legend(fontsize="small", title=title, title_fontsize="small")
    ratio = np.divide(d, m * scale, out=np.full_like(d, np.nan), where=m > 0)
    err = np.divide(np.sqrt(dv), m * scale, out=np.full_like(d, np.nan), where=m > 0)
    bottom.errorbar(centres, ratio, err, fmt="o", ms=2, color="black")
    bottom.axhline(1, color="grey", ls="--")
    bottom.set_ylim(0.5, 1.5)
    bottom.set_xlim(60, 120)
    bottom.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    bottom.set_ylabel("data/MC")
    h4l_style.label(top, data=True, lumi_fb=lumi)
    h4l_style.save(fig, path)
    window = (centres > 60) & (centres < 120)
    return {"data_60_120": float(d[window].sum()), "mc_60_120_normalized": float((m * scale)[window].sum()),
            "mc_normalization": float(scale)}


def collect(files: list, prefixes: tuple) -> dict:
    """Histograms whose names start with one of the prefixes, summed over the
    files, each file read once: {name: [values, variances, edges]}."""
    out = {}
    for path in files:
        with uproot.open(path) as f:
            for name in f.keys(recursive=False, cycle=False):
                if not name.startswith(prefixes):
                    continue
                h = f[name]
                values, edges = h.to_numpy()
                variances = h.variances()
                if name in out:
                    out[name][0] = out[name][0] + values
                    out[name][1] = out[name][1] + variances
                else:
                    out[name] = [values.copy(), variances.copy(), edges]
    return out


def kinematics(hist_dir: Path, family: str, tag: str, path: Path) -> int:
    """Per category leg pT and pair pT, data against MC (normalized)."""
    prefixes = (f"lpt{family}_{tag}_", f"ppt{family}_{tag}_")
    data = collect(sorted(str(p) for p in hist_dir.glob("data_*.root")), prefixes)
    mc = collect(sorted(str(p) for p in hist_dir.glob("mc_*.root")), prefixes)
    names = sorted({n[len(prefixes[0]):] for n in data if n.startswith(prefixes[0])},
                   key=lambda c: tuple(int(x) for x in c.split("_")))
    count = 0
    with PdfPages(path) as pdf:
        for start in range(0, len(names), 8):
            fig, axes = plt.subplots(4, 4, figsize=(16, 12))
            for slot, category in enumerate(names[start:start + 8]):
                for column, prefix in enumerate(prefixes):
                    ax = axes[slot // 2, 2 * (slot % 2) + column]
                    d, m = data.get(prefix + category), mc.get(prefix + category)
                    if d is None or m is None or not m[0].sum() > 0:
                        continue
                    scale = d[0].sum() / m[0].sum()
                    ax.stairs(m[0] * scale, m[2], color="#5790fc", fill=True, alpha=0.5, label="MC")
                    ax.errorbar(0.5 * (d[2][1:] + d[2][:-1]), d[0], np.sqrt(d[1]), fmt="o", ms=1.5, color="black", label="data")
                    ax.set_title(f"{family} {category} {'leg pT' if prefix.startswith('lpt') else 'pair pT'}", fontsize=8)
                    ax.tick_params(labelsize=7)
                    count += 1
            pdf.savefig(fig)
            plt.close(fig)
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--plot-dir", default="plots")
    args = parser.parse_args()
    base = PRODUCTION / "calibration" / args.version / args.run
    summary_path = base / "summary_final.json" if (base / "summary_final.json").exists() else base / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not (summary.get("converged") or summary.get("finalized_by_average")):
        fail(f"{base} neither converged nor was finalized")
    final = json.loads((base / "payload.json").read_text(encoding="utf-8"))
    config = json.loads((base / "config.json").read_text(encoding="utf-8"))
    run = json.loads((base / "run.json").read_text(encoding="utf-8"))
    scan = json.loads(Path(run["extract_scan"]).read_text(encoding="utf-8"))
    lumi = scan["samples"]["pseudo_data"]["lumi_fb"]
    last = summary["converged_iteration"]
    plots = base / args.plot_dir
    if plots.exists():
        fail(f"{plots} exists; choose a new --plot-dir")
    plots.mkdir(parents=True)
    point = config["report_point"]
    models = {tag: Model(config["model"][name]) for tag, name in FLAVOURS.items()}
    report = {"run": args.run, "converged_iteration": last, "report_point": point, "flavours": {}}

    # Convergence.
    iterations = [h["iteration"] for h in summary["iterations"]]
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    axes[0].semilogy(iterations, [h["worst_residual_sigma"] for h in summary["iterations"]], "o-")
    axes[0].axhline(config["iteration"]["converged_sigma"], color="grey", ls="--")
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel("worst residual / sigma")
    for tag in FLAVOURS:
        for q, style in (("scale", "o-"), ("smear", "s--")):
            values = [h["chi2_ndf"][tag][q][0] / max(h["chi2_ndf"][tag][q][1], 1) for h in summary["iterations"]]
            axes[1].semilogy(iterations, values, style, label=f"{LABELS[tag]} {q}")
    axes[1].axhline(1, color="grey", ls="--")
    axes[1].set_xlabel("iteration")
    axes[1].set_ylabel("per-lepton chi2/ndf")
    axes[1].legend(fontsize="small")
    for tag, name in FLAVOURS.items():
        s_hist, r_hist = [], []
        for n in iterations:
            payload = json.loads((base / f"iter_{n:02d}" / "payload.json").read_text(encoding="utf-8"))[name]
            u, v, _ = evaluate(models[tag], payload, point["pt"], point["abs_eta"])
            s_hist.append(100 * math.expm1(u))
            r_hist.append(100 * math.sqrt(max(v, 0.0)))
        axes[2].plot(np.array(iterations) + 1, s_hist, "o-", label=LABELS[tag])
        axes[3].plot(np.array(iterations) + 1, r_hist, "o-", label=LABELS[tag])
    axes[2].set_xlabel("payload after iteration")
    axes[2].set_ylabel("scale shift at (45 GeV, 1.2) [%]")
    axes[3].set_xlabel("payload after iteration")
    axes[3].set_ylabel("smear at (45 GeV, 1.2) [%]")
    for ax in axes[2:]:
        ax.legend()
    h4l_style.save(fig, plots / "convergence")

    final_dir = base / f"iter_{last:02d}"
    solution = json.loads((final_dir / "solution.json").read_text(encoding="utf-8"))
    response_path = base / "iter_01" / "response.json"
    response = json.loads(response_path.read_text(encoding="utf-8")) if response_path.exists() else None
    for tag, name in FLAVOURS.items():
        model = models[tag]
        payload = final["flavours"][name]
        cov_u = full_covariance(solution[tag]["scale"], model.n_par)
        cov_v = full_covariance(solution[tag]["smear"], model.n_par)
        # Report point: the applied payload plus the final (vanishing) residual.
        u, v, grad = evaluate(model, payload, point["pt"], point["abs_eta"])
        du = float(grad @ np.array(solution[tag]["scale"]["value"]))
        dv = float(grad @ np.array(solution[tag]["smear"]["value"]))
        u_err = math.sqrt(float(grad @ cov_u @ grad))
        v_err = math.sqrt(float(grad @ cov_v @ grad))
        rms = final.get("iteration_rms", {}).get(tag)
        u_rms = math.sqrt(float(np.sum((grad * np.array(rms["scale"])) ** 2))) if rms else 0.0
        v_rms = math.sqrt(float(np.sum((grad * np.array(rms["smear"])) ** 2))) if rms else 0.0
        s_value = math.expm1(u + du)
        variance = v + dv
        r_value = math.sqrt(max(variance, 0.0))
        report["flavours"][name] = {
            "scale_shift": {"value": s_value, "stat": (1 + s_value) * u_err, "iteration": (1 + s_value) * u_rms,
                            "applied": math.expm1(u), "residual_u": du},
            "smear": {"value": r_value, "stat": v_err / (2 * r_value) if r_value > 0 else None,
                      "iteration": v_rms / (2 * r_value) if r_value > 0 else None},
            "smear_variance": {"value": variance, "stat": v_err, "iteration": v_rms, "applied": v, "residual": dv},
            "gradient_parameters": {str(k): float(g) for k, g in enumerate(grad) if g},
            "payload": payload,
            "least_squares": {q: {k: solution[tag][q][k] for k in ("chi2", "ndf", "inflation", "categories")}
                              for q in ("scale", "smear")}}
        # |eta| profile at the reference pT node, and the pT terms.
        ref_pt = payload["node_pt"][0][model.ref]
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        centres = [0.5 * (model.eta[e] + (model.eta[e + 1] if e + 1 < model.n_eta else model.eta[e] + 0.3))
                   for e in range(model.n_eta)]
        errors_a = [math.sqrt(cov_u[e, e]) for e in range(model.n_eta)]
        errors_c = [math.sqrt(cov_v[e, e]) for e in range(model.n_eta)]
        axes[0].errorbar(centres, [100 * math.expm1(x) for x in payload["a"]], [100 * e for e in errors_a], fmt="o-")
        axes[0].set_ylabel(f"scale shift at the {ref_pt:.0f} GeV node [%]")
        axes[1].errorbar(centres, [100 * math.sqrt(max(x, 0)) for x in payload["c"]],
                         [100 * e / (2 * math.sqrt(x)) if x > 0 else 0 for e, x in zip(errors_c, payload["c"])], fmt="o-")
        axes[1].set_ylabel("extra per-lepton smear at the reference node [%]")
        for ax in axes:
            ax.set_xlabel(r"$|\eta|$" if not model.node["eta_variable"].endswith("sc") else r"$|\eta_{SC}|$")
            ax.set_title(LABELS[tag])
        h4l_style.save(fig, plots / f"profile_eta_{tag}")
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        for r in range(model.n_regions):
            nodes = payload["node_pt"][r]
            eb = [math.sqrt(cov_u[model.b_index[(r, p)], model.b_index[(r, p)]]) if p != model.ref else 0 for p in range(model.n_pt)]
            ed = [math.sqrt(cov_v[model.b_index[(r, p)], model.b_index[(r, p)]]) if p != model.ref else 0 for p in range(model.n_pt)]
            label = f"region {model.regions[r]:g}+"
            axes[0].errorbar(nodes, [1e3 * x for x in payload["b"][r]], [1e3 * e for e in eb], fmt="o-", label=label)
            axes[1].errorbar(nodes, [1e5 * x for x in payload["d"][r]], [1e5 * e for e in ed], fmt="o-", label=label)
        axes[0].set_ylabel(r"$b_R(p_T)$ = pT term of ln(1+s) $\times 10^3$")
        axes[1].set_ylabel(r"$d_R(p_T)$ = pT term of $r^2$ $\times 10^5$")
        for ax in axes:
            ax.set_xscale("log")
            ax.set_xlabel(r"node $p_{\mathrm{T}}$ [GeV]")
            ax.axhline(0, color="grey", lw=0.5)
            ax.legend(fontsize="small", title=LABELS[tag])
        h4l_style.save(fig, plots / f"profile_pt_{tag}")
        # Pulls.
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for ax, quantity in zip(axes, ("scale", "smear")):
            pulls = solution[tag][quantity]["pulls"]
            for family, colour in (("A", "#e42536"), ("B", "#5790fc")):
                values = np.array([p for c, p in pulls.items() if c.startswith(family + "_")])
                if values.size:
                    ax.hist(np.clip(values, -10, 10), bins=40, range=(-10, 10), histtype="step", color=colour,
                            label=f"family {family} (rms {values.std():.2f})")
            ax.set_xlabel(f"category residual / sigma ({quantity})")
            ax.set_title(f"{LABELS[tag]}: chi2/ndf {solution[tag][quantity]['chi2']:.0f}/{solution[tag][quantity]['ndf']}")
            ax.legend(fontsize="small")
        h4l_style.save(fig, plots / f"pulls_{tag}")
        if response:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            for ax, quantity in zip(axes, ("scale", "smear")):
                for family, colour in (("A", "#e42536"), ("B", "#5790fc")):
                    values = np.array([v[quantity] for c, v in response["response"][tag].items()
                                       if c.startswith(family + "_") and quantity in v])
                    if values.size:
                        ax.hist(np.clip(values, -1, 2.5), bins=70, range=(-1, 2.5), histtype="step", color=colour,
                                label=f"family {family}")
                for limit in response["limits"]:
                    ax.axvline(limit, color="grey", ls="--")
                ax.set_xlabel(f"category response ({quantity})")
                ax.set_title(LABELS[tag])
                ax.legend(fontsize="small")
            h4l_style.save(fig, plots / f"response_{tag}")
        counts = {}
        for family in FAMILIES:
            for n in sorted({0, last}):
                counts[f"{family}_iter{n}"] = gallery(load_fits(base / f"iter_{n:02d}" / "main" / "fits", family, tag),
                                                      f"{LABELS[tag]} family {family} template fits, iteration {n}",
                                                      plots / f"gallery_{family}_{tag}_iter{n:02d}.pdf")
            counts[f"kinematics_{family}"] = kinematics(final_dir / "main" / "hist", family, tag,
                                                        plots / f"kinematics_{family}_{tag}.pdf")
        report["flavours"][name]["plotted"] = counts
        numbers = {}
        for label, n in (("before", 0), ("after", last)):
            hist_dir = base / f"iter_{n:02d}" / "main" / "hist"
            data = summed(sorted(str(p) for p in hist_dir.glob("data_*.root")), f"mB_{tag}_")
            mc = summed(sorted(str(p) for p in hist_dir.glob("mc_*.root")), f"mB_{tag}_")
            numbers[label] = zpeak_plot(data, mc, f"{LABELS[tag]}, {'uncorrected' if n == 0 else 'calibrated'}", lumi,
                                        plots / f"zpeak_{tag}_{label}")
        report["flavours"][name]["zpeak_totals"] = numbers
    (plots / "report.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({name: {k: v for k, v in f.items() if k in ("scale_shift", "smear")}
                      for name, f in report["flavours"].items()}, indent=1))
    print(f"[plots] {plots}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
