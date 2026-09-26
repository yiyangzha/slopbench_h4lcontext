"""Tag-and-probe plots, fit galleries and the scale-factor report (stage 3c).

    pixi run py -- analysis_v3/tnp/scripts/plot_tnp.py --extract v2 --run v3 --fit-label fits_i1

Reads production_v3/tnp/<extract>/<run>/<fit-label>/sf.json and fits/fit_<f>_<step>.json
and writes <run>/<fit-label>/plots/:
  * eff_<f>_<step>.png/pdf: data and MC tag-and-probe efficiencies against
    the probe pT, one panel per |eta| bin, with the data/MC scale factor and
    its statistical and fit-model uncertainties (the MC generator truth is
    not shown: it never enters a scale factor);
  * sf_<f>_full.png/pdf: the full single-lepton SF (the product of the chain
    id | loose, sip | id, iso | id and sip) against the directly measured full
    step;
  * gallery/<f>_<step>_<role>_pNN.png: every fit, 8 per page, the pass
    (black) and the fail (red) spectrum side by side on a linear scale with
    the fitted total (solid) and background (dashed) and the pulls below,
    labelled with the bin ranges, the efficiency of the nominal and the
    alternative models and of the counting, the status and the flags
    (parameters at limits, relaxed width bound, shared resolution,
    override, binomial error floor, template sources) -- the pages that are
    inspected one by one;
  * gallery_alternatives/alt_<f>_<step>_<role>_pNN.png: every fit whose
    alternative model deviates from the nominal by more than max(3 sigma,
    0.003), nominal and alternative side by side (the alternatives entering
    the fit-model systematic are inspected too);
  * report.json: the full SF at pT = 45 GeV, |eta| = 1.2 (bilinear between
    bin centres) and the list of fits by status.
Data bin errors must be sqrt(N); the script stops otherwise.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import multiprocessing
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402

FLAVOURS = {"mm": "muon", "ee": "electron"}
STEPS = {"id": "ID | loose", "sip": "SIP | ID", "iso": "isolation | ID, SIP", "full": "ID, SIP, iso | loose"}
ETA_NAME = {"mm": "|eta|", "ee": "|eta_SC|"}
ETA_TOP = {"mm": 2.4, "ee": 2.5}


def bin_label(tag: str, bins: dict, b: str) -> str:
    if b == "all":
        return "inclusive"
    k = int(b)
    n_eta = len(bins["eta_edges"])
    ip, ie = divmod(k, n_eta)
    pt, eta = bins["pt_edges"], bins["eta_edges"]
    pt_hi = f"{pt[ip + 1]:g}" if ip + 1 < len(pt) else "inf"
    eta_hi = f"{eta[ie + 1]:g}" if ie + 1 < len(eta) else f"{ETA_TOP[tag]:g}"
    return f"pT {pt[ip]:g}-{pt_hi}, {ETA_NAME[tag]} {eta[ie]:g}-{eta_hi}"


def efficiency_plots(sf: dict, tag: str, name: str, step: str, plots: Path) -> None:
    flavour = sf["flavours"][name]
    bins = flavour["bins"]
    pt_edges, eta_edges = bins["pt_edges"], bins["eta_edges"]
    n_eta = len(eta_edges)
    rows = {row["bin"]: row for row in flavour["steps"][step]}
    upper = pt_edges[1:] + [pt_edges[-1] * 1.5]
    x = np.array([0.5 * (a + b) for a, b in zip(pt_edges, upper)])
    xerr = np.array([0.5 * (b - a) for a, b in zip(pt_edges, upper)])
    ncol = min(n_eta, 6)
    nrow = math.ceil(n_eta / ncol)
    fig = plt.figure(figsize=(3.4 * ncol, 5.2 * nrow))
    grid = fig.add_gridspec(2 * nrow, ncol, height_ratios=[2, 1] * nrow, hspace=0.35, wspace=0.35)
    for ie in range(n_eta):
        r, c = divmod(ie, ncol)
        top = fig.add_subplot(grid[2 * r, c])
        bottom = fig.add_subplot(grid[2 * r + 1, c], sharex=top)
        eta_hi = eta_edges[ie + 1] if ie + 1 < n_eta else ETA_TOP[tag]
        for role, colour, marker in (("data", "black", "o"), ("mc", "#5790fc", "s")):
            values, errors = [], []
            for ip in range(len(pt_edges)):
                entry = rows[str(ip * n_eta + ie)][role]["nominal"]
                values.append(entry["efficiency"] if entry["efficiency"] is not None else np.nan)
                errors.append(entry["error"] if entry["error"] is not None else 0.0)
            top.errorbar(x, values, errors, xerr, fmt=marker, ms=2.5, color=colour, lw=0.8, label=role if ie == 0 else None)
        top.set_title(f"{ETA_NAME[tag]} {eta_edges[ie]:g}-{eta_hi:g}", fontsize=8)
        top.set_xscale("log")
        top.tick_params(labelsize=7, labelbottom=False)
        if c == 0:
            top.set_ylabel("efficiency", fontsize=8)
        sfv, sfe, sfm = [], [], []
        for ip in range(len(pt_edges)):
            s = rows[str(ip * n_eta + ie)].get("sf")
            sfv.append(s["value"] if s else np.nan)
            sfe.append(s["stat"] if s else 0.0)
            sfm.append(math.hypot(s["stat"], s["fit_model"]) if s else 0.0)
        sfv, sfe, sfm = np.array(sfv), np.array(sfe), np.array(sfm)
        bottom.fill_between(x, sfv - sfm, sfv + sfm, color="#f89c20", alpha=0.3, label="stat (+) fit model" if ie == 0 else None)
        bottom.errorbar(x, sfv, sfe, xerr, fmt="o", ms=2.5, color="black", lw=0.8, label="data / MC (stat)" if ie == 0 else None)
        bottom.axhline(1, color="grey", lw=0.6)
        bottom.set_xlabel("probe pT [GeV]", fontsize=8)
        if c == 0:
            bottom.set_ylabel("SF", fontsize=8)
        bottom.tick_params(labelsize=7)
        if ie == 0:
            top.legend(fontsize=6)
            bottom.legend(fontsize=6)
    fig.suptitle(f"{name} {STEPS[step]} tag-and-probe efficiency", fontsize=11)
    h4l_style.save(fig, plots / f"eff_{tag}_{step}")
    plt.close(fig)


def full_sf_plot(sf: dict, tag: str, name: str, plots: Path) -> None:
    flavour = sf["flavours"][name]
    bins = flavour["bins"]
    pt_edges, eta_edges = bins["pt_edges"], bins["eta_edges"]
    n_eta = len(eta_edges)
    product = {row["bin"]: row for row in flavour["product"]}
    direct = {row["bin"]: row for row in flavour["steps"]["full"]}
    upper = pt_edges[1:] + [pt_edges[-1] * 1.5]
    x = np.array([0.5 * (a + b) for a, b in zip(pt_edges, upper)])
    xerr = np.array([0.5 * (b - a) for a, b in zip(pt_edges, upper)])
    ncol = min(n_eta, 6)
    nrow = math.ceil(n_eta / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.2 * nrow), squeeze=False)
    for ie in range(n_eta):
        ax = axes[divmod(ie, ncol)]
        pv, pe, dv, de = [], [], [], []
        for ip in range(len(pt_edges)):
            b = str(ip * n_eta + ie)
            p = product[b]
            pv.append(p["value"] if p.get("value") is not None else np.nan)
            pe.append(p.get("total", 0.0) if p.get("value") is not None else 0.0)
            d = direct[b].get("sf")
            dv.append(d["value"] if d else np.nan)
            de.append(math.hypot(d["stat"], d["fit_model"]) if d else 0.0)
        ax.errorbar(x, pv, pe, xerr, fmt="o", ms=2.5, color="black", lw=0.8, label="product of the chain")
        ax.errorbar(x * 1.04, dv, de, fmt="s", ms=2.5, color="#e42536", lw=0.8, label="full step (direct)")
        ax.axhline(1, color="grey", lw=0.6)
        eta_hi = eta_edges[ie + 1] if ie + 1 < n_eta else ETA_TOP[tag]
        ax.set_title(f"{ETA_NAME[tag]} {eta_edges[ie]:g}-{eta_hi:g}", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("probe pT [GeV]", fontsize=8)
        ax.set_ylabel("SF (ID, SIP, iso)", fontsize=8)
        ax.tick_params(labelsize=7)
        if ie == 0:
            ax.legend(fontsize=6)
    for k in range(n_eta, nrow * ncol):
        axes[divmod(k, ncol)].axis("off")
    fig.suptitle(f"{name} full single-lepton selection scale factor", fontsize=11)
    fig.tight_layout()
    h4l_style.save(fig, plots / f"sf_{tag}_full")
    plt.close(fig)


def gallery(fits: dict, bins: dict, tag: str, role: str, out_dir: Path) -> int:
    """Eight fits per page; per fit the pass and the fail spectrum side by side (linear scale) with pulls."""
    entries = [e for e in fits["results"] if e["role"] == role and "fits" in e and e["fits"].get("nominal", {}).get("curve")]
    out_dir.mkdir(parents=True, exist_ok=True)
    per_page = 8
    for page, start in enumerate(range(0, len(entries), per_page)):
        fig = plt.figure(figsize=(18, 20))
        outer = fig.add_gridspec(4, 2, top=0.95, bottom=0.03, left=0.04, right=0.99, hspace=0.38, wspace=0.12)
        for slot, entry in enumerate(entries[start:start + per_page]):
            r, c = divmod(slot, 2)
            cell = outer[r, c].subgridspec(2, 2, height_ratios=[3, 1], hspace=0.06, wspace=0.16)
            fit = entry["fits"]["nominal"]
            curve = fit["curve"]
            x = np.array(curve["x"])
            flags = []
            if fit.get("at_limit"):
                flags.append("limit:" + ",".join(fit["at_limit"]))
            if fit.get("width_bound_relaxed"):
                flags.append("relaxed")
            if fit.get("shared_resolution"):
                flags.append("shared-res")
            if fit.get("override"):
                flags.append("override")
            if fit.get("error_floor_applied"):
                flags.append("err=binomial-floor")
            if fit.get("templates"):
                sources = fit["templates"]
                flags.append(f"tmpl {sources['pass']['source']}/{sources['fail']['source']}")
            err = fit["efficiency_error"]
            err_text = f"{err:.4f}" if isinstance(err, (int, float)) and err == err else "nan"
            alts = entry["fits"]
            alt_text = (f"alt_sig {alts['alt_signal']['efficiency']:.4f}, alt_bkg {alts['alt_background']['efficiency']:.4f}, "
                        f"count {entry['counting']['efficiency']:.4f}")
            for side, kind, colour in ((0, "pass", "black"), (1, "fail", "#e42536")):
                top = fig.add_subplot(cell[0, side])
                bottom = fig.add_subplot(cell[1, side], sharex=top)
                y, ey = np.array(curve[kind]), np.array(curve[f"{kind}_error"])
                if role == "data" and not np.allclose(ey ** 2, y, rtol=1e-6, atol=1e-6):
                    raise SystemExit(f"data bin errors are not sqrt(N) in {tag} {entry['step']} bin {entry['bin']} ({kind})")
                model, bkg = np.array(curve[f"model_{kind}"]), np.array(curve[f"background_{kind}"])
                top.errorbar(x, y, ey, fmt="o", ms=1.8, color=colour, lw=0.6)
                top.plot(x, model, color="#3f90da", lw=1.1)
                top.plot(x, bkg, color="#3f90da", lw=0.8, ls="--")
                top.set_ylim(0, 1.15 * max(float(np.max(y + ey)) if len(y) else 1.0, float(np.max(model)) if len(model) else 1.0, 1.0))
                pull = np.divide(y - model, ey, out=np.zeros_like(y), where=ey > 0)
                bottom.plot(x, pull, ".", ms=2.5, color=colour)
                bottom.axhline(0, color="grey", lw=0.5)
                bottom.set_ylim(-5, 5)
                top.tick_params(labelsize=6, labelbottom=False)
                bottom.tick_params(labelsize=6)
                raw = entry.get(f"raw_{kind}", float("nan"))
                if side == 0:
                    top.set_title(f"{entry['step']} bin {entry['bin']}: {bin_label(tag, bins, entry['bin'])} [{entry['status']}]\n"
                                  f"eps {fit['efficiency']:.4f} +- {err_text}; {alt_text}\npass: {raw:.0f} entries, "
                                  f"{fit['bin_width']:g} GeV bins", fontsize=6.5, loc="left")
                else:
                    top.set_title(f"{' '.join(flags)}\n\nfail: {raw:.0f} entries", fontsize=6.5, loc="left")
        fig.suptitle(f"{FLAVOURS[tag]} {entries[start]['step']} {role}: pass (black) | fail (red) per fit; "
                     f"model solid, background dashed; x = m_ll [GeV], pulls below", fontsize=10, y=0.99)
        fig.savefig(out_dir / f"{tag}_{entries[start]['step']}_{role}_p{page:02d}.png", dpi=100)
        plt.close(fig)
    return len(entries)


def draw_fit_panels(fig, cell_row, column0: int, entry: dict, fit: dict, role: str, tag: str, label: str) -> None:
    """Pass and fail panels (with pulls) of one model of one fit into columns column0, column0 + 1 of cell_row."""
    curve = fit["curve"]
    x = np.array(curve["x"])
    for side, kind, colour in ((0, "pass", "black"), (1, "fail", "#e42536")):
        top = fig.add_subplot(cell_row[0, column0 + side])
        bottom = fig.add_subplot(cell_row[1, column0 + side], sharex=top)
        y, ey = np.array(curve[kind]), np.array(curve[f"{kind}_error"])
        model, bkg = np.array(curve[f"model_{kind}"]), np.array(curve[f"background_{kind}"])
        top.errorbar(x, y, ey, fmt="o", ms=1.6, color=colour, lw=0.5)
        top.plot(x, model, color="#3f90da", lw=1.1)
        top.plot(x, bkg, color="#3f90da", lw=0.8, ls="--")
        top.set_ylim(0, 1.15 * max(float(np.max(y + ey)) if len(y) else 1.0, float(np.max(model)) if len(model) else 1.0, 1.0))
        pull = np.divide(y - model, ey, out=np.zeros_like(y), where=ey > 0)
        bottom.plot(x, pull, ".", ms=2.2, color=colour)
        bottom.axhline(0, color="grey", lw=0.5)
        bottom.set_ylim(-5, 5)
        top.tick_params(labelsize=6, labelbottom=False)
        bottom.tick_params(labelsize=6)
        if side == 0:
            top.set_title(f"{label}: eps {fit['efficiency']:.4f} +- {fit['efficiency_error']:.4f}"
                          f"{' [unusable]' if not fit.get('usable') else ''}", fontsize=6.5, loc="left")


def alternatives_gallery(fits: dict, bins: dict, tag: str, role: str, out_dir: Path) -> int:
    """Four fits per page, each row nominal (pass | fail) and the deviating alternative (pass | fail): every fit
    whose alternative efficiency differs from the nominal by more than max(3 sigma, 0.003) -- the alternatives
    that enter the fit-model systematic are inspected like the nominal fits."""
    rows = []
    for entry in fits["results"]:
        if entry["role"] != role or "fits" not in entry:
            continue
        nominal = entry["fits"]["nominal"]
        for model in ("alt_signal", "alt_background"):
            alt = entry["fits"][model]
            if not alt.get("curve") or not nominal.get("curve"):
                continue
            # A nominal fit without a valid error (null in the report) is compared with the fixed 0.003 only.
            error = nominal["efficiency_error"] if isinstance(nominal["efficiency_error"], (int, float)) else float("nan")
            threshold = max(3.0 * error, 0.003) if error == error else 0.003
            if abs(alt["efficiency"] - nominal["efficiency"]) > threshold:
                rows.append((entry, model))
    out_dir.mkdir(parents=True, exist_ok=True)
    per_page = 4
    for page, start in enumerate(range(0, len(rows), per_page)):
        fig = plt.figure(figsize=(18, 20))
        outer = fig.add_gridspec(4, 1, top=0.95, bottom=0.03, left=0.04, right=0.99, hspace=0.40)
        for slot, (entry, model) in enumerate(rows[start:start + per_page]):
            cell = outer[slot, 0].subgridspec(2, 4, height_ratios=[3, 1], hspace=0.06, wspace=0.16)
            head = f"{entry['step']} bin {entry['bin']}: {bin_label(tag, bins, entry['bin'])}"
            draw_fit_panels(fig, cell, 0, entry, entry["fits"]["nominal"], role, tag, f"{head} | nominal")
            draw_fit_panels(fig, cell, 2, entry, entry["fits"][model], role, tag, model)
        fig.suptitle(f"{FLAVOURS[tag]} {rows[start][0]['step']} {role}: nominal (left) and deviating alternative (right); "
                     f"pass black, fail red; model solid, background dashed", fontsize=10, y=0.99)
        fig.savefig(out_dir / f"alt_{tag}_{rows[start][0]['step']}_{role}_p{page:02d}.png", dpi=100)
        plt.close(fig)
    return len(rows)


def galleries(task: tuple) -> str:
    """The fit and alternative-model galleries of one flavour, step and role (one worker)."""
    fits_path, bins, tag, role, plots = task
    fits = json.loads(Path(fits_path).read_text(encoding="utf-8"))
    gallery(fits, bins, tag, role, plots / "gallery")
    alternatives_gallery(fits, bins, tag, role, plots / "gallery_alternatives")
    return f"{tag}_{fits['step']}_{role}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--extract", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--fit-label", default="fits")
    parser.add_argument("--no-gallery", action="store_true")
    parser.add_argument("--workers", type=int, default=16, help="parallel gallery workers (one flavour, step and role each)")
    args = parser.parse_args()
    base = PRODUCTION / "tnp" / args.extract / args.run / args.fit_label
    sf = json.loads((base / "sf.json").read_text(encoding="utf-8"))
    plots = base / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    report = {"schema": "h4l_v3_tnp_report/2", "run": str(base), "flavours": {}, "status": {}}
    gallery_tasks = []
    for tag, name in FLAVOURS.items():
        bins = sf["flavours"][name]["bins"]
        for step in STEPS:
            efficiency_plots(sf, tag, name, step, plots)
            fits_path = base / "fits" / f"fit_{tag}_{step}.json"
            fits = json.loads(fits_path.read_text(encoding="utf-8"))
            for entry in fits["results"]:
                report["status"].setdefault(f"{tag}_{step}_{entry['role']}", {}).setdefault(entry["status"], []).append(entry["bin"])
            if not args.no_gallery:
                gallery_tasks += [(str(fits_path), bins, tag, role, plots) for role in ("data", "mc")]
        full_sf_plot(sf, tag, name, plots)
        rp = sf["flavours"][name]["report_point"]
        inclusive = next(r for r in sf["flavours"][name]["product"] if r["bin"] == "all")
        report["flavours"][name] = {"sf_full_at_45GeV_eta1p2": rp, "sf_full_inclusive": inclusive}
        value = rp["sf_full"].get("value")
        print(f"[tnp] {name}: SF(full) at 45 GeV, |eta| 1.2: {value} +- {rp['sf_full'].get('error')}; inclusive {inclusive.get('value')}")
    if gallery_tasks:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("fork")) as pool:
            for done in pool.map(galleries, gallery_tasks):
                print(f"[tnp] galleries {done}", flush=True)
    (plots / "report.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    counts = {k: {s: len(v) for s, v in d.items()} for k, d in report["status"].items()}
    print(f"[tnp] fit status counts: {json.dumps(counts)}")
    print(f"[tnp] plots {plots}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
