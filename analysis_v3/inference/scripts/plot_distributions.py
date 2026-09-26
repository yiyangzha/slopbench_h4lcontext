"""Yields and distributions of the selected four-lepton events (JHEP 11 (2017) 047 sections 5-7; stage 8).

    pixi run py -- analysis_v3/inference/scripts/plot_distributions.py --select v5 --zx zx_v1 \
        --sf production_v3/tnp/v2/<run>/<label>/sf.json --categories analysis_v3/reconstruction/config/categories_ul16_v2.json \
        --label dist_v1

Expected contributions (pre-fit, the SM at m_H = 125 GeV): signal (ggH, VBF, VH MC) and qqZZ, ggZZ MC with
the genWeight normalization and every lepton weighted by its tag-and-probe scale factor; Z+X from the combined
OS/SS estimate per final state (m4l > 70 GeV) with its fitted m4l shape.  Written:
  * Table 1: per final state, m4l > 70 GeV (and 118-130 GeV): signal, qqZZ, ggZZ, Z+X, total, observed;
  * Table 2: per category in 118-130 GeV: signal per production mode, backgrounds, total, observed;
  * m4l distributions: 70-870 GeV (10 GeV bins) and 70-170 GeV (2 GeV bins), all and per final state, data
    against the stacked expectation; 105-140 GeV per category;
  * m_Z2 against m_Z1 in 118-130 GeV (data points over the signal MC density);
  * D_bkg^kin against m4l (data with their per-event mass uncertainty, over the qqZZ MC density);
  * the production discriminants D_2jet, D_1jet, D_WH, D_ZH (categories configuration) in 118-130 GeV.
Writes production_v3/results/<select>/<label>/ (yields.json, yields.md, plots).
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
import h4l_categories as cats  # noqa: E402
import h4l_select_io as io  # noqa: E402
import h4l_sf  # noqa: E402
import h4l_shapes as shapes  # noqa: E402

FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}
SIGNAL = {"GluGluToHToZZ_M125": "ggH", "VBF_HToZZ_M125": "VBF", "VHToZZ_M125": "VH"}
BACKGROUNDS = {"qqZZ": ["ZZTo4L"], "ggZZ": ["GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"]}
BRANCHES = ["final_state", "m4l", "m4l_refit", "m4l_err", "m4l_refit_err", "mz1", "mz2", "d_bkg_kin", "l_pdg", "l_pt", "l_eta",
            "l_eta_sc"] + cats.CATEGORY_BRANCHES
COLOURS = {"signal": "#e42536", "qqZZ": "#5790fc", "ggZZ": "#3f90da", "zx": "#7a21dd"}


def zx_expectation(zx: dict, fs: str, lo: float, hi: float) -> float:
    c = zx["combined"][fs]
    shape = zx["shapes"][fs]
    grid = np.linspace(70.0, 870.0, 16001)
    dens = shapes.zx_pdf(grid, shape, 70.0, 870.0)
    inside = (grid > lo) & (grid < hi)
    return float(c["value"] * np.trapezoid(dens[inside], grid[inside]))


def zx_hist(zx: dict, fs: str, edges: np.ndarray) -> np.ndarray:
    return np.array([zx_expectation(zx, fs, a, b) for a, b in zip(edges[:-1], edges[1:])])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--zx", required=True)
    parser.add_argument("--sf", type=Path, required=True)
    parser.add_argument("--categories", type=Path, default=REPO / "analysis_v3/reconstruction/config/categories_ul16_v2.json")
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    out_dir = PRODUCTION / "results" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    scan = io.load_scan(args.select)
    lumi = io.lumi_fb(scan)
    zx = json.loads((PRODUCTION / "backgrounds" / args.select / args.zx / "zx.json").read_text(encoding="utf-8"))
    categories = json.loads((args.categories if args.categories.is_absolute() else REPO / args.categories).read_text(encoding="utf-8"))
    order = categories["order"]
    sf = h4l_sf.LeptonSF(args.sf if args.sf.is_absolute() else REPO / args.sf)
    above70 = lambda a: a["m4l"] > 70.0  # noqa: E731
    data = io.read(scan, ["pseudo_data"], "SR", BRANCHES, cut=above70)
    mc = {}
    for sample, name in SIGNAL.items():
        mc[name] = io.read(scan, [sample], "SR", BRANCHES, cut=above70)
    for name, samples in BACKGROUNDS.items():
        mc[name] = io.read(scan, samples, "SR", BRANCHES, cut=above70)
    for rows in list(mc.values()) + [data]:
        rows["category"] = cats.category_index(rows, categories)
    for rows in mc.values():
        rows["weight_sf"] = rows["w"] * sf.event(rows)["sf"]
    out_dir.mkdir(parents=True)
    plots = out_dir / "plots"
    plots.mkdir()

    # ---- Tables.
    def expected(name, sel_fn):
        r = mc[name]
        return float(np.sum(r["weight_sf"][sel_fn(r)]))

    table1 = {}
    for window_name, (lo, hi) in (("m4l_gt_70", (70.0, np.inf)), ("m4l_118_130", (118.0, 130.0))):
        table1[window_name] = {}
        for fs, name in list(FINAL_STATES.items()) + [(None, "4l")]:
            def sel(r, fs=fs, lo=lo, hi=hi):
                m = (r["m4l"] > lo) & (r["m4l"] < hi)
                return m if fs is None else m & (r["final_state"] == fs)
            entry = {mode: expected(mode, sel) for mode in SIGNAL.values()}
            entry["signal"] = sum(entry[m] for m in SIGNAL.values())
            for b in BACKGROUNDS:
                entry[b] = expected(b, sel)
            zx_fs = [name] if fs is not None else list(FINAL_STATES.values())
            entry["zx"] = sum(zx_expectation(zx, f, lo, min(hi, 870.0)) for f in zx_fs)
            entry["background"] = entry["qqZZ"] + entry["ggZZ"] + entry["zx"]
            entry["total"] = entry["signal"] + entry["background"]
            entry["observed"] = int(np.sum(sel(data)))
            table1[window_name][name] = entry
    table2 = {}
    for c, cat in enumerate(order):
        def sel(r, c=c):
            return (r["m4l_refit"] > 118.0) & (r["m4l_refit"] < 130.0) & (r["category"] == c)
        entry = {mode: expected(mode, sel) for mode in SIGNAL.values()}
        entry["signal"] = sum(entry[m] for m in SIGNAL.values())
        for b in BACKGROUNDS:
            entry[b] = expected(b, sel)
        entry["observed"] = int(np.sum(sel(data)))
        table2[cat] = entry
    # Z+X per category: its category fractions follow the OS-method control rows in build_model; the table gives the
    # final-state totals scaled by the categories' share of the qqZZ-like expectation only as a display proxy.
    tot_bkg = sum(table2[c]["qqZZ"] for c in order)
    zx_window = table1["m4l_118_130"]["4l"]["zx"]
    for c in order:
        table2[c]["zx_proxy"] = zx_window * table2[c]["qqZZ"] / tot_bkg if tot_bkg > 0 else 0.0
    report = {"schema": "h4l_v3_yields/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
              "luminosity_fb": lumi, "sf_payload": str(args.sf), "zx": str(PRODUCTION / "backgrounds" / args.select / args.zx / "zx.json"),
              "categories": str(args.categories), "table1": table1, "table2": table2,
              "note": "pre-fit SM expectation at m_H = 125 GeV; Z+X per category in table2 is a display proxy (qqZZ shares)"}
    (out_dir / "yields.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    lines = [f"# Yields ({lumi} fb^-1, pre-fit SM at m_H = 125 GeV)", ""]
    for window_name, t in table1.items():
        lines += [f"## Table 1, {window_name.replace('_', ' ')}", "", "| | 4mu | 4e | 2e2mu | 4l |", "|---|---|---|---|---|"]
        for key in ("ggH", "VBF", "VH", "signal", "qqZZ", "ggZZ", "zx", "background", "total", "observed"):
            lines.append(f"| {key} | " + " | ".join(f"{t[n][key]:.2f}" if key != "observed" else str(t[n][key])
                                                  for n in ("4mu", "4e", "2e2mu", "4l")) + " |")
        lines.append("")
    lines += ["## Table 2, 118-130 GeV per category", "", "| category | ggH | VBF | VH | signal | qqZZ | ggZZ | Z+X (proxy) | observed |",
              "|---|---|---|---|---|---|---|---|---|"]
    for cat in order:
        t = table2[cat]
        lines.append(f"| {cat} | {t['ggH']:.2f} | {t['VBF']:.2f} | {t['VH']:.2f} | {t['signal']:.2f} | {t['qqZZ']:.2f} | {t['ggZZ']:.2f} | "
                     f"{t['zx_proxy']:.2f} | {t['observed']} |")
    (out_dir / "yields.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- m4l distributions.
    def stack_plot(edges, fs_list, title, stem, variable="m4l", cat=None, log=False):
        centres = 0.5 * (edges[1:] + edges[:-1])
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        bottoms = np.zeros(len(centres))
        comps = []
        zx_h = np.zeros(len(centres))
        for fs in fs_list:
            zx_h += zx_hist(zx, FINAL_STATES[fs], edges)
        if cat is not None:
            # Display proxy: the category's share of the qqZZ expectation in the window (the likelihood uses the
            # category fractions of the OS-method control rows).
            q = mc["qqZZ"]
            window = (q["m4l_refit"] > edges[0]) & (q["m4l_refit"] < edges[-1])
            total = float(np.sum(q["weight_sf"][window]))
            zx_h *= float(np.sum(q["weight_sf"][window & (q["category"] == cat)])) / total if total > 0 else 0.0
        comps.append(("Z+X", zx_h, COLOURS["zx"]))
        for name, label in (("ggZZ", "gg->ZZ"), ("qqZZ", "qq->ZZ")):
            r = mc[name]
            sel = np.isin(r["final_state"], fs_list) & (True if cat is None else (r["category"] == cat))
            h, _ = np.histogram(r[variable][sel], bins=edges, weights=r["weight_sf"][sel])
            comps.append((label, h, COLOURS[name]))
        sig = np.zeros(len(centres))
        for mode in SIGNAL.values():
            r = mc[mode]
            sel = np.isin(r["final_state"], fs_list) & (True if cat is None else (r["category"] == cat))
            h, _ = np.histogram(r[variable][sel], bins=edges, weights=r["weight_sf"][sel])
            sig += h
        comps.append(("H(125)", sig, COLOURS["signal"]))
        for label, h, colour in comps:
            ax.bar(centres, h, width=np.diff(edges), bottom=bottoms, color=colour, alpha=0.75, label=label, align="center")
            bottoms = bottoms + h
        sel = np.isin(data["final_state"], fs_list) & (True if cat is None else (data["category"] == cat))
        hd, _ = np.histogram(data[variable][sel], bins=edges)
        ax.errorbar(centres, hd, yerr=np.sqrt(hd), fmt="ko", ms=3, label=f"data ({int(hd.sum())})")
        ax.set_xlabel(f"{variable} [GeV]")
        ax.set_ylabel(f"events / {edges[1] - edges[0]:g} GeV")
        if log:
            ax.set_yscale("log")
            ax.set_ylim(0.05, None)
        ax.set_title(f"{title}, {lumi} fb-1 (pre-fit SM, m_H = 125 GeV)", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(plots / f"{stem}.png", dpi=120)
        fig.savefig(plots / f"{stem}.pdf")
        plt.close(fig)
        return float(hd.sum()), float(bottoms.sum())

    checks = {}
    all_fs = list(FINAL_STATES)
    checks["m4l_full"] = stack_plot(np.arange(70.0, 870.0 + 1e-9, 10.0), all_fs, "4l, 70-870 GeV", "m4l_full", log=True)
    checks["m4l_low"] = stack_plot(np.arange(70.0, 170.0 + 1e-9, 2.0), all_fs, "4l, 70-170 GeV", "m4l_low")
    for fs, name in FINAL_STATES.items():
        checks[f"m4l_low_{name}"] = stack_plot(np.arange(70.0, 170.0 + 1e-9, 2.0), [fs], f"{name}, 70-170 GeV", f"m4l_low_{name}")
    for c, cat in enumerate(order):
        checks[f"m4l_{cat}"] = stack_plot(np.arange(105.0, 140.0 + 1e-9, 2.5), all_fs, f"{cat}, 105-140 GeV", f"m4l_refit_{cat}",
                                          variable="m4l_refit", cat=c)

    # ---- m_Z2 vs m_Z1 in 118-130 GeV.
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    sig = {k: np.concatenate([mc[m][k] for m in SIGNAL.values()]) for k in ("mz1", "mz2", "m4l", "weight_sf")}
    s_sel = (sig["m4l"] > 118) & (sig["m4l"] < 130)
    ax.hist2d(sig["mz1"][s_sel], sig["mz2"][s_sel], bins=[np.linspace(40, 120, 41), np.linspace(10, 120, 45)],
              weights=sig["weight_sf"][s_sel], cmap="Reds")
    d_sel = (data["m4l"] > 118) & (data["m4l"] < 130)
    ax.plot(data["mz1"][d_sel], data["mz2"][d_sel], "ko", ms=3, label=f"data ({int(d_sel.sum())})")
    ax.set_xlabel("m_Z1 [GeV]")
    ax.set_ylabel("m_Z2 [GeV]")
    ax.set_title("118 < m4l < 130 GeV: data over the signal MC", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plots / "mz2_vs_mz1.png", dpi=120)
    fig.savefig(plots / "mz2_vs_mz1.pdf")
    plt.close(fig)

    # ---- D_bkg^kin vs m4l with per-event errors.
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    q = mc["qqZZ"]
    q_sel = (q["m4l"] > 70) & (q["m4l"] < 170)
    ax.hist2d(q["m4l"][q_sel], q["d_bkg_kin"][q_sel], bins=[np.linspace(70, 170, 51), np.linspace(0, 1, 21)], weights=q["weight_sf"][q_sel],
              cmap="Blues")
    d_sel = (data["m4l"] > 70) & (data["m4l"] < 170)
    ax.errorbar(data["m4l"][d_sel], data["d_bkg_kin"][d_sel], xerr=data["m4l_err"][d_sel], fmt="ko", ms=2.5, lw=0.6,
                label="data (per-event m4l uncertainty)")
    ax.set_xlabel("m4l [GeV]")
    ax.set_ylabel("D_bkg^kin")
    ax.set_title("D_bkg^kin against m4l; qqZZ MC density", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plots / "dkin_vs_m4l.png", dpi=120)
    fig.savefig(plots / "dkin_vs_m4l.pdf")
    plt.close(fig)

    # ---- Production discriminants in 118-130 GeV.
    for key, label in (("d_2jet", "D_2jet"), ("d_1jet", "D_1jet"), ("d_wh_mela", "D_WH"), ("d_zh_mela", "D_ZH")):
        fig, ax = plt.subplots(figsize=(5.6, 4.2))
        edges = np.linspace(0, 1, 21)
        for name, colour in (("ggH", "#e42536"), ("VBF", "#f89c20"), ("VH", "#964a8b"), ("qqZZ", "#5790fc")):
            r = mc[name]
            disc = cats.production_discriminants(r, categories)[key] if "discriminants" in categories else r[key]
            sel = (r["m4l"] > 118) & (r["m4l"] < 130) & (disc >= 0)
            h, _ = np.histogram(disc[sel], bins=edges, weights=r["weight_sf"][sel])
            ax.step(0.5 * (edges[1:] + edges[:-1]), h / max(h.sum(), 1e-12), where="mid", color=colour, label=name)
        ax.set_xlabel(label + " (categories configuration)")
        ax.set_ylabel("normalized")
        ax.set_title(f"{label}, 118 < m4l < 130 GeV (MC, unit area)", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(plots / f"disc_{label}.png", dpi=120)
        fig.savefig(plots / f"disc_{label}.pdf")
        plt.close(fig)

    report["plot_checks"] = {k: {"data_events": v[0], "expected": v[1]} for k, v in checks.items()}
    (out_dir / "yields.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    t = table1["m4l_gt_70"]["4l"]
    print(f"[dist] m4l > 70: signal {t['signal']:.2f}, qqZZ {t['qqZZ']:.2f}, ggZZ {t['ggZZ']:.2f}, Z+X {t['zx']:.2f}, total {t['total']:.2f}, "
          f"observed {t['observed']}\n[dist] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
