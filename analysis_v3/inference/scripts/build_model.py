"""Inputs of the H -> 4l likelihood (stage 7): events, yields, shapes and templates per channel.

    pixi run py -- analysis_v3/inference/scripts/build_model.py --select v5 --signal sm_v1 --zx zx_v1 \
        --sf production_v3/tnp/v2/<run>/<label>/sf.json --label model_v1

Channels (--channels): final state (4mu, 4e, 2e2mu) x category (the signal model's categories;
default), x nothing ("inclusive": the fiducial measurement without categorization), or x the
reconstruction-level bin of a differential observable ("pt4l", "njets", "ptj1"; bins of
h4l_fiducial.DIFFERENTIAL); the window [lo, hi] of the signal model (105-140 GeV).  Written per
channel:
  * the data events: m4l_refit, its per-event error, m4l, its error, D_bkg^kin;
  * expected yields (MC normalized with the genWeight convention, every lepton weighted by
    its tag-and-probe scale factor): signal per production mode at m_H = 125 GeV (the
    m_H dependence from the signal model), qqZZ (ZZTo4L), ggZZ (GGZZ*); Z+X from the
    combined OS/SS estimate of the final state times its window fraction (zx.json)
    times the category fraction of the OS-method control rows (signed weights; a
    negative fraction is set to zero and the fractions renormalized to preserve the
    final state's total, recorded);
  * the scale-factor variations of every MC yield (muon and electron SFs moved by their
    total errors) and eff_correction (the yield with SFs over the yield without);
  * the MC statistical uncertainty of every MC yield;
  * per signal mode the generator-level composition of its yield ("breakdown": production class
    ggH / VBF / VH_had / VH_lep; stage-0 class x |y_H| < 2.5 or above; fiducial final state or
    non-fiducial; fiducial bin of every differential observable or non-fiducial), from which the
    likelihood builds the per-mode, stage-0 and fiducial parameterizations.
Shapes per final state (shared by the categories):
  * signal: the DCB of the signal model (refit and non-refit, fixed and per-event width),
    the VH non-resonant Landau;
  * qqZZ and ggZZ: a cubic Bernstein polynomial fitted (weighted, unbinned) to the MC in
    the window;
  * Z+X: the Landau-plus-exponential shape of the Z+X estimate (its window fraction from the predicted rows);
  * D_bkg^kin: templates P(D bin | m4l bin) in 5 m4l bins x 5 D bins (the quintiles of the
    signal D in the window, all final states), from the signal MC, the qqZZ and ggZZ MC
    and, for Z+X, the same-sign control rows (weighted by f3 f4; their MELA has one
    charge flipped); a small floor keeps every bin positive;
  * per-event relative mass error: templates P(dm/m bin) in 5 bins (the quintiles of the
    signal of the final state) for the 3D fit, same sources.
Writes production_v3/inference/<select>/<label>/model.json and events.npz.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
sys.path.insert(0, str(REPO / "analysis_v3/backgrounds/scripts"))
import h4l_categories as cats  # noqa: E402
import h4l_fiducial as fid  # noqa: E402
import h4l_select_io as io  # noqa: E402
import h4l_sf  # noqa: E402
import h4l_shapes as shapes  # noqa: E402
from iminuit import Minuit  # noqa: E402

FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}
SIGNAL = {"GluGluToHToZZ_M125": "ggH", "VBF_HToZZ_M125": "VBF", "VHToZZ_M125": "VH"}
BACKGROUNDS = {"qqZZ": ["ZZTo4L"], "ggZZ": ["GGZZ4Mu", "GGZZ4E", "GGZZ2E2Mu"]}
EVENT_BRANCHES = ["final_state", "m4l", "m4l_err", "m4l_refit", "m4l_refit_err", "d_bkg_kin", "l_pdg", "l_pt", "l_eta", "l_eta_sc",
                  "mz1", "mz2", "pt4l", "jet_pt", "jet_eta"] + cats.CATEGORY_BRANCHES
SIGNAL_TRUTH_BRANCHES = ["g_vh_class", "g_h_y", "g_fid_pass", "g_fid_final_state", "g_fid_pt4l", "g_fid_njets", "g_fid_jet1_pt"]
CHANNEL_KINDS = ["categories", "inclusive"] + list(fid.DIFFERENTIAL)
M_BINS = 5
D_BINS = 5
E_BINS = 5


def fit_bernstein(x, w, lo, hi, order=3) -> list[float]:
    sel = (x > lo) & (x < hi)
    xx, ww = x[sel], w[sel]
    norm = np.sum(ww)

    def nll(*c):
        coeff = np.concatenate([[1.0], np.abs(c)])
        dens = shapes.bernstein_pdf(xx, coeff, lo, hi)
        return -np.sum(ww * np.log(np.maximum(dens, 1e-300))) / norm * len(xx)

    m = Minuit(nll, *([1.0] * order))
    m.errordef = Minuit.LIKELIHOOD
    m.migrad()
    return [1.0] + [float(abs(v)) for v in m.values]


def template_2d(m, d, w, m_edges, d_edges, floor=1e-3) -> list[list[float]]:
    """P(D bin | m4l bin), each m4l row normalized to 1 (a floor keeps every bin positive)."""
    h, _, _ = np.histogram2d(m, d, bins=[m_edges, d_edges], weights=w)
    h = np.maximum(h, 0.0)
    out = []
    for row in h:
        total = row.sum()
        p = row / total if total > 0 else np.full(len(row), 1.0 / len(row))
        p = np.maximum(p, floor)
        out.append((p / p.sum()).tolist())
    return out


def template_1d(x, w, edges, floor=1e-3) -> list[float]:
    h, _ = np.histogram(x, bins=edges, weights=w)
    h = np.maximum(h, 0.0)
    p = h / h.sum() if h.sum() > 0 else np.full(len(h), 1.0 / len(h))
    p = np.maximum(p, floor)
    return (p / p.sum()).tolist()


def breakdown(mode: str, rows: dict, sel: np.ndarray, w: np.ndarray) -> dict:
    """Generator-level composition (fractions of the weighted yield) of the selected signal rows of one mode."""
    total = float(np.sum(w))
    out = {}
    if total == 0:
        return out
    prod = fid.production_class(mode, rows["g_vh_class"][sel])
    central = np.abs(rows["g_h_y"][sel]) < fid.STXS_ABS_Y
    in_fid = rows["g_fid_pass"][sel].astype(bool)

    def fractions(keys):
        fr = {}
        for k in np.unique(keys):
            fr[str(k)] = float(np.sum(w[keys == k]) / total)
        return fr
    out["prod"] = fractions(prod)
    out["stxs"] = fractions(np.array([f"{c}_{'central' if y else 'forward'}" for c, y in zip(prod, central)], dtype=object))
    gen_fs = np.array([f"fid_{fid.FINAL_STATES.get(int(f), 'none')}" if ok else "nonfid"
                       for f, ok in zip(rows["g_fid_final_state"][sel], in_fid)], dtype=object)
    out["fid_fs"] = fractions(gen_fs)
    for obs, edges in fid.DIFFERENTIAL.items():
        idx = fid.bin_of(fid.gen_observable(obs, {k: v[sel] for k, v in rows.items() if k.startswith("g_")}, prefix="g_"), edges)
        out[f"fid_{obs}"] = fractions(np.array([f"bin{k}" if ok else "nonfid" for k, ok in zip(idx, in_fid)], dtype=object))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--select", required=True)
    parser.add_argument("--signal", required=True, help="signal-model label (production_v3/signal_model/<select>/<label>)")
    parser.add_argument("--zx", required=True, help="Z+X label (production_v3/backgrounds/<select>/<label>)")
    parser.add_argument("--sf", type=Path, required=True, help="tag-and-probe sf.json")
    parser.add_argument("--label", required=True)
    parser.add_argument("--channels", default="categories", choices=CHANNEL_KINDS,
                        help="final state x category (default), x nothing, or x reconstruction-level bin of an observable")
    args = parser.parse_args()
    out_dir = PRODUCTION / "inference" / args.select / args.label
    if out_dir.exists():
        raise SystemExit(f"{out_dir} exists; choose a new --label")
    scan = io.load_scan(args.select)
    signal_model = json.loads((PRODUCTION / "signal_model" / args.select / args.signal / "signal_model.json").read_text(encoding="utf-8"))
    zx = json.loads((PRODUCTION / "backgrounds" / args.select / args.zx / "zx.json").read_text(encoding="utf-8"))
    categories = signal_model["categories"]
    lo, hi = signal_model["window"]
    sf = h4l_sf.LeptonSF(args.sf if args.sf.is_absolute() else REPO / args.sf)
    if args.channels == "categories":
        order = categories["order"]
    elif args.channels == "inclusive":
        order = ["Inclusive"]
    else:
        order = [f"{args.channels}_{k}" for k in range(len(fid.DIFFERENTIAL[args.channels]) - 1)]
    near = lambda a: (a["m4l_refit"] > lo - 15.0) & (a["m4l_refit"] < hi + 15.0)  # noqa: E731

    def channel_bin(rows):
        """The channel index within the final state (the "category" key of the rows)."""
        if args.channels == "categories":
            return cats.category_index(rows, categories)
        if args.channels == "inclusive":
            return np.zeros(len(rows["m4l"]), dtype=np.int8)
        return fid.bin_of(fid.reco_observable(args.channels, rows), fid.DIFFERENTIAL[args.channels]).astype(np.int8)

    def load(samples, truth=False):
        rows = io.read(scan, samples, "SR", EVENT_BRANCHES + (SIGNAL_TRUTH_BRANCHES if truth else []), cut=near)
        rows["category"] = channel_bin(rows)
        return rows

    data = load(["pseudo_data"])
    mc = {name: load([sample], truth=True) for sample, name in SIGNAL.items()}
    for name, samples in BACKGROUNDS.items():
        mc[name] = load(samples)
    for name, rows in mc.items():
        rows.update(sf.event(rows))
    in_window = lambda r: (r["m4l_refit"] > lo) & (r["m4l_refit"] < hi)  # noqa: E731
    # Binnings: m4l bins uniform in the window; D bins the quintiles of the signal D in the window; error bins the
    # quintiles of the signal relative error per final state.
    m_edges = np.linspace(lo, hi, M_BINS + 1)
    sig_all = {k: np.concatenate([mc[m][k] for m in ("ggH", "VBF", "VH")]) for k in ("m4l_refit", "d_bkg_kin", "w", "final_state",
                                                                                    "m4l_refit_err", "m4l_err", "m4l")}
    sw = in_window(sig_all)
    d_edges = np.quantile(sig_all["d_bkg_kin"][sw], np.linspace(0, 1, D_BINS + 1))
    d_edges[0], d_edges[-1] = -1.5, 1.0 + 1e-9
    model = {"schema": "h4l_v3_model/1", "select_scan": str(PRODUCTION / "h4l_select" / args.select / "scan.json"),
             "select_plan_sha256": scan["plan_sha256"], "signal_model": str(PRODUCTION / "signal_model" / args.select / args.signal / "signal_model.json"),
             "zx": str(PRODUCTION / "backgrounds" / args.select / args.zx / "zx.json"), "sf_payload": sf.path, "sf_fallbacks": sf.fallbacks,
             "window": [lo, hi], "categories": order, "channel_kind": args.channels, "luminosity_fb": io.lumi_fb(scan),
             "m_edges": m_edges.tolist(), "d_edges": d_edges.tolist(), "final_states": {}, "channels": {}, "eff_correction": {}}
    if args.channels in fid.DIFFERENTIAL:
        model["observable_edges"] = [e if np.isfinite(e) else None for e in fid.DIFFERENTIAL[args.channels]]
    events = {}
    # Z+X category fractions from the OS-method control rows (signed weights), per final state.
    cr = io.read(scan, ["pseudo_data"], "CR", ["cr_type", "region", "final_state", "l_pdg", "l_pt", "l_eta", "l_eta_sc", "l_pass",
                                                "m4l", "m4l_refit", "d_bkg_kin", "m4l_refit_err", "pt4l", "jet_pt", "jet_eta"] +
                 cats.CATEGORY_BRANCHES)
    cr["category"] = channel_bin(cr)
    import zx_estimate as zxe
    # The fake-rate bins of the Z+X estimate (they follow its selection).
    zxe.FR_PT = {int(k): [x if x is not None else math.inf for x in v] for k, v in zx["fake_rate_bins"].items()}
    rates_os = {int(k): np.array(v) for k, v in zx["fake_rates_os"].items()}
    rates_ss = {int(k): np.array(v) for k, v in zx["fake_rates_ss"].items()}
    osr = {k: v[(cr["region"] == 0) & (cr["cr_type"] <= 1)] for k, v in cr.items()}
    f3 = zxe.lookup(rates_os, osr["l_pdg"][:, 2], osr["l_pt"][:, 2], zxe.abs_eta_of(osr["l_pdg"][:, 2], osr["l_eta"][:, 2], osr["l_eta_sc"][:, 2]))
    f4 = zxe.lookup(rates_os, osr["l_pdg"][:, 3], osr["l_pt"][:, 3], zxe.abs_eta_of(osr["l_pdg"][:, 3], osr["l_eta"][:, 3], osr["l_eta_sc"][:, 3]))
    fail3 = osr["l_pass"][:, 2] == 0
    ff = np.where(fail3, f3, f4)
    osw = np.where(osr["cr_type"] == 1, ff / (1 - ff), -f3 * f4 / ((1 - f3) * (1 - f4)))
    ssr = {k: v[(cr["region"] == 0) & (cr["cr_type"] == 2)] for k, v in cr.items()}
    s3 = zxe.lookup(rates_ss, ssr["l_pdg"][:, 2], ssr["l_pt"][:, 2], zxe.abs_eta_of(ssr["l_pdg"][:, 2], ssr["l_eta"][:, 2], ssr["l_eta_sc"][:, 2]))
    s4 = zxe.lookup(rates_ss, ssr["l_pdg"][:, 3], ssr["l_pt"][:, 3], zxe.abs_eta_of(ssr["l_pdg"][:, 3], ssr["l_eta"][:, 3], ssr["l_eta_sc"][:, 3]))
    ssw = s3 * s4
    total_sf, total_nosf = {}, {}

    def acceptance_ratios(fs_name, mode, channel_label):
        """A x eff(m_H) / A x eff(125) of a channel: the signal model's category value, or (for the inclusive and the
        observable channels) the yield-weighted mean over the signal model's categories of the final state."""
        ys = signal_model["final_states"][fs_name]["yields"][mode]
        if args.channels == "categories":
            return ys[channel_label]["acceptance_ratio_scale"], ys[channel_label]["acceptance_ratio_shift"]
        total = sum(v["yield_125"] for v in ys.values())
        out = []
        for key in ("acceptance_ratio_scale", "acceptance_ratio_shift"):
            out.append((sum(v["yield_125"] * np.array(v[key]) for v in ys.values()) / total).tolist() if total > 0 else
                       [1.0] * len(signal_model["mh_grid"]))
        return out[0], out[1]

    for fs, fs_name in FINAL_STATES.items():
        sm = signal_model["final_states"][fs_name]
        fs_entry = {"signal": {v: sm[v] for v in ("m4l_refit", "m4l")}, "backgrounds": {}}
        for name in BACKGROUNDS:
            r = mc[name]
            sel = r["final_state"] == fs
            fs_entry["backgrounds"][name] = {"bernstein_m4l_refit": fit_bernstein(r["m4l_refit"][sel], (r["w"] * r["sf"])[sel], lo, hi),
                                             "bernstein_m4l": fit_bernstein(r["m4l"][sel], (r["w"] * r["sf"])[sel], lo, hi)}
        zshape = zx["shapes"][fs_name]
        zc = zx["combined"][fs_name]
        if list(zc["window"]) != [lo, hi]:
            raise SystemExit(f"the Z+X window {zc['window']} differs from the model window {[lo, hi]}")
        zx_window = zc["value"] * zc["window_fraction"]
        fs_entry["zx"] = {"shape": {k: zshape[k] for k in ("form", "range", "mpv", "width", "exp_fraction", "exp_slope")},
                          "yield_window": zx_window, "fraction_window": zc["window_fraction"],
                          "kappa_low": zc["kappa_low"], "kappa_high": zc["kappa_high"]}
        # D templates and relative-error templates.
        sfs = sw & (sig_all["final_state"] == fs)
        e_edges, e_edges_raw = {}, {}
        for variable, error in (("m4l_refit", "m4l_refit_err"), ("m4l", "m4l_err")):
            rel = sig_all[error][sfs] / sig_all[variable][sfs]
            edges = np.quantile(rel, np.linspace(0, 1, E_BINS + 1))
            # The raw quantile edges (the 0.5 and 99.5 % points at the ends) serve the generation of
            # pseudo-experiments; the lookup edges are open at both ends.
            raw = edges.copy()
            raw[0], raw[-1] = np.quantile(rel, 0.005), np.quantile(rel, 0.995)
            e_edges_raw[variable] = raw.tolist()
            edges[0], edges[-1] = 0.0, 10.0
            e_edges[variable] = edges.tolist()
        fs_entry["e_edges"] = e_edges
        fs_entry["e_edges_raw"] = e_edges_raw
        dt = {"signal": template_2d(sig_all["m4l_refit"][sfs], sig_all["d_bkg_kin"][sfs], sig_all["w"][sfs], m_edges, d_edges)}
        et = {"signal": {v: template_1d((sig_all[e] / sig_all[v])[sfs], sig_all["w"][sfs], np.array(e_edges[v]))
                         for v, e in (("m4l_refit", "m4l_refit_err"), ("m4l", "m4l_err"))}}
        for name in BACKGROUNDS:
            r = mc[name]
            sel = in_window(r) & (r["final_state"] == fs)
            dt[name] = template_2d(r["m4l_refit"][sel], r["d_bkg_kin"][sel], (r["w"] * r["sf"])[sel], m_edges, d_edges)
            et[name] = {v: template_1d((r[e] / r[v])[sel], (r["w"] * r["sf"])[sel], np.array(e_edges[v]))
                        for v, e in (("m4l_refit", "m4l_refit_err"), ("m4l", "m4l_err"))}
        ssel = (ssr["final_state"] == fs) & (ssr["m4l_refit"] > lo - 20) & (ssr["m4l_refit"] < hi + 40)
        mm = np.clip(ssr["m4l_refit"][ssel], lo + 1e-6, hi - 1e-6)
        dt["zx"] = template_2d(mm, ssr["d_bkg_kin"][ssel], ssw[ssel], m_edges, d_edges)
        et["zx"] = {"m4l_refit": template_1d((ssr["m4l_refit_err"] / ssr["m4l_refit"])[ssel], ssw[ssel], np.array(e_edges["m4l_refit"])),
                    "m4l": template_1d((ssr["m4l_refit_err"] / ssr["m4l_refit"])[ssel], ssw[ssel], np.array(e_edges["m4l"]))}
        fs_entry["d_templates"] = dt
        fs_entry["e_templates"] = et
        model["final_states"][fs_name] = fs_entry
        # Channels.
        oin = osr["final_state"] == fs
        os_total = float(np.sum(osw[oin]))
        # Z+X category fractions (signed OS-method weights); a negative fraction is set to zero and the fractions are
        # renormalized, so the final state's Z+X total is preserved (AGENTS.md signed-template rule).
        raw_frac = [float(np.sum(osw[oin & (osr["category"] == c)])) / os_total if os_total != 0 else (1.0 if cat == "Untagged" else 0.0)
                    for c, cat in enumerate(order)]
        clipped = [max(f, 0.0) for f in raw_frac]
        norm = sum(clipped)
        final_frac = [f / norm for f in clipped] if norm > 0 else [1.0 if cat == "Untagged" else 0.0 for cat in order]
        for c, cat in enumerate(order):
            ch = f"{fs_name}_{cat}"
            entry = {"final_state": fs_name, "category": cat, "signal": {}, "backgrounds": {}}
            for name in list(SIGNAL.values()) + list(BACKGROUNDS):
                r = mc[name]
                sel = in_window(r) & (r["final_state"] == fs) & (r["category"] == c)
                w = r["w"][sel]
                y = float(np.sum(w * r["sf"][sel]))
                y_nosf = float(np.sum(w))
                # Efficiency nuisances: the coherent fit-model part in quadrature with the bin-by-bin independent
                # statistical part (h4l_sf.efficiency_variation).
                mu_up, mu_down, mu_parts = h4l_sf.efficiency_variation(w, r, sel, "muon")
                e_up, e_down, e_parts = h4l_sf.efficiency_variation(w, r, sel, "electron")
                block = {"yield": y, "yield_no_sf": y_nosf, "mc_stat": float(np.sqrt(np.sum((w * r["sf"][sel]) ** 2))),
                         "sf_muon": [mu_up, mu_down], "sf_electron": [e_up, e_down],
                         "sf_parts": {"muon": mu_parts, "electron": e_parts}}
                if name in SIGNAL.values():
                    block["acceptance_ratio_scale"], block["acceptance_ratio_shift"] = acceptance_ratios(fs_name, name, cat)
                    block["breakdown"] = breakdown(name, r, sel, w * r["sf"][sel])
                    entry["signal"][name] = block
                    total_sf[fs_name] = total_sf.get(fs_name, 0.0) + y
                    total_nosf[fs_name] = total_nosf.get(fs_name, 0.0) + y_nosf
                else:
                    entry["backgrounds"][name] = block
            entry["backgrounds"]["zx"] = {"yield": final_frac[c] * zx_window, "category_fraction": final_frac[c],
                                          "category_fraction_signed": raw_frac[c], "negative_fraction_zeroed": raw_frac[c] < 0}
            dsel = in_window(data) & (data["final_state"] == fs) & (data["category"] == c)
            for v in ("m4l_refit", "m4l_refit_err", "m4l", "m4l_err", "d_bkg_kin"):
                events[f"{ch}__{v}"] = data[v][dsel]
            entry["observed"] = int(dsel.sum())
            model["channels"][ch] = entry
        model["eff_correction"][fs_name] = total_sf[fs_name] / total_nosf[fs_name] if total_nosf.get(fs_name) else None
    model["eff_correction"]["inclusive"] = sum(total_sf.values()) / sum(total_nosf.values())
    out_dir.mkdir(parents=True)
    np.savez(out_dir / "events.npz", **events)
    (out_dir / "model.json").write_text(json.dumps(model, indent=1) + "\n", encoding="utf-8")
    for fs_name in FINAL_STATES.values():
        chans = [c for c in model["channels"].values() if c["final_state"] == fs_name]
        s = sum(sum(b["yield"] for b in c["signal"].values()) for c in chans)
        b = {n: sum(c["backgrounds"][n]["yield"] for c in chans) for n in ("qqZZ", "ggZZ", "zx")}
        n = sum(c["observed"] for c in chans)
        print(f"[model] {fs_name}: signal {s:.2f}, qqZZ {b['qqZZ']:.2f}, ggZZ {b['ggZZ']:.2f}, Z+X {b['zx']:.2f}; observed {n}; "
              f"eff_correction {model['eff_correction'][fs_name]:.4f}", flush=True)
    print(f"[model] eff_correction inclusive {model['eff_correction']['inclusive']:.4f}\n[model] {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
