"""Inputs of the Z1 kinematic refit and of the per-event mass uncertainty.

    pixi run py -- analysis_v3/reconstruction/scripts/make_refit_inputs.py --skims v2 --version v1

1. The true Z1 line shape of AN-16-442 5.4 ("true gen level Z1 lineshape from
   the SM Higgs boson sample"): from the GenTable of the ggH (m_H = 125 GeV)
   skim, the mass of the dressed (photons within the fiducial dressing cone)
   same-flavour opposite-sign pair of Higgs-decay leptons closest to m_Z,
   weighted with genWeight, estimated in 30-125 GeV with an adaptive Gaussian
   kernel density (Abramson: h_i = h0 (f_pilot(m_i) / g)^(-1/2), h0 = 0.25 GeV,
   clipped to 0.15-4 GeV, f_pilot the 1 GeV-smoothed histogram, g its
   geometric mean), tabulated on a 0.05 GeV grid.  v1 used the histogram
   smoothed by a fixed 0.2 GeV (rough in the off-shell tail); a parametric
   BW x exp(Chebyshev) could not describe the peak and the kinematic edge
   m_Z1 < m_H - m_Z2 (chi2 46,225 / 374).
2. The FSR photon pT resolution ("a parametrization obtained from
   simulation", JHEP 11 (2017) 047 section 5): every AN FSR photon of the
   signal skims (pT > 2, |eta| < 2.4, relIso < 1.8, dR/ET^2 < 0.012,
   dR < 0.5 to its muon) matched to the closest status-1 generator photon
   within dR < 0.1 and 0.5 < pT_gen / pT_reco < 2; the relative resolution
   in bins of reco pT is the half width of the central 68.3 % interval of
   (pT_reco - pT_gen) / pT_gen, fitted with sigma / pT = a / sqrt(pT) (+) b.
Writes production_v3/h4l_reco/refit/<version>/refit_inputs.json and plots.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_style  # noqa: E402

Z_MASS = 91.1876


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def dressed_z1_masses(path: str) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(path) as f:
        a = f["GenTable"].arrays(["weight", "d_pt", "d_eta", "d_phi", "d_mass", "d_pdg", "d_from_h"], library="np")
    masses, weights = [], []
    for i in range(len(a["weight"])):
        sel = np.nonzero(a["d_from_h"][i])[0]
        if len(sel) != 4:
            continue
        pt, eta, phi, mass, pdg = (a[k][i][sel].astype(float) for k in ("d_pt", "d_eta", "d_phi", "d_mass", "d_pdg"))
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        e = np.sqrt(px ** 2 + py ** 2 + pz ** 2 + mass ** 2)
        best = None
        for x in range(4):
            for y in range(x + 1, 4):
                if pdg[x] != -pdg[y]:
                    continue
                m2 = (e[x] + e[y]) ** 2 - (px[x] + px[y]) ** 2 - (py[x] + py[y]) ** 2 - (pz[x] + pz[y]) ** 2
                m = math.sqrt(max(m2, 0.0))
                if best is None or abs(m - Z_MASS) < abs(best - Z_MASS):
                    best = m
        if best is not None:
            masses.append(best)
            weights.append(a["weight"][i])
    return np.array(masses), np.array(weights)


def photon_pairs(path: str) -> tuple[np.ndarray, np.ndarray]:
    branches = ["FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi", "FsrPhoton_relIso03", "FsrPhoton_dROverEt2",
                "FsrPhoton_muonIdx", "Muon_eta", "Muon_phi", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_pdgId",
                "GenPart_status"]
    reco, gen = [], []
    with uproot.open(path) as f:
        for a in f["Events"].iterate(branches, library="np", step_size=20000):
            for i in range(len(a["FsrPhoton_pt"])):
                ph_pt, ph_eta, ph_phi = a["FsrPhoton_pt"][i], a["FsrPhoton_eta"][i], a["FsrPhoton_phi"][i]
                if len(ph_pt) == 0:
                    continue
                g = (a["GenPart_pdgId"][i] == 22) & (a["GenPart_status"][i] == 1)
                g_pt, g_eta, g_phi = a["GenPart_pt"][i][g], a["GenPart_eta"][i][g], a["GenPart_phi"][i][g]
                for k in range(len(ph_pt)):
                    m = a["FsrPhoton_muonIdx"][i][k]
                    if m < 0 or m >= len(a["Muon_eta"][i]):
                        continue
                    dphi_mu = (ph_phi[k] - a["Muon_phi"][i][m] + np.pi) % (2 * np.pi) - np.pi
                    dr_mu = math.hypot(ph_eta[k] - a["Muon_eta"][i][m], dphi_mu)
                    if not (ph_pt[k] > 2 and abs(ph_eta[k]) < 2.4 and a["FsrPhoton_relIso03"][i][k] < 1.8
                            and a["FsrPhoton_dROverEt2"][i][k] < 0.012 and dr_mu < 0.5):
                        continue
                    if len(g_pt) == 0:
                        continue
                    dphi = (g_phi - ph_phi[k] + np.pi) % (2 * np.pi) - np.pi
                    dr = np.hypot(g_eta - ph_eta[k], dphi)
                    ok = (dr < 0.1) & (g_pt > 0.5 * ph_pt[k]) & (g_pt < 2.0 * ph_pt[k])
                    if not ok.any():
                        continue
                    j = np.argmin(np.where(ok, dr, np.inf))
                    reco.append(ph_pt[k])
                    gen.append(g_pt[j])
    return np.array(reco), np.array(gen)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skims", required=True)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    out_dir = PRODUCTION / "h4l_reco" / "refit" / args.version
    target = out_dir / "refit_inputs.json"
    if target.exists():
        fail(f"{target} exists; choose a new --version")
    plan = json.loads((PRODUCTION / "skims" / args.skims / "plan.json").read_text(encoding="utf-8"))
    signal = [t for t in plan["tasks"] if t["config"]["role"] == "signal"]
    ggh = [t for t in signal if t["config"].get("mode") == "ggH"]
    if len(ggh) != 1:
        fail("expected exactly one ggH skim task")
    masses, weights = dressed_z1_masses(ggh[0]["outputs"]["root"])
    keep = (masses > 30.0) & (masses < 125.0)
    masses, weights = masses[keep], weights[keep]
    fit_edges = np.round(np.arange(30.0, 125.0 + 1e-9, 0.25), 6)
    hist, _ = np.histogram(masses, fit_edges, weights=weights)
    fit_centres = 0.5 * (fit_edges[1:] + fit_edges[:-1])
    # Pilot density: the 0.25 GeV histogram smoothed by 1 GeV.
    kx = np.arange(-16, 17) * 0.25
    kernel = np.exp(-0.5 * (kx / 1.0) ** 2)
    pilot_hist = np.convolve(hist, kernel / kernel.sum(), mode="same") / (hist.sum() * 0.25)
    pilot = np.interp(masses, fit_centres, np.maximum(pilot_hist, 1e-6))
    g = math.exp(np.average(np.log(pilot), weights=weights))
    bandwidth = np.clip(0.25 * np.sqrt(g / pilot), 0.15, 4.0)
    edges = np.round(np.arange(30.0, 125.0 + 1e-9, 0.05), 6)
    centres = 0.5 * (edges[1:] + edges[:-1])
    density = np.zeros_like(centres)
    for start in range(0, len(masses), 4000):
        m, w, h = masses[start:start + 4000], weights[start:start + 4000], bandwidth[start:start + 4000]
        z = (centres[None, :] - m[:, None]) / h[:, None]
        density += np.sum(w[:, None] * np.exp(-0.5 * z * z) / (h[:, None] * math.sqrt(2 * math.pi)), axis=0)
    density = density / (density.sum() * 0.05)
    mu_fit = np.interp(fit_centres, centres, density) * hist.sum() * 0.25
    chi2 = float(np.sum((hist - mu_fit) ** 2 / np.maximum(mu_fit, 1e-9)))
    coefficients = np.array([])
    # FSR photons of all signal samples.
    reco, gen = [], []
    for task in signal:
        r, g = photon_pairs(task["outputs"]["root"])
        reco.append(r)
        gen.append(g)
    reco, gen = np.concatenate(reco), np.concatenate(gen)
    rel = (reco - gen) / gen
    pt_edges = np.array([2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100])
    points = []
    for lo, hi in zip(pt_edges[:-1], pt_edges[1:]):
        sel = (reco >= lo) & (reco < hi)
        if sel.sum() < 100:
            continue
        q16, q50, q84 = np.percentile(rel[sel], [15.865, 50, 84.135])
        points.append({"pt_low": float(lo), "pt_high": float(hi), "mean_pt": float(reco[sel].mean()), "n": int(sel.sum()),
                       "median": float(q50), "sigma": float(0.5 * (q84 - q16))})
    x = np.array([p["mean_pt"] for p in points])
    y = np.array([p["sigma"] for p in points])
    # sigma^2 = a^2 / pT + b^2: linear least squares in (1/pT, 1).
    design = np.vstack([1.0 / x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(design, y ** 2, rcond=None)
    a2, b2 = max(coef[0], 0.0), max(coef[1], 0.0)
    result = {"schema": "h4l_v3_refit_inputs/1", "version": args.version,
              "z1_lineshape": {"source": ggh[0]["outputs"]["root"], "events": int(len(masses)),
                               "definition": "dressed H-decay SF OS pair closest to m_Z, genWeight; adaptive Gaussian KDE (h0 0.25 GeV, 0.15-4 GeV) in 30-125 GeV",
                               "chi2_vs_histogram": chi2, "histogram_bins": int(len(hist)),
                               "edges": edges.tolist(), "density": density.tolist()},
              "fsr_photon_resolution": {"sources": [t["outputs"]["root"] for t in signal], "photons": int(len(reco)),
                                        "form": "sigma(pT)/pT = sqrt(a^2 / pT + b^2), pT in GeV",
                                        "a": math.sqrt(a2), "b": math.sqrt(b2), "points": points}}
    out_dir.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    temporary.rename(target)
    # Plots.
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    norm = hist.sum() * 0.25
    ax.errorbar(fit_centres, hist / norm, np.sqrt(np.maximum(hist, 1)) / norm, fmt="o", ms=2, color="black",
                label="generator (dressed H-decay leptons)")
    ax.plot(fit_centres, mu_fit / norm, color="#e42536", lw=1.2, label=f"adaptive KDE, chi2 {chi2:.0f}/{len(hist)}")
    ax.set_yscale("log")
    ax.set_ylabel("density [1/GeV]")
    ax.legend(fontsize=8)
    pull = (hist - mu_fit) / np.sqrt(np.maximum(mu_fit, 1e-9))
    bx.plot(fit_centres, pull, ".", color="black", ms=2)
    bx.axhline(0, color="grey", lw=0.5)
    bx.set_ylim(-5, 5)
    bx.set_xlabel("m(Z1) generator [GeV]")
    bx.set_ylabel("pull")
    h4l_style.save(fig, out_dir / "z1_lineshape")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(x, y, fmt="o", color="black", label="FSR photons (signal MC)")
    grid = np.linspace(2, 100, 200)
    ax.plot(grid, np.sqrt(a2 / grid + b2), color="#e42536",
            label=f"sqrt(a^2/pT + b^2), a = {math.sqrt(a2):.3f}, b = {math.sqrt(b2):.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("FSR photon pT [GeV]")
    ax.set_ylabel("sigma(pT)/pT")
    ax.legend()
    h4l_style.save(fig, out_dir / "fsr_photon_resolution")
    plt.close(fig)
    print(f"[refit] {target}: {len(masses)} ggH events, {len(reco)} matched FSR photons, a = {math.sqrt(a2):.4f}, "
          f"b = {math.sqrt(b2):.4f}")
    for p in points:
        print(f"  pT {p['pt_low']:g}-{p['pt_high']:g}: n {p['n']}, median {p['median']:+.4f}, sigma {p['sigma']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
