"""The analysis note (TeX, pdflatex) generated from the result files (stage 10).

    pixi run py -- analysis_v3/documentation/make_an.py --select v5 --tag r1 --model-cat model_cat_r0 \
        --tnp v7/fits_i1c --calibration v4/nominal --zx zx_v2 --signal sm_v2 --dist dist_r1 \
        --toys <model>/<label> ... --out deliverables/AN_h4l_ul16_pfnano_v1 [--truth truth.json]

Every number of the note is read from the result files (results_full.json, RESULT.json, eval_selfreport.json,
yields.json, zx.json, the calibration and tag-and-probe reports, the closures, the toys); every figure is
copied from production_v3 into <out>/figures (never referenced in place).  The tag-and-probe galleries are
included page by page (downscaled JPEG copies), the calibration galleries as PDF pages.  The note is compiled
twice with /usr/bin/pdflatex.  With --truth (after the user's unblinding) the truth-comparison section is added.
Never overwrites an existing output directory.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
PRODUCTION = REPO / "production_v3"
FINAL_STATES = ("4mu", "4e", "2e2mu")
FS_TEX = {"4mu": r"$4\mu$", "4e": r"$4e$", "2e2mu": r"$2e2\mu$"}


def tex_escape(s: str) -> str:
    return (str(s).replace("\\", r"\textbackslash{}").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
            .replace("#", r"\#").replace("$", r"\$").replace("{", r"\{").replace("}", r"\}").replace("~", r"\~{}")
            .replace("^", r"\^{}").replace("<", r"$<$").replace(">", r"$>$"))


def asym(p: dict, digits: int = 2) -> str:
    """value^{+up}_{-down} of a result block with a 'total' interval."""
    lo, hi = p["total"]
    return f"${p['value']:.{digits}f}^{{+{hi:.{digits}f}}}_{{-{abs(lo):.{digits}f}}}$"


def asym_stat_syst(p: dict, digits: int = 2) -> str:
    lo, hi = p["stat"]
    return (f"${p['value']:.{digits}f}^{{+{hi:.{digits}f}}}_{{-{abs(lo):.{digits}f}}}\\,(\\text{{stat}})"
            f"\\pm{p.get('syst', float('nan')):.{digits}f}\\,(\\text{{syst}})$")


class Figures:
    def __init__(self, out: Path):
        self.dir = out / "figures"
        self.dir.mkdir(parents=True)
        self.count = 0

    def add(self, source: Path, name: str | None = None) -> str | None:
        if not source.exists():
            return None
        target = self.dir / (name or f"f{self.count:04d}_{source.name}")
        shutil.copy2(source, target)
        self.count += 1
        return f"figures/{target.name}"

    def add_jpeg(self, source: Path, scale: float = 0.5) -> str | None:
        """A downscaled JPEG copy of a PNG gallery page (the note stays of a manageable size)."""
        if not source.exists():
            return None
        from PIL import Image
        image = Image.open(source).convert("RGB")
        image = image.resize((int(image.width * scale), int(image.height * scale)))
        target = self.dir / f"g{self.count:04d}_{source.stem}.jpg"
        image.save(target, quality=82)
        self.count += 1
        return f"figures/{target.name}"


def figure(path: str | None, caption: str, width: str = "0.85\\textwidth") -> str:
    if path is None:
        return f"\\noindent\\textit{{[figure not available: {tex_escape(caption[:60])}]}}\n\n"
    return (f"\\begin{{figure}}[H]\\centering\\includegraphics[width={width}]{{{path}}}\n"
            f"\\caption{{{caption}}}\\end{{figure}}\n")


def figures_grid(paths: list, caption: str, width: str = "0.48\\textwidth") -> str:
    items = [p for p in paths if p]
    if not items:
        return ""
    body = "\n".join(f"\\includegraphics[width={width}]{{{p}}}" for p in items)
    return f"\\begin{{figure}}[H]\\centering\n{body}\n\\caption{{{caption}}}\\end{{figure}}\n"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("--select", "--tag", "--model-cat", "--tnp", "--calibration", "--zx", "--signal", "--dist", "--out"):
        parser.add_argument(name, required=True)
    parser.add_argument("--nm1", default=None)
    parser.add_argument("--toys", nargs="*", default=[], help="<model>/<label> of run_toys outputs (validation)")
    parser.add_argument("--paired", default=None, help="<model>/<label> of the paired-injection toys")
    parser.add_argument("--tnp-closure", nargs="*", default=[], help="production_v3/tnp/v2/<run>/<label> of the efficiency closures")
    parser.add_argument("--truth", type=Path, default=None)
    parser.add_argument("--pyhf-check", type=Path, default=None)
    args = parser.parse_args()
    out = REPO / args.out
    if out.exists():
        raise SystemExit(f"{out} exists; the note is never overwritten (choose a new version)")
    out.mkdir(parents=True)
    fig = Figures(out)
    res_dir = PRODUCTION / "results" / args.select / args.tag
    full = load(res_dir / "results_full.json")
    result = load(res_dir / "RESULT.json")
    selfreport = load(res_dir / "eval_selfreport.json") or full.get("declared_cuts_and_coverage")
    yields = load(PRODUCTION / "results" / args.select / args.dist / "yields.json")
    zx = load(PRODUCTION / "backgrounds" / args.select / args.zx / "zx.json")
    calib_dir = PRODUCTION / "calibration" / args.calibration
    tnp_dir = PRODUCTION / "tnp" / "v2" / args.tnp
    sf = load(tnp_dir / "sf.json")
    model = load(PRODUCTION / "inference" / args.select / args.model_cat / "model.json")
    manifests = {p.stem[3:]: load(p) for p in sorted((PRODUCTION / "manifests" / "v2").glob("mc_*.json"))}
    tex = []
    add = tex.append

    # ------------------------------------------------------------------ front matter
    add(r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=2.2cm]{geometry}
\usepackage{graphicx,amsmath,booktabs,longtable,float,caption,xcolor,pdfpages}
\usepackage[hidelinks]{hyperref}
\captionsetup{font=small}
\setlength{\parskip}{0.4em}\setlength{\parindent}{0pt}
\title{Measurement of Higgs boson properties in the $H\to ZZ^*\to 4\ell$ channel\\
with 20~fb$^{-1}$ of UL16 PFNano pseudo-data: a reference analysis}
\author{Private work}
\date{\today}
\begin{document}\maketitle
\begin{abstract}
""")
    b = full["inclusive"]["set_b"]
    add(f"The $H\\to ZZ^*\\to 4\\ell$ ($\\ell = e,\\mu$) analysis of AN-16-442 and JHEP 11 (2017) 047 is carried out on "
        f"{full['luminosity_fb']:g}~fb$^{{-1}}$ of UL16 PFNano pseudo-data with the delivered simulation. Every detector effect "
        f"is measured in situ: the lepton momentum scale and resolution from $Z\\to\\ell\\ell$ template fits, the lepton "
        f"efficiencies from tag-and-probe, the reducible background from control regions. The signal strength is "
        f"$\\mu = {b['data']['mu']['value']:.2f}^{{+{b['data']['mu']['total'][1]:.2f}}}_{{-{abs(b['data']['mu']['total'][0]):.2f}}}$ and the "
        f"Higgs boson mass $m_H = {b['data']['mH']['value']:.2f}^{{+{b['data']['mH']['total'][1]:.2f}}}_{{-{abs(b['data']['mH']['total'][0]):.2f}}}$~GeV "
        f"(benchmark result set, $m_H$ floating), with a local significance of {b['significance']['observed']:.1f} standard "
        f"deviations. Signal strengths per final state, category and production mode, stage-0 simplified template cross "
        f"sections, fiducial cross sections, the mass in 1D, 2D and 3D fits, the width and the $Z\\to4\\ell$ mass are "
        f"measured, and the whole chain is validated with pseudo-experiments and closure tests.\n\\end{{abstract}}\n\\tableofcontents\\newpage\n")

    # ------------------------------------------------------------------ 1 introduction
    add(r"""\section{Introduction}
This note documents a complete reference analysis of the $H\to ZZ^*\to4\ell$ decay on a pseudo-data set of 2016
(ultra-legacy) conditions in the PFNano format.  The method follows the CMS analysis note AN-16-442 (v8) and the
paper JHEP 11 (2017) 047 as closely as the inputs allow; where an identical implementation is impossible a documented
substitute is used (Section~\ref{sec:subst}).  Two result sets are reported: set (b), the benchmark-facing one, contains
only the uncertainties that are real for this pseudo-data (statistics, lepton scale and resolution, efficiency scale
factors, the reducible-background method, the $A\times\epsilon(m_H)$ extrapolation); set (a) adds the theory and
luminosity uncertainties of the paper.  The data are described only by their integrated luminosity.
""")

    # ------------------------------------------------------------------ 2 samples
    add(r"\section{Data and simulated samples}" + "\n")
    add(f"The data correspond to {full['luminosity_fb']:g}~fb$^{{-1}}$ at $\\sqrt{{s}}=13$~TeV (83 files, treated exactly as collision data). "
        "The simulation is the delivered MC; its normalization uses the generator weights: "
        "$w = w_{\\text{gen}}\\,\\sigma_{\\text{eff}} L / \\sum w_{\\text{gen}}$ summed over exactly the files present "
        "(Table~\\ref{tab:mc}); the effective cross sections include the decays and filters of the samples.\n")
    add(r"\begin{table}[H]\centering\small\caption{Simulated samples, effective cross sections and file counts.}\label{tab:mc}"
        r"\begin{tabular}{lrrl}\toprule sample & $\sigma_{\text{eff}}$ [pb] & files & role\\\midrule" + "\n")
    roles = {"ZZTo4L": r"$q\bar q\to ZZ\to4\ell$", "GGZZ4Mu": r"$gg\to ZZ\to4\mu$", "GGZZ4E": r"$gg\to ZZ\to4e$",
             "GGZZ2E2Mu": r"$gg\to ZZ\to2e2\mu$", "DYJetsToLL": r"$Z/\gamma^*$+jets (Z+X closure, T\&P)",
             "TTBar": r"$t\bar t\to2\ell2\nu$ (Z+X closure)", "GluGluToHToZZ_M125": "ggH, $m_H=125$ GeV",
             "VBF_HToZZ_M125": "VBF, $m_H=125$ GeV", "VHToZZ_M125": "ZH + WH, $m_H=125$ GeV"}
    for name, m in manifests.items():
        add(f"{tex_escape(name)} & {m['normalization']['sigma_eff_pb']:g} & {len(m['files'])} & {roles.get(name, '')}\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")

    # ------------------------------------------------------------------ 3 objects
    add(r"""\section{Objects and event selection}
\textbf{Leptons.} Loose muons: global or tracker muons with $p_T>5$~GeV, $|\eta|<2.4$, $|d_{xy}|<0.5$~cm,
$|d_z|<1$~cm; loose electrons: $p_T>7$~GeV, $|\eta|<2.5$ with the same impact-parameter requirements.  Tight: muons
passing the PF identification (the tracker high-$p_T$ identification above 200~GeV), electrons passing
mvaFall17V2noIso WPL (the AN BDT is not available: substitution), $|\mathrm{SIP}_{3D}|<4$ and the FSR-subtracted
relative isolation ($\Delta R=0.3$) below 0.35.  FSR photons (muons only, the NanoAOD FSR collection) are attached as
in the AN ($\Delta R/E_T^2<0.012$, relative isolation $<1.8$) and removed from the isolation of every loose lepton.
Ghost and cross cleaning follow the AN (substitute: $\Delta R$-based ghost removal).

\textbf{Candidates.} $Z_1$: the opposite-sign same-flavour pair closest to $m_Z$ with $40<m_{Z_1}<120$~GeV; $Z_2$:
the other pair with $10<m_{Z_2}<120$~GeV (optimized on MC, the paper uses 12~GeV); leading and subleading lepton
$p_T>20$ and $10$~GeV; $m_{\ell\ell}>4$~GeV for every opposite-sign pair; $\Delta R>0.02$ between leptons; the smart
cut; $m_{4\ell}>70$~GeV; among several candidates the one with $Z_1$ closest to $m_Z$ and the largest
$D_{\text{bkg}}^{\text{kin}}$ (MELA).  The $Z\to4\ell$ region uses $m_{Z_2}>4$~GeV.

\textbf{Jets and categories.} Jets with $p_T>30$~GeV, $|\eta|<4.7$, cleaned against the leptons and FSR photons;
DeepJet medium b tagging.  Seven categories (untagged, VBF-1jet, VBF-2jet, VH-hadronic, VH-leptonic, VH-MET, ttH) with
the paper-convention discriminants $D'=1/(1+c\,P_{\text{bkg}}/P_{\text{sig}})$ from the MELA probabilities
($c_{2\text{jet}}=0.1$, $c_{1\text{jet}}=0.0014$, $c_{WH}=0.16$, $c_{ZH}=0.4$, working points 0.5, VH-MET with
$E_T^{\text{miss}}>90$~GeV; constants optimized on MC).
""")
    if yields:
        add(yields_tables(yields))
    dist = PRODUCTION / "results" / args.select / args.dist / "plots"
    add(figures_grid([fig.add(dist / n) for n in ("m4l_full.png", "m4l_low.png")],
                     "The $m_{4\\ell}$ distribution, full and low-mass range: data against the expectation (SM at 125 GeV)."))
    add(figures_grid([fig.add(dist / f"m4l_low_{fs}.png") for fs in FINAL_STATES], "The $m_{4\\ell}$ distribution per final state.", "0.32\\textwidth"))
    add(figures_grid([fig.add(p) for p in sorted(dist.glob("m4l_refit_*.png"))], "The $m_{4\\ell}$ distribution per category (105-140 GeV).", "0.32\\textwidth"))
    add(figures_grid([fig.add(dist / n) for n in ("mz2_vs_mz1.png", "dkin_vs_m4l.png")],
                     "Left: $m_{Z_2}$ against $m_{Z_1}$ in 118-130 GeV. Right: $D_{\\text{bkg}}^{\\text{kin}}$ against $m_{4\\ell}$ with the per-event mass uncertainties."))
    add(figures_grid([fig.add(p) for p in sorted(dist.glob("disc_*.png"))], "The production discriminants in 118-130 GeV.", "0.32\\textwidth"))

    # ------------------------------------------------------------------ 4 calibration
    report = load(calib_dir / "plots" / "report.json")
    add(r"""\section{Lepton momentum calibration}\label{sec:calib}
\textbf{Control pairs.} $Z\to\mu\mu$ and $Z\to ee$ pairs of tight leptons (the analysis selection) with one leg matched
to a trigger object; the pair mass in 60-120~GeV.  \textbf{Bins.} Per leg the calibrated $p_T$ and $|\eta|$
($|\eta_{SC}|$ for electrons); the pair categories are the leg-bin combinations.  \textbf{Method (user decision).} The
per-lepton scale $s(p_T,|\eta|)$ and extra resolution $r(p_T,|\eta|)$ are determined by MC template fits of the data
in the leg-bin categories: the MC dilepton spectra are smeared per lepton, $p_T \to p_T(1 + r N)$ with one Gaussian
deviate per lepton, and the data corrected, $p_T\to p_T\,e^{-u}$, iterated to the fixed point (10 iterations, the
final values the mean of the last four).  BW$\otimes$DCB fits with a background validate the result and calibrate the
per-event mass uncertainties.  \textbf{Definitions.} Scale shift = data/MC $-1$ of the lepton momentum; smear = the
extra relative per-lepton $p_T$ width of the data (never the dilepton width); both quoted at $p_T=45$~GeV,
$|\eta|=1.2$ and in every bin.  \textbf{Uncertainties.} Statistical (the fit), iteration (spread of the last
iterations) and the method accuracy from the closures on independent MC halves with injected effects (the largest
closure-profile rms below 100~GeV for the scale; the largest smear-variance deviation for the smear).
\textbf{Application.} Data leptons are corrected and every MC lepton is smeared with its own deviate in every stage
(skims, tag-and-probe, selection).
""")
    add(r"\begin{table}[H]\centering\small\caption{Calibration at $p_T=45$~GeV, $|\eta|=1.2$.}\begin{tabular}{lcc}\toprule "
        r" & muons & electrons\\\midrule" + "\n")
    c = full.get("calibration_report_point") or full["calibration"]
    add(f"scale shift & ${c['muon']['scale_shift']['value']:.5f}\\pm{c['muon']['scale_shift']['unc']:.5f}$ & "
        f"${c['electron']['scale_shift']['value']:.5f}\\pm{c['electron']['scale_shift']['unc']:.5f}$\\\\\n")
    add(f"smear & ${c['muon']['smear']['value']:.4f}\\pm{c['muon']['smear']['unc']:.4f}$ & "
        f"${c['electron']['smear']['value']:.4f}\\pm{c['electron']['smear']['unc']:.4f}$\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")
    cp = calib_dir / "plots"
    add(figures_grid([fig.add(cp / f"zpeak_{f}_{w}.png") for f in ("mm", "ee") for w in ("before", "after")],
                     "Dilepton mass before and after the calibration (muons top, electrons bottom)."))
    add(figures_grid([fig.add(cp / f"profile_{v}_{f}.png") for f in ("mm", "ee") for v in ("pt", "eta")],
                     "Scale and smear profiles against $p_T$ and $|\\eta|$ (muons top, electrons bottom)."))
    add(figures_grid([fig.add(cp / n) for n in ("convergence.png", "pulls_mm.png", "pulls_ee.png")], "Convergence of the iteration and the fit pulls.", "0.32\\textwidth"))
    diag = calib_dir / "diag_final" / "plots"
    add(figures_grid([fig.add(diag / f"charge_eta_phi_{f}.png") for f in ("mm", "ee")], "Residual data/MC differences after the calibration in charge x $\\eta$ and $\\phi$ slices."))
    for name in ("closure_uniform_halves", "closure_sloped_halves"):
        cl = PRODUCTION / "calibration" / "v4" / name / "closure"
        add(figures_grid([fig.add(cl / f"closure_{f}.png") for f in ("mm", "ee")],
                         f"Calibration closure ({tex_escape(name)}): independent MC halves, the injected scale and smear recovered."))
    lam = calib_dir / "lambda_v2_r2" / "plots"
    add(figures_grid([fig.add(p) for p in sorted(lam.glob("*.png"))][:6], "Calibration of the per-event mass uncertainty ($\\lambda$).", "0.32\\textwidth"))

    # ------------------------------------------------------------------ 5 efficiencies
    add(r"""\section{Lepton efficiencies: tag-and-probe}
Tag: a tight, trigger-matched lepton; probe: the other leg passing the loose selection; the pair mass in 60-120~GeV.
The efficiency chain follows the selection (user decision): identification given loose, SIP given identification,
isolation given identification and SIP; the full selection is also measured directly as a cross-check.  The pass and
fail spectra are fitted simultaneously per probe bin ($p_T\times|\eta|$, the binning of efficiency\_tnp\_an): signal =
the DY-MC template of truth-matched prompt pairs of the bin convolved with a Gaussian, the fail template mixed with a
free pass-like fraction; background = CMSShape with the EGM ranges; alternatives (fit-model systematic): an analytic
DSCB signal and a Bernstein background.  $\mathrm{SF}=\epsilon_{\text{data}}/\epsilon_{\text{MC}}$, both from the fits;
generator-truth efficiencies are never used.  Every fit was inspected visually (the galleries of
Appendix~\ref{app:tnp}).
""")
    add(r"\begin{table}[H]\centering\small\caption{Full single-lepton SF at $p_T=45$~GeV, $|\eta|=1.2$ and the resulting "
        r"efficiency correction of the signal expectation.}\begin{tabular}{lc}\toprule & value\\\midrule" + "\n")
    for f in ("muon", "electron"):
        e = (full.get("efficiency_sf_report_point") or full["efficiency_sf"])[f]
        add(f"SF {f} & ${e['value']:.4f}\\pm{e['unc']:.4f}$\\\\\n")
    for k, v in full["eff_correction"].items():
        add(f"eff\\_correction {tex_escape(k)} & {v:.4f}\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")
    tp = tnp_dir / "plots"
    add(figures_grid([fig.add(tp / f"eff_{t}_{s}.png") for t in ("mm", "ee") for s in ("id", "sip", "iso", "full")],
                     "Tag-and-probe efficiencies (data and MC) per step (muons first, electrons second).", "0.24\\textwidth"))
    add(figures_grid([fig.add(p) for p in sorted(tp.glob("sf_full_*.png"))], "The full single-lepton scale factors."))
    for cdir in args.tnp_closure:
        cl = load(REPO / cdir / "closure.json")
        if cl:
            add(tnp_closure_text(cl, cdir))

    # ------------------------------------------------------------------ 6 backgrounds
    add(r"""\section{Backgrounds}
The irreducible $q\bar q\to ZZ$ and $gg\to ZZ$ backgrounds come from the simulation (with the scale factors).  The
reducible background (Z+X) is estimated with the opposite-sign (OS: 2P2F and 3P1F control regions, prompt ZZ
contribution subtracted) and same-sign (SS, with the conversion correction for electrons) methods of the AN; fake rates
from $Z+1$ loose-lepton events ($|m_{Z_1}-m_Z|<7$~GeV, $E_T^{\text{miss}}<25$~GeV), denominator the loose lepton
passing SIP (user decision).  The systematic uncertainty is the combined-final-state MC closure of both methods
(user decision): $\max(|1-r|,\sigma_r)$ = """ + f"{zx['relative_systematic']:.3f}" + r""".  The two methods are combined
per final state (inverse-variance weights, the envelope as the uncertainty); the shape is a Landau plus exponential fit
to the average of the normalized OS and SS distributions.
""")
    add(r"\begin{table}[H]\centering\small\caption{Z+X yields for $m_{4\ell}>70$~GeV.}\begin{tabular}{lccc}\toprule "
        r"final state & OS & SS & combined\\\midrule" + "\n")
    for fs in FINAL_STATES:
        cz = zx["combined"][fs]
        add(f"{FS_TEX[fs]} & ${cz['os']['value']:.2f}\\pm{cz['os']['total']:.2f}$ & ${cz['ss']['value']:.2f}\\pm{cz['ss']['total']:.2f}$ & "
            f"${cz['value']:.2f}$ [{cz['envelope'][0]:.2f}, {cz['envelope'][1]:.2f}]\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")
    zp = PRODUCTION / "backgrounds" / args.select / args.zx / "plots"
    add(figures_grid([fig.add(zp / f"fake_rate_{f}.png") for f in ("muon", "electron")], "Fake rates (OS window and SS method)."))
    add(figures_grid([fig.add(zp / f"zx_shape_{fs}.png") for fs in FINAL_STATES], "Z+X $m_{4\\ell}$ shapes per final state.", "0.32\\textwidth"))

    # ------------------------------------------------------------------ 7 signal model
    add(r"""\section{Signal model}
Per final state a double-sided Crystal Ball is fitted to the ggH + VBF simulation at $m_H=125$~GeV (with and without
the $Z_1$ refit); VH adds a non-resonant Landau part.  For the 3D fit the density is conditional on the relative
per-event mass uncertainty $D_{\text{mass}}=\delta m/m_{4\ell}$ with width $s\,D_{\text{mass}}\,m_H$.  The $m_H$
dependence (only 125~GeV samples exist) is obtained by rest-frame scaling of the four-lepton system, $k=m_H/125$,
which also gives the acceptance ratios; the alternative shift morphing defines the morphing nuisance.  The expected
yields are $\sigma_{\text{eff}}(m_H)\,L\,A\epsilon(m_H)$ with $\sigma_{\text{eff}}(m_H)$ from the YR4 table.
""")
    sp = PRODUCTION / "signal_model" / args.select / args.signal / "plots"
    add(figures_grid([fig.add(sp / f"signal_dcb_{fs}.png") for fs in FINAL_STATES], "Signal shapes per final state.", "0.32\\textwidth"))
    dm = PRODUCTION / "results" / args.select / "dmass_v1"
    for fs in FINAL_STATES:
        add(figure(fig.add(dm / f"dmass_{fs}.png"), f"Validation of the $D_{{\\text{{mass}}}}$-conditional signal model, {FS_TEX[fs]}: signal MC "
                   f"against the model per $D_{{\\text{{mass}}}}$ bin.", "\\textwidth"))

    # ------------------------------------------------------------------ 8 statistical model and systematics
    syst = load(res_dir / "results_full.json")["inclusive"]["set_b"]["impacts"]
    add(r"""\section{Statistical model and systematic uncertainties}
Extended unbinned likelihood per channel (final state $\times$ category) in $m_{4\ell}$ (with the $Z_1$ refit),
$D_{\text{bkg}}^{\text{kin}}$ (templates conditional on $m_{4\ell}$) and $D_{\text{mass}}$; $\mu$ unbounded, $m_H$
floating in 110-140~GeV; nuisances with unit-Gaussian constraints.  Stat/syst separation by explicit fixed-nuisance
impacts (the statistical interval: the same likelihood with the nuisances fixed at their best fit).
""")
    add(r"\begin{table}[H]\centering\small\caption{Impacts on $\mu$ and $m_H$ (set b, data; $+1\sigma$ / $-1\sigma$).}"
        r"\begin{tabular}{lcc}\toprule nuisance & $\Delta\mu$ & $\Delta m_H$ [GeV]\\\midrule" + "\n")
    for n, v in syst.items():
        add(f"{tex_escape(n)} & ${v['mu'][0]:+.4f}/{v['mu'][1]:+.4f}$ & ${v['mH'][0]:+.4f}/{v['mH'][1]:+.4f}$\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")
    cov = selfreport["systematics"]
    add(r"\begin{table}[H]\centering\small\caption{Coverage of the twelve systematic groups (set a).}\begin{tabular}{llp{9cm}}\toprule "
        r"group & covered & implementation\\\midrule" + "\n")
    for k, v in cov.items():
        add(f"{tex_escape(k)} & {v} & {tex_escape(selfreport['notes']['systematics_evidence'][k])}\\\\\n")
    add(r"\bottomrule\end{tabular}\end{table}" + "\n")

    # ------------------------------------------------------------------ 9 results
    add(results_section(full, result, fig, args))

    # ------------------------------------------------------------------ 10 validation
    add(validation_section(full, fig, args))

    # ------------------------------------------------------------------ 11 paper comparison
    add(r"\section{Comparison with the paper}" + "\n" + "The data sets differ (20~fb$^{-1}$ of pseudo-data here, 35.9~fb$^{-1}$ of 2016 data in the paper): "
        "the comparison shows the consistency of the methods and the expected scaling of the precision, not an agreement of the values.\n")
    add(figure(fig.add(res_dir / "plots" / "paper_comparison.png"), "Principal results of this analysis against JHEP 11 (2017) 047.", "\\textwidth"))

    # ------------------------------------------------------------------ 12 substitutions
    add(substitutions_section())

    # ------------------------------------------------------------------ 13 truth comparison
    if args.truth:
        add(truth_section(json.loads(args.truth.read_text(encoding="utf-8")), result, selfreport))

    # ------------------------------------------------------------------ appendices: galleries
    add(r"\appendix" + "\n" + r"\section{Tag-and-probe fit galleries}\label{app:tnp}" + "\n"
        "Every nominal fit (pass left, fail right; data, total model, signal and background) and the alternative models.\n")
    for kind in ("gallery", "gallery_alternatives"):
        for page in sorted((tp / kind).glob("*.png")):
            path = fig.add_jpeg(page, 0.5)
            add(f"\\begin{{figure}}[H]\\centering\\includegraphics[width=\\textwidth,height=0.85\\textheight,keepaspectratio]{{{path}}}"
                f"\\caption*{{{tex_escape(page.stem)}}}\\end{{figure}}\\clearpage\n")
    add(r"\section{Calibration fit galleries}" + "\n")
    for name in ("gallery_A_mm_iter10.pdf", "gallery_A_ee_iter10.pdf", "gallery_B_mm_iter10.pdf", "gallery_B_ee_iter10.pdf"):
        path = fig.add(cp / name)
        if path:
            add(f"\\includepdf[pages=-,fitpaper=false,scale=0.9,pagecommand={{\\thispagestyle{{empty}}}}]{{{path}}}\n")
    add(r"\end{document}" + "\n")
    (out / "main.tex").write_text("".join(tex), encoding="utf-8")
    for _ in range(2):
        run = subprocess.run(["/usr/bin/pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], cwd=out,
                             capture_output=True, text=True)
        if run.returncode != 0:
            (out / "pdflatex_error.log").write_text(run.stdout[-20000:], encoding="utf-8")
            raise SystemExit(f"pdflatex failed; see {out / 'pdflatex_error.log'}")
    print(f"[an] {out / 'main.pdf'}")
    return 0


def yields_tables(y: dict) -> str:
    out = []
    cols = ["ggH", "VBF", "VH", "signal", "qqZZ", "ggZZ", "zx", "total", "observed"]
    for window, label in (("m4l_gt_70", "$m_{4\\ell}>70$~GeV (paper Table 1)"), ("m4l_118_130", "$118<m_{4\\ell}<130$~GeV")):
        t1 = y.get("table1", {}).get(window)
        if not t1:
            continue
        out.append(r"\begin{table}[H]\centering\small\caption{Expected (SM, $m_H=125$~GeV) and observed yields, " + label + r".}"
                   r"\begin{tabular}{l" + "c" * len(cols) + r"}\toprule & " + " & ".join(tex_escape(c) for c in cols) + r"\\\midrule" + "\n")
        for fs, row in t1.items():
            out.append(f"{FS_TEX.get(fs, tex_escape(fs))} & " + " & ".join((f"{row[c]:.0f}" if c == "observed" else f"{row[c]:.2f}")
                       if isinstance(row.get(c), (int, float)) else "-" for c in cols) + "\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    t2 = y.get("table2")
    if t2:
        cats_order = list(t2)
        keys = sorted({k for row in t2.values() for k in row if isinstance(row[k], (int, float))})
        out.append(r"\begin{table}[H]\centering\scriptsize\caption{Expected and observed yields per category in 118-130~GeV (paper Table 2).}"
                   r"\begin{tabular}{l" + "c" * len(keys) + r"}\toprule & " + " & ".join(tex_escape(k) for k in keys) + r"\\\midrule" + "\n")
        for cat in cats_order:
            row = t2[cat]
            out.append(f"{tex_escape(cat)} & " + " & ".join(f"{row[k]:.2f}" if isinstance(row.get(k), (int, float)) else "-" for k in keys) + "\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    return "".join(out)


def tnp_closure_text(cl: dict, cdir: str) -> str:
    rows = []
    for f, v in cl["flavours"].items():
        rp = v["report_point_chain"]
        rows.append(f"{f} & {v['injected']:.3f} & ${rp['value']:.5f}\\pm{rp['error']:.5f}$ & {v['report_point_pull']:+.2f}\\\\")
    return (r"\begin{table}[H]\centering\small\caption{Efficiency closure (" + tex_escape(cdir.split('/')[-2]) +
            r"): MC halves, the data half thinned (identification kept with probability keep\_id); the recovered full SF at "
            r"the report point.}\begin{tabular}{lccc}\toprule flavour & injected & recovered & pull\\\midrule" + "\n" +
            "\n".join(rows) + "\n" + r"\bottomrule\end{tabular}\end{table}" + "\n")


def results_section(full: dict, result: dict, fig: Figures, args) -> str:
    b, a = full["inclusive"]["set_b"], full["inclusive"]["set_a"]
    out = [r"\section{Results}" + "\n", r"\subsection{Signal strength and mass}" + "\n"]
    out.append(r"\begin{table}[H]\centering\small\caption{Inclusive results (3D fit with the $Z_1$ refit, $m_H$ floating).}"
               r"\begin{tabular}{lcc}\toprule & $\mu$ & $m_H$ [GeV]\\\midrule" + "\n")
    out.append(f"set (b), observed & {asym_stat_syst(b['data']['mu'])} & {asym_stat_syst(b['data']['mH'])}\\\\\n")
    out.append(f"set (b), expected (pre-fit Asimov) & {asym(b['asimov_prefit']['mu'])} & {asym(b['asimov_prefit']['mH'])}\\\\\n")
    out.append(f"set (b), expected (post-fit Asimov) & {asym(b['asimov_postfit']['mu'])} & {asym(b['asimov_postfit']['mH'])}\\\\\n")
    out.append(f"set (a), observed & {asym_stat_syst(a['data']['mu'])} & {asym_stat_syst(a['data']['mH'])}\\\\\n")
    out.append(f"set (a), expected & {asym(a['asimov_prefit']['mu'])} & {asym(a['asimov_prefit']['mH'])}\\\\\n")
    p125 = full["paper_mu_125p09"]["data"]["mu"]
    out.append(f"$\\mu$ at $m_H=125.09$~GeV (set b) & {asym_stat_syst(p125)} & --\\\\\n")
    out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    out.append(f"The local significance is {b['significance']['observed']:.2f} (expected {b['significance']['expected']:.2f}); "
               f"discovery (local significance $\\geq5$): {result['discovery']}.  The saturated goodness of fit gives "
               f"$q={full['gof']['q_saturated']:.1f}$ over {full['gof']['n_bins']} bins, $p={full['gof']['p_value']:.3f}$.\n")
    out.append(r"\begin{table}[H]\centering\small\caption{$m_H$ in the fits of paper Table 6 (1D $L(m_{4\ell})$, 2D $L(m_{4\ell},D_{\text{mass}})$, "
               r"3D $L(m_{4\ell},D_{\text{mass}},D_{\text{bkg}}^{\text{kin}})$, with and without the $Z_1$ refit; 2D\_Dkin = $L(m_{4\ell},D_{\text{bkg}}^{\text{kin}})$ "
               r"as an extra), observed and expected (pre-fit Asimov) uncertainty.}"
               r"\begin{tabular}{lcc}\toprule fit & observed $m_H$ [GeV] & expected $\sigma$\\\midrule" + "\n")
    for k, v in full["table6"].items():
        out.append(f"{tex_escape(k)} & {asym(v['data']['mH'])} & $\\pm{v['asimov']['mH']['total_sym']:.3f}$\\\\\n")
    out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    mfs = full["mH_final_state"]
    out.append(r"\begin{table}[H]\centering\small\caption{Per final state: $\mu$ (one fit, three strengths) and $m_H$ (separate fits).}"
               r"\begin{tabular}{lcc}\toprule final state & $\mu$ & $m_H$ [GeV]\\\midrule" + "\n")
    for fs in FINAL_STATES:
        mu = full["mu_final_state"]["data"].get(f"mu_{fs}")
        out.append(f"{FS_TEX[fs]} & {asym(mu) if mu else '--'} & {asym(mfs['data']['separate_fits'][fs]['mH'])}\\\\\n")
    out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    out.append(f"Mutual compatibility of the final-state masses (one mass against three, $\\mu$ and nuisances profiled): "
               f"$p={mfs['compatibility']['p_value']:.3f}$; of the signal strengths: $p={full['mu_final_state']['compatibility']['p_value']:.3f}$.\n")
    out.append("The per-final-state, per-category, per-mode, $(\\mu_F,\\mu_V)$ and stage-0 results follow the paper (Eq.~10.1): "
               "$L(m_{4\\ell},D_{\\text{bkg}}^{\\text{kin}})$ at $m_H=125.09$~GeV with the signal strengths in $[0,20]$ (user decision); "
               "every reconstructed signal event of the delivered MC has $|y_H|<2.5$, so the stage-0 values equal the per-mode ones.\n")
    for block, title in (("mu_category", "per category"), ("mu_mode", "per production mode"), ("mu_fv", "$(\\mu_F,\\mu_V)$"),
                         ("stxs0", "stage-0 STXS ($|y_H|<2.5$), $\\sigma/\\sigma_{SM}$")):
        out.append(r"\begin{table}[H]\centering\small\caption{Signal strengths " + title + r".}\begin{tabular}{lcc}\toprule "
                   r" & observed & expected\\\midrule" + "\n")
        for k, v in full[block]["data"].items():
            e = full[block]["asimov"].get(k)
            out.append(f"{tex_escape(k)} & {asym(v)} & {asym(e) if e else '--'}\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    fid = full["fiducial"]
    out.append(r"\subsection{Fiducial cross sections}" + "\n")
    out.append(f"The fiducial volume of the paper (Table 4, dressed leptons); the fit uses $m_{{4\\ell}}$ without the refit, no categories and no "
               f"$D_{{\\text{{bkg}}}}^{{\\text{{kin}}}}$, the final-state fractions floating.  SM: {fid['sm']['sigma_fid_fb']:.3f}~fb.\n")
    out.append(r"\begin{table}[H]\centering\small\caption{Fiducial cross sections.}\begin{tabular}{lc}\toprule & $\sigma_{\text{fid}}$ [fb]\\\midrule" + "\n")
    for name, label in (("fid_int_fix", "integrated, $m_H=125.09$~GeV"), ("fid_int_prof", "integrated, $m_H$ profiled")):
        d = fid[name]["data"]
        out.append(f"{label} & ${d['sigma_fid_fb']:.3f}^{{+{d['total_fb'][1]:.3f}}}_{{-{abs(d['total_fb'][0]):.3f}}}$\\\\\n")
    for fs, v in fid["final_state"]["data"].items():
        out.append(f"{FS_TEX[fs]} & ${v['sigma_fid_fb']:.3f}^{{+{v['total_fb'][1]:.3f}}}_{{-{abs(v['total_fb'][0]):.3f}}}$\\\\\n")
    out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    if fid.get("final_state_compatibility"):
        out.append(f"Compatibility of the three final states: $p={fid['final_state_compatibility']['p_value']:.3f}$.\n")
    if fid.get("table5"):
        out.append(r"\begin{table}[H]\centering\small\caption{Paper Table 5 for the delivered samples: acceptance (relative to the "
                   r"sample's $\sigma_{\text{eff}}$, LO Pythia8), reconstruction efficiency of fiducial events and non-fiducial fraction.}"
                   r"\begin{tabular}{lccc}\toprule mode & $A_{\text{fid}}$ & $\epsilon$ & $f_{\text{nonfid}}$\\\midrule" + "\n")
        for mode, v in fid["table5"].items():
            out.append(f"{tex_escape(mode)} & {v['A_fid']:.3f} & {v['epsilon']:.3f} & {v['f_nonfid']:.3f}\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    for obs, d in fid["differential"].items():
        edges = d["edges"]
        out.append(r"\begin{table}[H]\centering\small\caption{Differential fiducial cross section in " + tex_escape(obs) +
                   r" (response-matrix unfolding in the likelihood).}\begin{tabular}{lccc}\toprule bin & $\sigma_{\text{fid}}$ [fb] & SM [fb]\\\midrule" + "\n")
        for r in d["fits"]["data"]:
            out.append(f"{r['bin']} ({tex_escape(str(edges[r['bin']]))}--{tex_escape(str(edges[r['bin'] + 1]))}) & "
                       f"${r['sigma_fid_fb']:.3f}^{{+{r['total_fb'][1]:.3f}}}_{{-{abs(r['total_fb'][0]):.3f}}}$ & {r['sm_fb']:.3f}\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    w, z4 = full["width"], full["z4l"]
    out.append(r"\subsection{Width and $Z\to4\ell$}" + "\n")
    out.append(f"On-shell width (1D $m_{{4\\ell}}$, 105-140~GeV, $m_H$ profiled, no interference): $\\Gamma_H<{w['observed'].get('limit_95')}$~GeV "
               f"at 95\\% CL (expected $<{w['expected'].get('limit_95')}$~GeV).  $Z\\to4\\ell$: $m_Z={z4['combined']['m_z']:.3f}\\pm{z4['combined']['stat']:.3f}"
               f"\\,(\\text{{stat}})\\pm{z4['combined']['syst']:.3f}\\,(\\text{{syst}})$~GeV "
               + "; ".join(f"{FS_TEX[fs]} {v['m_z']:.2f}" for fs, v in z4["final_states"].items()) + ".\n")
    model_dir = PRODUCTION / "inference" / args.select / args.model_cat
    for label, cap in (("scan_mu_data", "Likelihood scan of $\\mu$ (data)."), ("scan_mh_data", "Likelihood scan of $m_H$ (data)."),
                       ("scan_fv_data", "68\\% and 95\\% CL contours of $(\\mu_F,\\mu_V)$."), ("width", "Likelihood scan of $\\Gamma_H$."),
                       ("gof", "Saturated goodness of fit: pseudo-experiments and the data.")):
        sub = "width_scan.png" if label == "width" else ("gof.png" if label == "gof" else "scan.png")
        out.append(figure(fig.add(model_dir / f"{args.tag}_{label}" / "plots" / sub), cap, "0.6\\textwidth"))
    z4dir = PRODUCTION / "inference" / args.select / f"{args.tag}_z4l" / "plots"
    out.append(figures_grid([fig.add(z4dir / f"z4l_{fs}.png") for fs in FINAL_STATES], "$Z\\to4\\ell$ fits per final state.", "0.32\\textwidth"))
    return "".join(out)


def validation_section(full: dict, fig: Figures, args) -> str:
    out = [r"\section{Validation}" + "\n",
           "Closure tests on independent MC halves validate the calibration (Section~\\ref{sec:calib}) and the efficiencies "
           "(thinned MC); the likelihood is validated with the Asimov datasets (they return the injected point) and with "
           "frequentist pseudo-experiments (global observables drawn per toy).\n"]
    for spec in args.toys:
        toys = load(PRODUCTION / "inference" / spec / "toys.json")
        if not toys:
            continue
        out.append(r"\begin{table}[H]\centering\scriptsize\caption{Pseudo-experiments (" + tex_escape(spec) +
                   r"): median fitted values, pull means and widths (asymmetric MINOS errors), false discoveries.}"
                   r"\begin{tabular}{lcccccccc}\toprule point & $N$ & $\tilde\mu$ & pull$_\mu$ & width$_\mu$ & $\tilde m_H$ & pull$_{m_H}$ & width$_{m_H}$ & $Z>3$\\\midrule" + "\n")
        for key, s in toys["summary"].items():
            out.append(f"{tex_escape(key)} & {s['n_valid']} & {s['mu_median']:.3f} & {s['pull_mu_mean']:+.2f} & {s['pull_mu_width']:.2f} & "
                       f"{s['mH_median']:.2f} & {s['pull_mH_mean']:+.2f} & {s['pull_mH_width']:.2f} & {s['frac_Z_gt_3']:.3f}\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
        if toys.get("slopes"):
            out.append("Slopes of the median measured against the injected value: " +
                       ", ".join(f"{tex_escape(k)} {v:.3f}" for k, v in toys["slopes"].items()) + ".\n")
    if args.paired:
        pr = load(PRODUCTION / "inference" / args.paired / "toys.json")
        if pr:
            out.append("Paired injection (background shared, $\\Delta\\mu_{\\text{true}}=1$): " +
                       "; ".join(f"{tex_escape(k)} $\\langle\\Delta\\hat\\mu\\rangle={s['dmu_mean']:.3f}\\pm{s['dmu_mean_error']:.3f}$"
                                 for k, s in pr["summary"].items()) + ".\n")
    mc = full.get("validation", {}).get("mc_closure")
    if mc:
        out.append(r"\begin{table}[H]\centering\small\caption{Closure on the simulated events themselves (weighted to the expectation; "
                   r"truth $\mu=1$, $m_H=125$~GeV): the per-event resolution and its correlations are those of the MC, not of the model.}"
                   r"\begin{tabular}{lcc}\toprule fit & $\hat\mu$ & $\hat m_H$ [GeV]\\\midrule" + "\n")
        for k, v in mc.items():
            out.append(f"{tex_escape(k)} & {v['mu']:.4f} & {v['mH']:.4f}\\\\\n")
        out.append(r"\bottomrule\end{tabular}\end{table}" + "\n")
    if args.pyhf_check and args.pyhf_check.exists():
        chk = load(args.pyhf_check)
        out.append(f"MODEL.json (pyhf) reproduces the reported $\\mu$: {tex_escape(json.dumps(chk.get('consistency', chk))[:300])}.\n")
    return "".join(out)


def substitutions_section() -> str:
    items = [
        "Pseudo-data instead of collision data: no run or event identifiers, so no luminosity-block masks, pileup reweighting or run-dependent conditions; seeds from the file key, the original entry and the object index.",
        "Signal MC only at $m_H=125$~GeV (LO Pythia8): the $m_H$ dependence by rest-frame morphing with a morphing nuisance; no QCD-scale or PDF weights, so the theory uncertainties of set (a) are the YR4 normalization values.",
        "No ttH sample ($\\mu_{ttH}$ not measurable), no WZ or $Z\\gamma$ sample (Z+X prompt subtraction and closure with ZZ, DY and $t\\bar t$ only), ZH and WH not separable at sample level.",
        "The MC carries a production preselection (two muons or two electrons, any trigger bit): the generator tables miss about 10\\% of the events, the fiducial acceptance is relative to the preselected $\\sigma_{\\text{eff}}$; the reconstruction efficiency is not measurable with tag-and-probe (probe objects biased by the preselection).",
        "Electron identification: mvaFall17V2noIso WPL instead of the AN BDT (not in the PFNano).",
        "The lepton reconstruction efficiency is the MC's (user decision; not measurable with tag-and-probe here); the data/MC Z yield ratio of the calibration pairs is reported as a cross-check only.",
        "FSR photons from the NanoAOD FSR collection; ghost cleaning with a $\\Delta R$ criterion instead of segment arbitration.",
        "The per-event mass uncertainty from the lepton $p_T$ errors only (no full covariance); the $Z_1$ refit with Gaussian $p_T$ constraints.",
        "Calibration from $Z\\to\\ell\\ell$ only (no $J/\\psi$, $\\Upsilon$).",
        "No jet-energy-scale or resolution variations (not available): no jet systematics on the categorization; no b-tagging scale factors.",
        "The fiducial fits use the $D_{\\text{bkg}}^{\\text{kin}}$-chosen candidate at reconstruction level (the paper uses $Z_1$ closest to $m_Z$ and the largest $Z_2$ $p_T$ sum there).",
        "The width fit has no signal-background interference term (no interference sample).",
        "Signal-region $m_{Z_2}$ threshold 10~GeV and the category constants optimized on MC (not the paper's values).",
    ]
    return (r"\section{Substitutions and unavailable items}\label{sec:subst}" + "\n\\begin{itemize}\n" +
            "".join(f"\\item {i}\n" for i in items) + "\\end{itemize}\n")


def truth_section(truth: dict, result: dict, selfreport: dict) -> str:
    """The truth comparison (after the user's unblinding): the POIs and the per-flavour calibration of RESULT.json against
    the injection record and its profile definitions (truth.json of make_truth.py)."""
    rows = []
    measured = {"mu": (result["pois"][1]["value"], result["pois"][1]["total"]), "mH": (result["pois"][0]["value"], result["pois"][0]["total"])}
    for name, per in result.get("calibration", {}).items():
        for f, v in per.items():
            measured[f"{name} {f}"] = (v["value"], v["unc"])
    for k, t in truth["compare"].items():
        if k in measured:
            m, e = measured[k]
            rows.append(f"{tex_escape(k)} & {t:.5g} & {m:.5g} & {e:.3g} & {(m - t) / e:+.2f}\\\\")
    notes = "".join(f"\\item {tex_escape(n)}\n" for n in truth.get("notes", []))
    return (r"\section{Comparison with the truth (after unblinding)}" + "\n" +
            tex_escape(truth["source"]) + "\n\n" +
            r"\begin{table}[H]\centering\small\begin{tabular}{lcccc}\toprule quantity & truth & measured & uncertainty & pull\\\midrule" + "\n" +
            "\n".join(rows) + "\n" + r"\bottomrule\end{tabular}\end{table}" + "\n" +
            ("\\begin{itemize}\n" + notes + "\\end{itemize}\n" if notes else ""))


if __name__ == "__main__":
    raise SystemExit(main())
