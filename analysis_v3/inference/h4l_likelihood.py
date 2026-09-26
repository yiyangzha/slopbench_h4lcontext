"""The H -> 4l likelihood (stage 7): an extended unbinned likelihood per channel (final state x category)
in m4l (Z1-refitted or not), D_bkg^kin (templates conditional on m4l) and, for the 3D fit, the
per-event mass uncertainty (conditional signal width and templates of the relative error).

Parameters: mu (unbounded), mH (floating), and nuisances with unit-Gaussian constraints 0.5 (theta - g)^2
around the global observables g (0 for the data; the fitted nuisance values for the post-fit Asimov
dataset; drawn from N(theta_true, 1) for frequentist pseudo-experiments):
  scale_mu, scale_e   lepton momentum scale: signal peak position k mean_fs (1 + f_mu dS_mu th + f_e dS_e th)
  res_mu, res_e       lepton resolution: signal width x (1 + g_mu dR_mu th + g_e dR_e th)
  eff_mu, eff_e       efficiency scale factors: signal, qqZZ and ggZZ yields (asymmetric, from the SF payload)
  zx_<fs>             Z+X normalization per final state (log-normal, asymmetric kappa)
  morph               A x eff(m_H) morphing (interpolation between the scale and shift morphings)
  mcstat_<channel...> not included per channel (the MC statistical uncertainties are below 1 %, recorded)
  set (a) only:       lumi (2.5 %), br (2 %), qqzz_theory, ggzz_k (10 %), sig_qcd_<mode>, sig_pdf_<mode>
The signal yield of a channel is mu sum_modes M(channel, mode) Y_mode(125) sigma_eff(mode, mH)/
sigma_eff(mode, 125) R_mode(mH) (1 + eff variations), the m4l shape per mode the resonant DCB (mean
and width scaled by k = mH/125) and, for VH, the non-resonant Landau fraction.  The multiplier M
comes from the POI scheme (poi_scheme; the global mu is fixed at 1 by the fits of every scheme with
its own POIs):
  inclusive    M = 1
  final_state  mu_4mu, mu_4e, mu_2e2mu of the channel's final state
  category     mu_<category> of the channel's category
  fv           mu_F (ggH) and mu_V (VBF, VH)
  mode         mu_ggH, mu_VBF, mu_VHhad, mu_VHlep weighted by the generator-level production classes
               of the mode's yield in the channel (build_model "breakdown"; VH_unknown at the SM)
  stxs0        r_ggH, r_VBF, r_VHhad, r_VHlep for the classes with |y_H| < 2.5, the forward parts at
               the SM
  fid_fs       r_fid_4mu, r_fid_4e, r_fid_2e2mu (sigma_fid / sigma_fid^SM per generator final state):
               M = F + f_nonfid F / f_fid with F = sum over the fiducial classes of fraction x r, so the
               non-fiducial signal is a fixed fraction of the fiducial signal of the same channel
               (JHEP 11 (2017) 047 Eq. 10.2); the VH non-resonant part stays at the SM (a background)
  fid_<obs>    r_<obs>_bin<j> per generator-level bin of pt4l, njets or ptj1, the same way (the
               channels are then the reconstruction-level bins: the response matrix).
  mh_fs        m_H per final state (the paper's mutual compatibility test): mH is the 4mu mass, mH_4e and
               mH_2e2mu those of the other final states; mu common
  fid_int      the integrated fiducial cross section with the final-state fractions floating (the
               paper's integrated fit): r_fid = sigma_fid / sigma_fid^SM with r_fid_4mu and r_fid_4e free
               and r_fid_2e2mu = (S r_fid - s_4mu r_fid_4mu - s_4e r_fid_4e) / s_2e2mu, s_f the SM fiducial
               cross sections (fiducial_sm) and S their sum; otherwise as fid_fs.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path("/eos/user/y/yiyangz/codex/jfc/analyses/ref_h4l")
sys.path.insert(0, str(REPO / "analysis_v3/common/python"))
import h4l_shapes as shapes  # noqa: E402

MODES = ("ggH", "VBF", "VH")
YR4_MODES = {"ggH": "GluGluToHToZZ", "VBF": "VBF_HToZZ", "VH": "VHToZZ"}


class SigmaEff:
    """sigma_eff(mode, m_H) / sigma_eff(mode, 125) from the YR4 table (log-quadratic outside it)."""

    def __init__(self, yr4: dict):
        self.funcs = {}
        for mode in MODES:
            rows = yr4["sigma_eff"][YR4_MODES[mode]]
            masses = np.array([r["mH"] for r in rows])
            values = np.array([r["sigma_eff_pb"] for r in rows])
            ref = float(np.interp(125.0, masses, values))
            coeff = np.polyfit(masses - 125.0, np.log(values / ref), 2)
            self.funcs[mode] = (masses, values / ref, coeff)

    def ratio(self, mode: str, mh: float) -> float:
        masses, ratios, coeff = self.funcs[mode]
        if masses[0] <= mh <= masses[-1]:
            return float(np.interp(mh, masses, ratios))
        return float(np.exp(np.polyval(coeff, mh - 125.0)))


class Model:
    """Evaluates the negative log-likelihood for a dataset (events per channel, optionally weighted)."""

    def __init__(self, model: dict, signal_model: dict, yr4: dict, dimension: str = "3D", refit: bool = True,
                 result_set: str = "b", systematics: dict | None = None, poi_scheme: str = "inclusive", width: bool = False,
                 fiducial_sm: dict | None = None):
        self.model = model
        self.fiducial_sm = fiducial_sm
        self.signal_model = signal_model
        # Dimensions: 1D m4l; 2D m4l x D_bkg^kin; 2Dmass m4l x D_mass (the paper's 2D mass fit, Table 6); 3D all three.
        if dimension not in ("1D", "2D", "2Dmass", "3D"):
            raise ValueError(f"unknown dimension {dimension}")
        self.dim = dimension
        self.use_d = dimension in ("2D", "3D")
        self.use_e = dimension in ("2Dmass", "3D")
        self.variable = "m4l_refit" if refit else "m4l"
        self.error = "m4l_refit_err" if refit else "m4l_err"
        self.lo, self.hi = model["window"]
        self.sigma_eff = SigmaEff(yr4)
        self.mh_grid = np.array(signal_model["mh_grid"])
        self.m_edges = np.array(model["m_edges"])
        self.d_edges = np.array(model["d_edges"])
        self.result_set = result_set
        self.syst = systematics or {}
        self.channels = list(model["channels"])
        # Nuisance list.
        self.nuisances = ["scale_mu", "scale_e", "res_mu", "res_e", "eff_mu", "eff_e", "morph"] + \
                         [f"zx_{fs}" for fs in model["final_states"]]
        if result_set == "a":
            self.nuisances += ["lumi", "br", "qqzz_theory", "ggzz_k"] + [f"sig_qcd_{m}" for m in MODES] + [f"sig_pdf_{m}" for m in MODES]
        # The centres of the nuisance constraints (the global observables).
        self.global_observables = np.zeros(len(self.nuisances))
        self.poi_scheme = poi_scheme
        self.extra_pois, self.terms = self._poi_terms(poi_scheme)
        # The on-shell width (JHEP 11 (2017) 10.4): the resonant DCB convolved with a Breit-Wigner of width GammaH
        # (no signal-background interference: substitution); the per-event-width 3D model is not convolved.
        self.width = width
        if width:
            if self.use_e:
                raise ValueError("the width fit uses the 1D or 2D model (no per-event width)")
            self.extra_pois = self.extra_pois + ["GammaH"]
        self.fixed_nonresonant = poi_scheme.startswith("fid")
        self.parameters = ["mu", "mH"] + self.extra_pois + self.nuisances
        self.index = {name: i for i, name in enumerate(self.parameters)}
        # Resolve the POI names of the multiplier terms into parameter indices.
        self.terms = {key: ([(self.index[name] if name else None, frac) for name, frac in terms], nonfid, fidtot)
                      for key, (terms, nonfid, fidtot) in self.terms.items()}
        # Per final state: the fractions of the m4l scale and resolution responses carried by muons and electrons.
        self.flavour_fraction = self.syst.get("flavour_fraction", {"4mu": (1.0, 0.0), "4e": (0.0, 1.0), "2e2mu": (0.5, 0.5)})
        # The resolution nuisances move the m4l width by the flavour's share of the width (variance share).
        self.resolution_fraction = self.syst.get("resolution_flavour_fraction", self.flavour_fraction)

    # ------------------------------------------------------------------ POI schemes
    def _poi_terms(self, scheme: str):
        """The extra POIs of the scheme and, per (channel, mode), the multiplier terms [(POI name or None for the SM,
        fraction)], the non-fiducial fraction and the fiducial fraction (fid schemes)."""
        channels = self.model["channels"]
        prod_poi = {"ggH": "mu_ggH", "VBF": "mu_VBF", "VH_had": "mu_VHhad", "VH_lep": "mu_VHlep"}
        stxs_poi = {"ggH": "r_ggH", "VBF": "r_VBF", "VH_had": "r_VHhad", "VH_lep": "r_VHlep"}
        if scheme == "inclusive":
            extra = []
        elif scheme == "mh_fs":
            extra = [f"mH_{fs}" for fs in self.model["final_states"] if fs != "4mu"]
        elif scheme == "final_state":
            extra = [f"mu_{fs}" for fs in self.model["final_states"]]
        elif scheme == "category":
            extra = [f"mu_{c}" for c in self.model["categories"]]
        elif scheme == "fv":
            extra = ["mu_F", "mu_V"]
        elif scheme == "mode":
            extra = list(prod_poi.values())
        elif scheme == "stxs0":
            extra = list(stxs_poi.values())
        elif scheme == "fid_fs":
            extra = [f"r_fid_{fs}" for fs in self.model["final_states"]]
        elif scheme == "fid_int":
            if not self.fiducial_sm:
                raise ValueError("scheme fid_int needs the SM fiducial cross sections per final state (fiducial_sm)")
            states = list(self.model["final_states"])
            derived = states[-1]
            extra = ["r_fid"] + [f"r_fid_{fs}" for fs in states[:-1]]
            s_total = sum(self.fiducial_sm[fs] for fs in states)
        elif scheme.startswith("fid_"):
            obs = scheme[4:]
            n_bins = len(self.model["observable_edges"]) - 1 if "observable_edges" in self.model else None
            if n_bins is None or self.model.get("channel_kind") != obs:
                raise ValueError(f"scheme {scheme} needs a model built with --channels {obs}")
            extra = [f"r_{obs}_bin{k}" for k in range(n_bins)]
        else:
            raise ValueError(f"unknown POI scheme {scheme}")
        terms = {}
        for ch_name, ch in channels.items():
            for mode in MODES:
                block = ch["signal"][mode]
                bd = block.get("breakdown", {})
                nonfid, fidtot = 0.0, 1.0
                if scheme in ("inclusive", "mh_fs"):
                    t = [(None, 1.0)]
                elif scheme == "final_state":
                    t = [(f"mu_{ch['final_state']}", 1.0)]
                elif scheme == "category":
                    t = [(f"mu_{ch['category']}", 1.0)]
                elif scheme == "fv":
                    t = [("mu_F" if mode == "ggH" else "mu_V", 1.0)]
                elif scheme == "mode":
                    t = [(prod_poi.get(c), f) for c, f in bd.get("prod", {mode: 1.0}).items()]
                elif scheme == "stxs0":
                    t = []
                    for key, f in bd.get("stxs", {}).items():
                        cls, region = key.rsplit("_", 1)
                        t.append((stxs_poi.get(cls) if region == "central" else None, f))
                else:
                    table = bd.get("fid_fs" if scheme in ("fid_fs", "fid_int") else f"fid_{scheme[4:]}", {})
                    t = []
                    for key, f in table.items():
                        if key == "nonfid":
                            nonfid = f
                        elif scheme == "fid_fs":
                            t.append((f"r_{key}", f) if key[4:] in self.model["final_states"] else (None, f))
                        elif scheme == "fid_int":
                            fs = key[4:]
                            if fs == derived:
                                s_d = self.fiducial_sm[derived]
                                t.append(("r_fid", f * s_total / s_d))
                                t += [(f"r_fid_{other}", -f * self.fiducial_sm[other] / s_d) for other in states[:-1]]
                            else:
                                t.append((f"r_{key}", f) if fs in self.model["final_states"] else (None, f))
                        else:
                            t.append((f"r_{scheme[4:]}_{key}", f))
                    fidtot = sum(f for _, f in t)
                terms[(ch_name, mode)] = (t, nonfid, fidtot)
        return extra, terms

    def mh_of(self, p: np.ndarray, fs: str) -> float:
        """The Higgs boson mass of a final state (its own parameter in the mh_fs scheme, else mH)."""
        name = f"mH_{fs}"
        return float(p[self.index[name]]) if name in self.index else float(p[1])

    def default_parameters(self, mu: float = 1.0, mh: float = 125.0) -> np.ndarray:
        """The parameter vector of the SM hypothesis at (mu, m_H): the scheme's strengths and ratios at 1, the
        per-final-state masses at m_H, Gamma_H at 4.1 MeV, the nuisances at 0."""
        p = np.zeros(len(self.parameters))
        p[0], p[1] = mu, mh
        for name in self.extra_pois:
            p[self.index[name]] = mh if name.startswith("mH_") else (0.0041 if name == "GammaH" else 1.0)
        return p

    def multiplier(self, p: np.ndarray, channel: str, mode: str) -> float:
        terms, nonfid, fidtot = self.terms[(channel, mode)]
        value = sum(f * (p[i] if i is not None else 1.0) for i, f in terms)
        if nonfid > 0:
            value += nonfid * (value / fidtot if fidtot > 0 else 1.0)
        return value

    # ------------------------------------------------------------------ yields
    def acceptance_ratio(self, block: dict, mh: float, morph: float) -> float:
        scale = float(np.interp(mh, self.mh_grid, block["acceptance_ratio_scale"]))
        shift = float(np.interp(mh, self.mh_grid, block["acceptance_ratio_shift"]))
        return scale + morph * (shift - scale)

    def eff_factor(self, block: dict, th_mu: float, th_e: float) -> float:
        f = 1.0
        for th, key in ((th_mu, "sf_muon"), (th_e, "sf_electron")):
            up, down = block[key]
            f *= 1.0 + (up * th if th >= 0 else -down * th)
        return max(f, 1e-6)

    def yields(self, p: np.ndarray, channel: str) -> dict:
        ch = self.model["channels"][channel]
        mu, mh = p[0], self.mh_of(p, ch["final_state"])
        th = {n: p[self.index[n]] for n in self.nuisances}
        out = {}
        for mode in MODES:
            block = ch["signal"][mode]
            y = block["yield"] * self.sigma_eff.ratio(mode, mh) * self.acceptance_ratio(block, mh, th["morph"])
            y *= self.eff_factor(block, th["eff_mu"], th["eff_e"])
            if self.result_set == "a":
                y *= (1 + 0.025 * th["lumi"]) * (1 + 0.02 * th["br"])
                y *= (1 + self.syst.get("sig_qcd", {}).get(mode, 0.0) * th[f"sig_qcd_{mode}"])
                y *= (1 + self.syst.get("sig_pdf", {}).get(mode, 0.0) * th[f"sig_pdf_{mode}"])
            out[mode] = mu * self.multiplier(p, channel, mode) * y
            if mode == "VH":
                # The VH non-resonant part scales like the resonant one, except in the fiducial schemes (a background
                # at the SM there).
                out["VH_nonresonant_base"] = mu * y if self.fixed_nonresonant else out[mode]
        for name in ("qqZZ", "ggZZ"):
            block = ch["backgrounds"][name]
            y = block["yield"] * self.eff_factor(block, th["eff_mu"], th["eff_e"])
            if self.result_set == "a":
                y *= 1 + 0.025 * th["lumi"]
                y *= 1 + (self.syst.get("qqzz_theory", 0.0) * th["qqzz_theory"] if name == "qqZZ" else 0.10 * th["ggzz_k"])
            out[name] = y
        fs = ch["final_state"]
        zx = self.model["final_states"][fs]["zx"]
        t = th[f"zx_{fs}"]
        # Asymmetric log-normal through the combination's envelope: +1 sigma at central (1 + kappa_high), -1 sigma at
        # central (1 - kappa_low) (kappa_low = 1 - low / central of zx_estimate).
        kappa = (1.0 + zx["kappa_high"]) ** t if t >= 0 else (1.0 - zx["kappa_low"]) ** (-t)
        out["zx"] = ch["backgrounds"]["zx"]["yield"] * kappa
        return out

    # ------------------------------------------------------------------ densities
    def d_density(self, component: str, fs: str, m: np.ndarray, d: np.ndarray) -> np.ndarray:
        table = np.array(self.model["final_states"][fs]["d_templates"][component])
        mi = np.clip(np.searchsorted(self.m_edges, m, side="right") - 1, 0, len(self.m_edges) - 2)
        di = np.clip(np.searchsorted(self.d_edges, d, side="right") - 1, 0, len(self.d_edges) - 2)
        return table[mi, di]

    def e_density(self, component: str, fs: str, e: np.ndarray) -> np.ndarray:
        fse = self.model["final_states"][fs]
        edges = np.array(fse["e_edges"][self.variable])
        table = np.array(fse["e_templates"][component][self.variable])
        ei = np.clip(np.searchsorted(edges, e, side="right") - 1, 0, len(edges) - 2)
        return table[ei]

    def resonant_density(self, p: np.ndarray, fs: str, m: np.ndarray, dm: np.ndarray) -> np.ndarray:
        """The resonant DCB of the final state (shared by ggH, VBF and the resonant part of VH)."""
        sm = self.model["final_states"][fs]["signal"][self.variable]
        th = {n: p[self.index[n]] for n in ("scale_mu", "scale_e", "res_mu", "res_e")}
        f_mu, f_e = self.flavour_fraction[fs]
        g_mu, g_e = self.resolution_fraction[fs]
        d_scale = self.syst.get("scale_uncertainty", {"mu": 0.0005, "e": 0.0015})
        d_res = self.syst.get("resolution_uncertainty", {"mu": 0.10, "e": 0.10})
        mh = self.mh_of(p, fs)
        k = mh / 125.0
        scale = 1.0 + f_mu * d_scale["mu"] * th["scale_mu"] + f_e * d_scale["e"] * th["scale_e"]
        res = max(1.0 + g_mu * d_res["mu"] * th["res_mu"] + g_e * d_res["e"] * th["res_e"], 0.2)
        if self.use_e:
            # Conditional on the relative error D_mass = dm / m (the third observable): width s x D_mass x m_H,
            # a density in m normalized at fixed D_mass.
            v = sm["dcb_per_event_width"]["values"]
            sigma = v["width"] * (dm / m) * mh * res
        else:
            v = sm["dcb"]["values"]
            sigma = v["width"] * k * res
        mean = v["mean"] * k * scale
        if self.width:
            return shapes.dcb_bw_pdf(m, mean, sigma, v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"], self.lo, self.hi,
                                     p[self.index["GammaH"]])
        return shapes.dcb_pdf(m, mean, sigma, v["alpha_l"], v["n_l"], v["alpha_r"], v["n_r"], self.lo, self.hi)

    def signal_m_density(self, p: np.ndarray, fs: str, m: np.ndarray, dm: np.ndarray, mode: str) -> np.ndarray:
        """The m4l density of one production mode (the resonant DCB; VH with its non-resonant Landau)."""
        dens = self.resonant_density(p, fs, m, dm)
        if mode == "VH":
            vh = self.model["final_states"][fs]["signal"][self.variable]["vh_nonresonant"]
            dens = vh["f_res"] * dens + (1 - vh["f_res"]) * shapes.landau_pdf(m, vh["landau_mpv"], vh["landau_width"], self.lo, self.hi)
        return dens

    def prepare(self, events: dict) -> dict:
        """Precompute the parameter-independent densities of every event (backgrounds, templates)."""
        prepared = {}
        for channel in self.channels:
            fs = self.model["channels"][channel]["final_state"]
            m = events[channel]["m"]
            if len(m) == 0:
                prepared[channel] = None
                continue
            d = events[channel]["d"]
            dm = events[channel]["dm"]
            w = events[channel].get("w", np.ones(len(m)))
            # Weighted (Asimov) cells without expected content carry no information.
            keep = w > 0
            if not np.all(keep):
                m, d, dm, w = m[keep], d[keep], dm[keep], w[keep]
                if len(m) == 0:
                    prepared[channel] = None
                    continue
            e = dm / m
            fse = self.model["final_states"][fs]
            bkg = {}
            for name in ("qqZZ", "ggZZ"):
                coeff = fse["backgrounds"][name][f"bernstein_{self.variable}"]
                bkg[name] = shapes.bernstein_pdf(m, coeff, self.lo, self.hi)
            bkg["zx"] = shapes.zx_pdf(m, fse["zx"]["shape"], self.lo, self.hi)
            for name in ("qqZZ", "ggZZ", "zx"):
                if self.use_d:
                    bkg[name] = bkg[name] * self.d_density(name, fs, m, d)
                if self.use_e:
                    bkg[name] = bkg[name] * self.e_density(name, fs, e)
            sig_extra = np.ones(len(m))
            if self.use_d:
                sig_extra = sig_extra * self.d_density("signal", fs, m, d)
            if self.use_e:
                sig_extra = sig_extra * self.e_density("signal", fs, e)
            vh = fse["signal"][self.variable]["vh_nonresonant"]
            nonres = shapes.landau_pdf(m, vh["landau_mpv"], vh["landau_width"], self.lo, self.hi) * sig_extra
            prepared[channel] = {"m": m, "dm": dm, "w": w, "fs": fs, "bkg": bkg, "sig_extra": sig_extra, "nonres": nonres,
                                 "f_res": vh["f_res"]}
        return prepared

    def expected_total(self, y: dict, f_res: float) -> float:
        """The expected number of events of a channel from its yields (the VH non-resonant part may be SM-fixed)."""
        return (y["ggH"] + y["VBF"] + f_res * y["VH"] + (1 - f_res) * y["VH_nonresonant_base"] + y["qqZZ"] + y["ggZZ"] + y["zx"])

    def nll(self, p: np.ndarray, prepared: dict) -> float:
        total = 0.0
        for channel in self.channels:
            y = self.yields(p, channel)
            fs = self.model["channels"][channel]["final_state"]
            f_res = self.model["final_states"][fs]["signal"][self.variable]["vh_nonresonant"]["f_res"]
            total += self.expected_total(y, f_res)
            ev = prepared.get(channel)
            if ev is None:
                continue
            resonant_yield = y["ggH"] + y["VBF"] + ev["f_res"] * y["VH"]
            dens = resonant_yield * self.resonant_density(p, ev["fs"], ev["m"], ev["dm"]) * ev["sig_extra"]
            dens = dens + (1 - ev["f_res"]) * y["VH_nonresonant_base"] * ev["nonres"]
            for name in ("qqZZ", "ggZZ", "zx"):
                dens += y[name] * ev["bkg"][name]
            if not np.all(dens > 0):
                # A non-positive (or undefined) density at an event: the point is outside the model's domain.
                return 1e30
            total -= float(np.sum(ev["w"] * np.log(dens)))
        for i, n in enumerate(self.nuisances):
            total += 0.5 * (p[self.index[n]] - self.global_observables[i]) ** 2
        return total
