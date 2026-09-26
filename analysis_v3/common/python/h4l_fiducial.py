"""Fiducial and simplified-template-cross-section definitions shared by the generator-level
cross sections (signal_model/scripts/fiducial_xsec.py), the model builder and the fits.

Generator level (the skim GenTable, JHEP 11 (2017) 047 Table 4 with dressed leptons, and the
h4l_select SR rows' g_* copies): fid_pass, fid_final_state (0 4mu, 1 4e, 2 2e2mu), fid_pt4l,
fid_njets and fid_jet1_pt (jets pT > 30 GeV, |eta| < 2.5, dR > 0.4 from the leptons; -1 = no
jet), stage0 (STXS stage 0: 10/11 ggH forward/central, 20/21 VBF, 22/23 VH hadronic, 30/31 WH
leptonic, 40/41 ZH leptonic), vh_class (1 W->lnu, 2 Z->ll, 3 Z->nunu, 4 W->qq, 5 Z->qq).
Reconstruction level: the same observables of the selected candidate (pt4l; jets of the
candidate record with pT > 30 GeV and |eta| < 2.5).

Differential bins (JHEP 11 (2017) 047 Fig. 10, the binning of CMS-HIG-16-041): pT(H) 0-15-30-
85-200-inf GeV; N(jets) 0, 1, 2, >= 3; pT of the leading jet: no jet (pT < 30 GeV), 30-55,
55-95, 95-200, > 200 GeV.  Internal edges are half-open [lo, hi), the last bin is unbounded.
"""

from __future__ import annotations

import numpy as np

FINAL_STATES = {0: "4mu", 1: "4e", 2: "2e2mu"}
DIFFERENTIAL = {"pt4l": [0.0, 15.0, 30.0, 85.0, 200.0, np.inf],
                "njets": [0.0, 1.0, 2.0, 3.0, np.inf],
                "ptj1": [0.0, 30.0, 55.0, 95.0, 200.0, np.inf]}
CENTRAL_JET = {"pt": 30.0, "abs_eta": 2.5}
STXS_ABS_Y = 2.5
# Production classes of the paper's per-mode signal strengths and stage-0 measurement.
PRODUCTION_CLASSES = ("ggH", "VBF", "VH_had", "VH_lep")


def bin_of(values: np.ndarray, edges: list[float]) -> np.ndarray:
    """Index of the half-open bin [lo, hi) of every value (values below the first edge go to bin 0)."""
    return np.clip(np.searchsorted(np.asarray(edges), values, side="right") - 1, 0, len(edges) - 2)


def production_class(mode: str, vh_class: np.ndarray) -> np.ndarray:
    """ggH, VBF, VH_had (W/Z -> qq), VH_lep (W -> l nu, Z -> ll, Z -> nu nu); VH events without a decay
    class are VH_unknown."""
    n = len(vh_class)
    if mode in ("ggH", "VBF"):
        return np.full(n, mode, dtype=object)
    out = np.full(n, "VH_unknown", dtype=object)
    out[np.isin(vh_class, (4, 5))] = "VH_had"
    out[np.isin(vh_class, (1, 2, 3))] = "VH_lep"
    return out


def gen_observable(name: str, rows: dict, prefix: str = "") -> np.ndarray:
    """The generator-level differential observable (GenTable names with prefix "", SR rows with "g_")."""
    if name == "pt4l":
        return np.asarray(rows[prefix + "fid_pt4l"], dtype=np.float64)
    if name == "njets":
        return np.asarray(rows[prefix + "fid_njets"], dtype=np.float64)
    if name == "ptj1":
        return np.maximum(np.asarray(rows[prefix + "fid_jet1_pt"], dtype=np.float64), 0.0)
    raise KeyError(name)


def reco_observable(name: str, rows: dict) -> np.ndarray:
    """The reconstruction-level differential observable of selected candidates (rows need pt4l, jet_pt and
    jet_eta)."""
    if name == "pt4l":
        return np.asarray(rows["pt4l"], dtype=np.float64)
    # h4l_select stores the four leading analysis jets (pT > 30 GeV, |eta| < 4.7; -99 = none): the central jets
    # are counted among them (an event with more than four jets may miss a central one: negligible).
    pt = np.asarray(rows["jet_pt"], dtype=np.float64).reshape(len(rows["pt4l"]), -1)
    eta = np.asarray(rows["jet_eta"], dtype=np.float64).reshape(pt.shape)
    central = (pt > CENTRAL_JET["pt"]) & (np.abs(eta) < CENTRAL_JET["abs_eta"])
    if name == "njets":
        return central.sum(axis=1).astype(np.float64)
    if name == "ptj1":
        return np.where(central, pt, 0.0).max(axis=1)
    raise KeyError(name)
