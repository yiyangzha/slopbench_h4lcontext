"""Event categories of AN-16-442 v8 section 6.1 (configuration categories_ul16_v<N>.json).

category_index(rows, config) returns, per selected candidate row, the index into
config["order"] (0 VBF-2jet, 1 VH-hadronic, 2 VH-leptonic, 3 ttH, 4 VH-MET, 5
VBF-1jet, 6 Untagged): an event goes to the first category whose requirements it
meets.  The rows need the branches CATEGORY_BRANCHES (h4l_select outputs; the stored
discriminants and the raw probabilities, so that either kind of configuration applies).

Production discriminants:
  * a configuration without a "discriminants" block (cat_v1) uses the stored
    MELA-normalized d_2jet, d_1jet, d_wh_mela, d_zh_mela (-1 = undefined for the
    jet multiplicity);
  * a configuration with it (cat_v2, user decision 2026-09-25 after the stage-4b MC
    optimization) recomputes them from the stored raw MELA probabilities with the
    paper convention D' = 1 / (1 + c(m4l) P_bkg / P_sig), c(m4l) = c0 (m4l /
    m_ref)^p:
        D_2jet: P_sig = p_JJVBF,                                      P_bkg = p_JJQCD (njets >= 2)
        D_1jet: P_sig = p_JVBF pAux_JVBF,                             P_bkg = p_JQCD  (njets == 1)
        D_WH:   P_sig = p_HadWH p_HadWH_mavjj / p_HadWH_mavjj_true,   P_bkg = p_JJQCD (njets >= 2)
        D_ZH:   P_sig = p_HadZH p_HadZH_mavjj / p_HadZH_mavjj_true,   P_bkg = p_JJQCD (njets >= 2)
    (-1 where undefined: wrong jet multiplicity or a non-positive probability).
The D_bkg^kin of the likelihood is always the stored one.
"""

from __future__ import annotations

import numpy as np

BASE_BRANCHES = ["njets", "nbjets", "n_extra_leptons", "n_extra_os_sf_pairs", "met", "m4l"]
STORED_DISCRIMINANTS = ["d_2jet", "d_1jet", "d_wh_mela", "d_zh_mela"]
RAW_PROBABILITIES = ["p_JJVBF", "p_JJQCD", "p_JVBF", "pAux_JVBF", "p_JQCD", "p_HadWH", "p_HadWH_mavjj", "p_HadWH_mavjj_true",
                     "p_HadZH", "p_HadZH_mavjj", "p_HadZH_mavjj_true"]
CATEGORY_BRANCHES = BASE_BRANCHES + STORED_DISCRIMINANTS + RAW_PROBABILITIES
# Discriminant name -> (row key used by category_index, required jet multiplicity: 2 = njets >= 2, 1 = njets == 1).
PRODUCTION = {"D_2jet": ("d_2jet", 2), "D_1jet": ("d_1jet", 1), "D_WH": ("d_wh_mela", 2), "D_ZH": ("d_zh_mela", 2)}


def _probabilities(rows: dict, name: str) -> tuple[np.ndarray, np.ndarray]:
    f8 = lambda key: np.asarray(rows[key], dtype=np.float64)  # noqa: E731

    def ratio(numerator, denominator):
        return np.where(denominator > 0, numerator / np.where(denominator > 0, denominator, 1.0), 0.0)

    if name == "D_2jet":
        return f8("p_JJVBF"), f8("p_JJQCD")
    if name == "D_1jet":
        return f8("p_JVBF") * f8("pAux_JVBF"), f8("p_JQCD")
    if name == "D_WH":
        return f8("p_HadWH") * ratio(f8("p_HadWH_mavjj"), f8("p_HadWH_mavjj_true")), f8("p_JJQCD")
    if name == "D_ZH":
        return f8("p_HadZH") * ratio(f8("p_HadZH_mavjj"), f8("p_HadZH_mavjj_true")), f8("p_JJQCD")
    raise KeyError(name)


def production_discriminants(rows: dict, config: dict) -> dict:
    """The four production discriminants D' of the configuration, under the row keys of category_index."""
    out = {}
    nj = np.asarray(rows["njets"])
    m4l = np.asarray(rows["m4l"], dtype=np.float64)
    for name, (key, multiplicity) in PRODUCTION.items():
        spec = config["discriminants"][name]
        p_sig, p_bkg = _probabilities(rows, name)
        defined = (p_sig > 0) & (p_bkg > 0) & np.isfinite(p_sig) & np.isfinite(p_bkg)
        defined &= (nj >= 2) if multiplicity == 2 else (nj == 1)
        c = spec["c0"] * (m4l / spec.get("m_ref", 125.0)) ** spec.get("p", 0.0)
        d = np.full(len(nj), -1.0)
        d[defined] = 1.0 / (1.0 + c[defined] * p_bkg[defined] / p_sig[defined])
        out[key] = d
    return out


def d2jet_threshold(m4l: np.ndarray, wp) -> np.ndarray:
    if isinstance(wp, dict):
        return wp["a"] - wp["b"] / (m4l + wp["c"])
    return np.full(len(m4l), float(wp))


def category_index(rows: dict, config: dict) -> np.ndarray:
    wp = config["working_points"]
    d = production_discriminants(rows, config) if "discriminants" in config else {k: rows[k] for k in STORED_DISCRIMINANTS}
    nj, nb = rows["njets"], rows["nbjets"]
    extra, os_pairs = rows["n_extra_leptons"], rows["n_extra_os_sf_pairs"]
    four = extra == 0
    two_three = (nj >= 2) & (nj <= 3)
    four_plus_nob = (nj >= 4) & (nb == 0)
    vbf2 = four & ((two_three & (nb <= 1)) | four_plus_nob) & (d["d_2jet"] > d2jet_threshold(rows["m4l"], wp["d_2jet"]))
    vhhad = four & (two_three | four_plus_nob) & ((d["d_wh_mela"] > wp["d_wh"]) | (d["d_zh_mela"] > wp["d_zh"]))
    vhlep = ((nj <= 3) & (nb == 0) & ((extra == 1) | ((extra == 2) & (os_pairs >= 1)))) | ((nj == 0) & (extra >= 1))
    tth = ((nj >= 4) & (nb >= 1)) | (extra >= 1)
    vhmet = four & (nj <= 1) & (rows["met"] > wp["met_vh"])
    vbf1 = four & (nj == 1) & (d["d_1jet"] > wp["d_1jet"])
    out = np.full(len(nj), 6, dtype=np.int8)
    assigned = np.zeros(len(nj), dtype=bool)
    for index, mask in enumerate((vbf2, vhhad, vhlep, tth, vhmet, vbf1)):
        take = mask & ~assigned
        out[take] = index
        assigned |= take
    return out
