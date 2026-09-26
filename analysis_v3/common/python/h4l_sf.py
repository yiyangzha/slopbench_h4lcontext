"""Lepton efficiency scale factors of the tag-and-probe payload (production_v3/tnp/<extract>/<run>/<label>/sf.json).

The per-lepton SF is the product of the chain id | loose, sip | id, iso | id and sip (the "product" rows
of sf.json) in the lepton's (pT, |eta|) bin (|eta_SC| for electrons); bins are defined by lower edges,
the internal boundaries are half-open and the terminal bins unbounded (payload convention).  A bin
without a product value falls back to the directly measured full step, else to 1 with the error of the
inclusive product (recorded).  The event SF is the product over the four leptons.  The uncertainty has
two parts (final review 2026-09-25): the fit-model part moves every SF of one flavour coherently
(sf_<flavour>_model_up/down), the statistical part is independent from bin to bin (per lepton its relative
statistical error and its global bin, relstat_<flavour> and bin_<flavour>, for the quadrature over bins).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FLAVOURS = {13: "muon", 11: "electron"}


class LeptonSF:
    def __init__(self, path: str | Path):
        self.path = str(path)
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        self.tables = {}
        self.fallbacks = {}
        for pdg, name in FLAVOURS.items():
            flavour = payload["flavours"][name]
            bins = flavour["bins"]
            pt_edges, eta_edges = np.array(bins["pt_edges"], dtype=float), np.array(bins["eta_edges"], dtype=float)
            n_pt, n_eta = len(pt_edges), len(eta_edges)
            value = np.ones(n_pt * n_eta)
            error = np.zeros(n_pt * n_eta)
            stat = np.zeros(n_pt * n_eta)
            model = np.zeros(n_pt * n_eta)
            product = flavour["product"]
            inclusive = next(r for r in product if r["bin"] == "all")
            fallback = []
            full = {r["bin"]: r.get("sf") for r in flavour["steps"]["full"]}
            for row in product:
                if row["bin"] == "all":
                    continue
                k = int(row["bin"])
                if row.get("value") is not None:
                    value[k], error[k], stat[k], model[k] = row["value"], row["total"], row["stat"], row["fit_model"]
                elif full.get(row["bin"]):
                    f = full[row["bin"]]
                    value[k], error[k] = f["value"], float(np.hypot(f["stat"], f["fit_model"]))
                    stat[k], model[k] = f["stat"], f["fit_model"]
                    fallback.append({"bin": k, "source": "full step"})
                else:
                    value[k], error[k] = 1.0, inclusive.get("total", 0.0)
                    stat[k], model[k] = inclusive.get("stat", 0.0), inclusive.get("fit_model", 0.0)
                    fallback.append({"bin": k, "source": "unity with the inclusive error"})
            self.tables[pdg] = (pt_edges, eta_edges, value, error)
            self.parts = getattr(self, "parts", {})
            self.parts[pdg] = (stat, model)
            self.fallbacks[name] = fallback

    def lepton(self, pdg, pt, eta, eta_sc):
        """Per-lepton SF and total error (arrays of the same shape as pt)."""
        pdg, pt = np.asarray(pdg), np.asarray(pt, dtype=float)
        value = np.ones(pt.shape)
        error = np.zeros(pt.shape)
        for flavour, (pt_edges, eta_edges, v, e) in self.tables.items():
            sel = np.abs(pdg) == flavour
            if not np.any(sel):
                continue
            abs_eta = np.abs(np.asarray(eta)[sel] if flavour == 13 else np.asarray(eta_sc)[sel])
            ip = np.clip(np.searchsorted(pt_edges, pt[sel], side="right") - 1, 0, len(pt_edges) - 1)
            ie = np.clip(np.searchsorted(eta_edges, abs_eta, side="right") - 1, 0, len(eta_edges) - 1)
            k = ip * len(eta_edges) + ie
            value[sel], error[sel] = v[k], e[k]
        return value, error

    def bins(self, pdg, pt, eta, eta_sc):
        """Per lepton: its bin index within its flavour's table (-1 for leptons of no table)."""
        pdg, pt = np.asarray(pdg), np.asarray(pt, dtype=float)
        out = np.full(pt.shape, -1, dtype=np.int64)
        for flavour, (pt_edges, eta_edges, _, _) in self.tables.items():
            sel = np.abs(pdg) == flavour
            if not np.any(sel):
                continue
            abs_eta = np.abs(np.asarray(eta)[sel] if flavour == 13 else np.asarray(eta_sc)[sel])
            ip = np.clip(np.searchsorted(pt_edges, pt[sel], side="right") - 1, 0, len(pt_edges) - 1)
            ie = np.clip(np.searchsorted(eta_edges, abs_eta, side="right") - 1, 0, len(eta_edges) - 1)
            out[sel] = ip * len(eta_edges) + ie
        return out

    def event(self, rows: dict) -> dict:
        """Event SF (product over the four leptons), the per-flavour +-1 sigma variations of the total error (kept for
        reference), the coherent fit-model variations and the per-lepton statistical parts."""
        value, error = self.lepton(rows["l_pdg"], rows["l_pt"], rows["l_eta"], rows["l_eta_sc"])
        bins = self.bins(rows["l_pdg"], rows["l_pt"], rows["l_eta"], rows["l_eta_sc"])
        out = {"sf": np.prod(value, axis=1)}
        for flavour, name in FLAVOURS.items():
            is_flavour = np.abs(rows["l_pdg"]) == flavour
            stat_table, model_table = self.parts[flavour]
            k = np.where(is_flavour, bins, 0)
            model = np.where(is_flavour, model_table[k], 0.0)
            for sign, label in ((1, "up"), (-1, "down")):
                out[f"sf_{name}_{label}"] = np.prod(np.where(is_flavour, value + sign * error, value), axis=1)
                out[f"sf_{name}_model_{label}"] = np.prod(np.where(is_flavour, value + sign * model, value), axis=1)
            out[f"relstat_{name}"] = np.where(is_flavour, stat_table[k] / value, 0.0)
            out[f"bin_{name}"] = np.where(is_flavour, bins, -1)
        return out


def efficiency_variation(w: np.ndarray, rows: dict, sel: np.ndarray, name: str) -> tuple:
    """Relative (up, down) yield variation of one flavour's SF uncertainty for the selected rows: the coherent
    fit-model shift in quadrature with the statistical part summed in quadrature over the independent bins (each bin
    moving all its leptons together; linear in the relative errors)."""
    sf = rows["sf"][sel]
    y = float(np.sum(w * sf))
    if y <= 0:
        return 0.0, 0.0, {"model": [0.0, 0.0], "stat": 0.0}
    model_up = float(np.sum(w * rows[f"sf_{name}_model_up"][sel])) / y - 1.0
    model_down = float(np.sum(w * rows[f"sf_{name}_model_down"][sel])) / y - 1.0
    bins, rel = rows[f"bin_{name}"][sel], rows[f"relstat_{name}"][sel]
    per_lepton = (w * sf)[:, None] * rel
    use = bins >= 0
    shifts = np.bincount(bins[use], weights=per_lepton[use]) if np.any(use) else np.zeros(1)
    stat = float(np.sqrt(np.sum(shifts ** 2))) / y
    return (float(np.hypot(model_up, stat)), -float(np.hypot(model_down, stat)),
            {"model": [model_up, model_down], "stat": stat})
