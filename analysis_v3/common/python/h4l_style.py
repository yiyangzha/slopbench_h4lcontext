"""Shared plot style of the v3 analysis: mplhep CMS style, PDF and PNG output.

Every figure carries the label "CMS Private work" (the data are called
pseudo-data), the luminosity where data enter, and is written as PDF and PNG
(dpi 200, tight bounding box).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402

hep.style.use("CMS")


def label(axis, data: bool, lumi_fb: float | None = None) -> None:
    """CMS label; data=False marks a simulation-only figure."""
    hep.cms.label(llabel="Private work" if data else "Simulation Private work", data=True, lumi=lumi_fb,
                  com=13, ax=axis)


def figure(ratio: bool = False):
    if ratio:
        fig, (main, lower) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                          gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0})
        return fig, (main, lower)
    fig, axis = plt.subplots(figsize=(10, 10))
    return fig, axis


def save(fig, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in (".pdf", ".png"):
        path = stem.with_suffix(suffix)
        fig.savefig(path, bbox_inches="tight", dpi=200)
        written.append(path)
    plt.close(fig)
    return written
