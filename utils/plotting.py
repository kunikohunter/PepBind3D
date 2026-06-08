"""Shared matplotlib style for validation figures."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt


def set_plot_style() -> None:
    """Apply a clean, publication-leaning matplotlib style."""
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "lines.linewidth": 1.5,
        }
    )
