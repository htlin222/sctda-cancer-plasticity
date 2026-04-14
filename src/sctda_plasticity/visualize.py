"""
Publication-quality visualization for TDA results.
"""

import logging
from pathlib import Path
from typing import Dict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Publication style settings
# ──────────────────────────────────────────────
NATURE_COMMS_STYLE = {
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "lines.linewidth": 0.8,
}

TIMEPOINT_COLORS = {
    "D0": "#2196F3",
    "D1": "#4CAF50",
    "D2": "#FF9800",
    "D4": "#FF5722",
    "D9": "#F44336",
    "D11": "#9C27B0",
}

PH_COLORS = {
    0: "#333333",  # H₀
    1: "#E91E63",  # H₁
    2: "#00BCD4",  # H₂
}


def set_publication_style():
    """Set matplotlib rcParams for publication-quality figures."""
    mpl.rcParams.update(NATURE_COMMS_STYLE)
    sns.set_style("ticks")


def save_figure(fig, name: str, output_dir: str = "results/figures", formats=("pdf", "png")):
    """Save figure in multiple formats."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        filepath = output_path / f"{name}.{fmt}"
        fig.savefig(filepath, format=fmt)
        logger.info(f"Saved: {filepath}")


def plot_persistence_barcodes(
    dgms_dict: Dict[str, list],
    max_dim: int = 1,
    title: str = "",
    figsize: tuple = (180 / 25.4, 60 / 25.4),
) -> plt.Figure:
    """Plot persistence barcodes for multiple time points side by side.

    Parameters
    ----------
    dgms_dict : dict
        {timepoint_name: [dgm_H0, dgm_H1, ...]}
    max_dim : int
        Maximum dimension to plot
    title : str
        Figure title
    figsize : tuple
        Figure size in inches (default: Nature Comms single column)
    """
    set_publication_style()
    n_tp = len(dgms_dict)
    fig, axes = plt.subplots(max_dim + 1, n_tp, figsize=figsize, squeeze=False)

    for col, (tp_name, dgms) in enumerate(dgms_dict.items()):
        for dim in range(max_dim + 1):
            ax = axes[dim, col]
            dgm = dgms[dim]

            # Filter finite bars
            finite = dgm[np.isfinite(dgm[:, 1])] if len(dgm) > 0 else dgm

            if len(finite) > 0:
                persistence = finite[:, 1] - finite[:, 0]
                sorted_idx = np.argsort(persistence)[::-1]
                for i, idx in enumerate(sorted_idx[:30]):  # top 30 bars
                    ax.barh(
                        i, persistence[idx],
                        left=finite[idx, 0],
                        height=0.8,
                        color=PH_COLORS.get(dim, "#666"),
                        alpha=0.7,
                    )

            ax.set_xlabel("Filtration" if dim == max_dim else "")
            if col == 0:
                ax.set_ylabel(f"$H_{dim}$")
            if dim == 0:
                ax.set_title(tp_name, color=TIMEPOINT_COLORS.get(tp_name, "#333"))
            ax.set_yticks([])

    fig.suptitle(title, fontsize=9, y=1.02)
    fig.tight_layout()
    return fig


def plot_persistence_comparison(
    dgms_dict: Dict[str, np.ndarray],
    dim: int = 1,
    figsize: tuple = (85 / 25.4, 75 / 25.4),
) -> plt.Figure:
    """Plot persistence diagrams from multiple conditions overlaid.

    Parameters
    ----------
    dgms_dict : dict
        {condition_name: dgm_array}
    dim : int
        Homology dimension
    """
    set_publication_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Auto-assign colors: check TIMEPOINT_COLORS first, then tab10 palette
    cmap = plt.get_cmap("tab10")
    all_max = 0
    for i, (name, dgm) in enumerate(dgms_dict.items()):
        finite = dgm[np.isfinite(dgm[:, 1])] if len(dgm) > 0 else dgm
        if len(finite) > 0:
            # Try direct match, then strip prefix (e.g., "erl_D0" → "D0")
            color = TIMEPOINT_COLORS.get(name)
            if color is None and "_" in name:
                stripped = name.rsplit("_", 1)[-1]
                color = TIMEPOINT_COLORS.get(stripped)
            if color is None:
                color = cmap(i % 10)
            ax.scatter(
                finite[:, 0], finite[:, 1],
                c=[color], label=name, s=15, alpha=0.7,
                edgecolors="white", linewidth=0.3,
            )
            all_max = max(all_max, finite.max())

    # Diagonal
    lims = [0, all_max * 1.1]
    ax.plot(lims, lims, "k--", lw=0.5, alpha=0.3)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Birth")
    ax.set_ylabel("Death")
    ax.set_title(f"$H_{dim}$ Persistence Diagram")
    ax.legend(frameon=False, loc="lower right")
    ax.set_aspect("equal")
    fig.tight_layout()
    return fig


def plot_permutation_test(
    null_distribution: np.ndarray,
    observed: float,
    p_value: float,
    figsize: tuple = (85 / 25.4, 60 / 25.4),
) -> plt.Figure:
    """Plot permutation test result."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=figsize)

    ax.hist(null_distribution, bins=30, color="#ccc", edgecolor="#999", alpha=0.8, label="Null")
    ax.axvline(observed, color="#E91E63", lw=1.5, label=f"Observed (p={p_value:.4f})")
    ax.set_xlabel(r"Max $H_1$ persistence")
    ax.set_ylabel("Count")
    ax.set_title(r"Permutation Test for $H_1$")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_wasserstein_heatmap(
    distances: np.ndarray,
    labels: list,
    figsize: tuple = (85 / 25.4, 75 / 25.4),
) -> plt.Figure:
    """Plot Wasserstein distance heatmap between persistence diagrams."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(distances, cmap="YlOrRd", aspect="equal")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    plt.colorbar(im, ax=ax, label="Wasserstein distance", shrink=0.8)
    ax.set_title(r"$H_1$ Diagram Distances")
    fig.tight_layout()
    return fig
