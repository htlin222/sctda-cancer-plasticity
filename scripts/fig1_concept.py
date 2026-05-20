#!/usr/bin/env python
"""
Fig 1 concept: Topology distinguishes shapes that geometry cannot.

Three panels:
  (a) A coffee mug morphs through an intermediate shape into a donut
      (homeomorphism).
  (b) Single cells lie on a noisy closed loop in state-space.
  (c) A persistence barcode summarises one robust H_1 loop above many
      short H_0 components.

Output: results/figures/fig1_concept_topology.{pdf,png}
Also mirrored to manuscript/figures/.

Run: python scripts/fig1_concept.py
"""

import logging
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "src")

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from sctda_plasticity.config import FIGURES_DIR, PROJECT_ROOT, SEED
from sctda_plasticity.visualize import save_figure, set_publication_style

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Editorial palette — kept local so this script doesn't pollute the global module
NAVY = "#1f3a68"
ORANGE = "#e98a3a"
GREY = "#9a9a9a"


def _draw_mug(ax, cx, cy, scale=1.0):
    """Stylised coffee mug centred at (cx, cy). Body 0.30 × 0.40 wide × tall."""
    s = scale
    # body — wider than tall reads as a mug, not a vase
    body = mpatches.Rectangle(
        (cx - 0.15 * s, cy - 0.20 * s), 0.30 * s, 0.40 * s,
        linewidth=1.8, edgecolor=NAVY, facecolor="white", zorder=2,
    )
    ax.add_patch(body)
    # rim ellipse at top (perspective)
    rim = mpatches.Ellipse(
        (cx, cy + 0.20 * s), 0.30 * s, 0.07 * s,
        linewidth=1.6, edgecolor=NAVY, facecolor="white", zorder=3,
    )
    ax.add_patch(rim)
    # handle: C-shape on the left side
    handle = mpatches.Arc(
        (cx - 0.15 * s, cy - 0.02 * s), 0.18 * s, 0.28 * s,
        angle=0, theta1=90, theta2=270,
        linewidth=1.8, edgecolor=NAVY, zorder=1,
    )
    ax.add_patch(handle)


def _draw_intermediate(ax, cx, cy, scale=1.0):
    """Mug-into-donut midstage: a thick-walled bowl seen in perspective,
    with the handle still attached. The handle and cavity are drawn
    explicitly so genus = 1 is visually preserved throughout the
    deformation (mug → bowl-with-handle → donut).
    """
    s = scale
    # Bowl body (wider and flatter than the mug, suggesting it has
    # been "squashed" toward becoming a torus). Filled navy with a
    # darker rim outline.
    body = mpatches.Ellipse(
        (cx, cy - 0.02 * s), 0.46 * s, 0.34 * s,
        linewidth=1.8, edgecolor=NAVY, facecolor=NAVY, alpha=0.85, zorder=2,
    )
    ax.add_patch(body)
    # Cavity / opening — a smaller upper ellipse in white showing the
    # bowl is concave (not a solid disc). Off-centre to the right so
    # the handle on the left reads as separate from the cavity.
    cavity = mpatches.Ellipse(
        (cx + 0.03 * s, cy + 0.05 * s), 0.22 * s, 0.10 * s,
        linewidth=1.2, edgecolor="white", facecolor="white", zorder=3,
    )
    ax.add_patch(cavity)
    # Handle preserved on the left side — C-shape arc, smaller than the
    # mug's but unambiguously a topological handle (one hole through it).
    handle = mpatches.Arc(
        (cx - 0.23 * s, cy - 0.02 * s), 0.14 * s, 0.22 * s,
        angle=0, theta1=90, theta2=270,
        linewidth=1.8, edgecolor=NAVY, zorder=4,
    )
    ax.add_patch(handle)


def _draw_donut(ax, cx, cy, scale=1.0):
    """Donut (annulus)."""
    s = scale
    outer = mpatches.Circle(
        (cx, cy), 0.24 * s,
        linewidth=1.6, edgecolor=NAVY, facecolor=NAVY, alpha=0.85, zorder=2,
    )
    ax.add_patch(outer)
    hole = mpatches.Circle(
        (cx, cy), 0.10 * s,
        linewidth=0, facecolor="white", zorder=3,
    )
    ax.add_patch(hole)


def _draw_arrow(ax, x0, x1, y):
    ax.annotate(
        "", xy=(x1, y), xytext=(x0, y),
        arrowprops=dict(arrowstyle="->", lw=0.8, color=GREY),
    )


def panel_mug_to_donut(ax):
    # Place three shapes side-by-side. Use data coords matching the panel
    # aspect (≈1.55:1 after label-space allowance) so equal-aspect doesn't
    # crush the layout into a central square.
    _draw_mug(ax, cx=0.20, cy=0.5)
    _draw_arrow(ax, 0.36, 0.51, 0.5)
    _draw_intermediate(ax, cx=0.72, cy=0.5)
    _draw_arrow(ax, 0.96, 1.11, 0.5)
    _draw_donut(ax, cx=1.30, cy=0.5)
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_loop_cells(ax, rng):
    """60 cells on a noisy unit circle, coloured by angular position."""
    n = 60
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r = 1.0 + rng.normal(0, 0.07, n)
    t = theta + rng.normal(0, 0.07, n)
    x = r * np.cos(t)
    y = r * np.sin(t)
    cmap = plt.get_cmap("twilight_shifted")
    colors = cmap(theta / (2 * np.pi))
    ax.scatter(x, y, c=colors, s=20, edgecolors="white", linewidths=0.4, zorder=3)
    # discreet axis lines as a visual guide
    ax.axhline(0, color=GREY, lw=0.4, zorder=1, alpha=0.6)
    ax.axvline(0, color=GREY, lw=0.4, zorder=1, alpha=0.6)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_barcode(ax, rng):
    """Eight short H_0 bars + one long H_1 bar (orange)."""
    n_short = 8
    starts = rng.uniform(0.04, 0.18, n_short)
    lengths = rng.uniform(0.10, 0.30, n_short)
    order = np.argsort(lengths)[::-1]  # longest at top of the H0 stack
    bar_dy = 0.075
    for i, idx in enumerate(order):
        y = 0.10 + i * bar_dy
        ax.plot(
            [starts[idx], starts[idx] + lengths[idx]], [y, y],
            color=NAVY, lw=2.4, solid_capstyle="round",
        )
    # one robust H_1 bar, sitting close above the H_0 stack
    y_long = 0.10 + n_short * bar_dy + 0.05
    ax.plot([0.06, 0.96], [y_long, y_long],
            color=ORANGE, lw=3.6, solid_capstyle="round")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, y_long + 0.10)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    set_publication_style()
    rng = np.random.default_rng(SEED)

    # Nature double-column ≈ 180 mm wide; concept figure is short and wide
    fig = plt.figure(figsize=(180 / 25.4, 55 / 25.4))
    gs = gridspec.GridSpec(
        1, 3, figure=fig,
        width_ratios=[1.6, 1.0, 1.1],
        wspace=0.18,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    panel_mug_to_donut(ax_a)
    ax_a.text(-0.02, 1.02, "a", transform=ax_a.transAxes,
              fontsize=10, fontweight="bold", va="top")

    ax_b = fig.add_subplot(gs[0, 1])
    panel_loop_cells(ax_b, rng)
    ax_b.text(-0.02, 1.02, "b", transform=ax_b.transAxes,
              fontsize=10, fontweight="bold", va="top")

    ax_c = fig.add_subplot(gs[0, 2])
    panel_barcode(ax_c, rng)
    ax_c.text(-0.02, 1.02, "c", transform=ax_c.transAxes,
              fontsize=10, fontweight="bold", va="top")

    save_figure(fig, "fig1_concept_topology", output_dir=str(FIGURES_DIR))

    manuscript_dir = PROJECT_ROOT / "manuscript" / "figures"
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        src = FIGURES_DIR / f"fig1_concept_topology.{ext}"
        dst = manuscript_dir / f"fig1_concept_topology.{ext}"
        shutil.copyfile(src, dst)
        logger.info(f"Mirrored: {dst}")

    plt.close(fig)


if __name__ == "__main__":
    main()
