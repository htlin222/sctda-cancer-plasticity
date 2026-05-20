#!/usr/bin/env python
"""
Fig 1 — Graphical abstract: drug-tolerant cancer cells cycle between
states, measurable as a rising max H_1 persistence across cohorts.

Three panels tell the biology-first story:
  (a) AI render: untreated cancer cells (uniform cluster) → drug arrow
      → cells cycling between states (loop of differently-coloured
      cells). The biological finding.
  (b) AI render: same cycling visible in scRNA-seq as a closed loop in
      state-space (kept from previous version).
  (c) Bar chart: paired untreated-vs-treated max H_1 across 4 cohorts
      (PC9, PDX YU-006, PDX YU-003, Maynard patient cohort), with
      fold-change labels. The headline result.

The mug-donut topological-equivalence pedagogy now lives in one
sentence in the Background prose, not Figure 1.

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
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from sctda_plasticity.config import FIGURES_DIR, PROJECT_ROOT, SEED
from sctda_plasticity.visualize import save_figure, set_publication_style

PANEL_A_PNG = FIGURES_DIR / "fig1_panel_a_biology.png"
PANEL_B_PNG = FIGURES_DIR / "fig1_panel_b_3d.png"

# Headline paired (untreated, treated) max H_1 values across cohorts.
COHORT_DATA = [
    ("PC9\nosi D0→D14", 1.41, 3.78, r"2.7$\times$"),
    ("PDX\nYU-006", 2.81, 6.08, r"2.2$\times$"),
    ("PDX\nYU-003", 2.79, 3.85, r"1.4$\times$"),
    ("Maynard\nTN→PD", 5.37, 8.05, r"1.5$\times$"),
]


def _imshow_alpha_panel(ax, png_path, subtitle, overlay_fn=None):
    """Inset an alpha PNG with italic subtitle. `overlay_fn(ax, w, h)` runs
    after imshow, in image-pixel coordinates, to add labels/arrows."""
    if not png_path.exists():
        ax.text(0.5, 0.5, f"Missing: {png_path}",
                ha="center", va="center", color="red",
                transform=ax.transAxes)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    img = mpimg.imread(str(png_path))
    ax.imshow(img, interpolation="bilinear", aspect="equal")
    h, w = img.shape[0], img.shape[1]
    if overlay_fn is not None:
        overlay_fn(ax, w, h)
    ax.text(
        w / 2, h * 1.10, subtitle,
        fontsize=7.5, color=NAVY, style="italic",
        ha="center", va="top",
    )
    # extra right-side room for overlay callouts (mostly used by b and c)
    ax.set_xlim(-w * 0.02, w * 1.55)
    ax.set_ylim(h * 1.28, -h * 0.06)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _overlay_panel_b(ax, w, h):
    """Add explanatory callouts on top of the cells-on-loop image."""
    # Callout to one sphere on the loop, labeling it as one cell.
    ax.annotate(
        "one dot =\none cell",
        xy=(w * 0.78, h * 0.22), xytext=(w * 1.02, h * 0.10),
        fontsize=7, color=NAVY, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", lw=0.6, color=NAVY,
                        shrinkA=0, shrinkB=2),
    )


def _draw_filtration_step(ax, cx, cy, loop_r, ball_r, n_cells=10):
    """Draw n_cells on a loop centred at (cx, cy) with filtration balls
    of radius ball_r drawn around each cell."""
    theta = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)
    cell_x = cx + loop_r * np.cos(theta)
    cell_y = cy + loop_r * np.sin(theta)
    for x, y in zip(cell_x, cell_y):
        # filled ball (subtle)
        ax.add_patch(mpatches.Circle(
            (x, y), ball_r,
            facecolor=NAVY, edgecolor="none", alpha=0.18, zorder=2,
        ))
        # ball outline
        ax.add_patch(mpatches.Circle(
            (x, y), ball_r,
            facecolor="none", edgecolor=NAVY, lw=0.5, alpha=0.55, zorder=3,
        ))
    # cell dots on top
    ax.scatter(cell_x, cell_y, c=[NAVY], s=14, zorder=10,
               edgecolors="white", linewidths=0.4)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Editorial palette
NAVY = "#1f3a68"
ORANGE = "#e98a3a"
GREY = "#7a7a7a"
LIGHT_GREY = "#b9b9b9"


def panel_mug_to_donut(ax):
    """(a) AI biology schematic: drug pressure → cells cycle between states.

    Replaces the previous topology pedagogy (mug ≅ torus ≅ donut). The
    new image shows a cluster of identical untreated cancer cells on the
    left, an EGFR-TKI drug arrow, and a loop of differently-coloured
    cells on the right showing reversible state cycling.
    """
    _imshow_alpha_panel(
        ax, PANEL_A_PNG,
        "Drug pressure makes cancer cells cycle between states",
    )


def panel_loop_cells(ax, rng):
    """(b) AI 3D render: cells on a loop in scRNA-seq + overlay labels."""
    del rng
    _imshow_alpha_panel(
        ax, PANEL_B_PNG,
        "scRNA-seq sees the cycle as a closed loop",
        overlay_fn=_overlay_panel_b,
    )


def panel_barcode(ax, rng):
    """(c) Bar chart: max H_1 rises with drug across 4 cohorts.

    Paired (untreated, drug-treated) bars for PC9, PDX YU-006/YU-003,
    and the Maynard 14-patient cohort. Fold-change labels above each
    pair make the headline finding immediately readable.
    """
    del rng
    cohorts = [c[0] for c in COHORT_DATA]
    pre = np.array([c[1] for c in COHORT_DATA])
    post = np.array([c[2] for c in COHORT_DATA])
    folds = [c[3] for c in COHORT_DATA]

    n = len(cohorts)
    x = np.arange(n)
    w = 0.36

    ax.bar(x - w / 2, pre, w,
           label="untreated", color=NAVY,
           edgecolor="white", linewidth=0.5, zorder=3)
    ax.bar(x + w / 2, post, w,
           label="drug-treated", color=ORANGE,
           edgecolor="white", linewidth=0.5, zorder=3)

    # Fold-change labels above each cohort pair
    ymax = float(post.max())
    for i, (p, q, fold) in enumerate(zip(pre, post, folds)):
        top = max(p, q)
        ax.text(i, top + 0.30, fold,
                ha="center", va="bottom",
                fontsize=9, color=ORANGE, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, fontsize=7.5)
    ax.set_ylabel("max $H_1$ persistence", fontsize=9, color=NAVY)
    ax.set_ylim(0, ymax * 1.20)
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(frameon=False, loc="upper left", fontsize=8, handlelength=1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GREY)
    ax.spines["bottom"].set_color(GREY)
    ax.yaxis.set_tick_params(color=GREY)
    ax.grid(axis="y", color=GREY, linestyle=":", linewidth=0.4, alpha=0.5, zorder=0)

    # Subtitle (in italic) below the axis
    ax.text(
        0.5, -0.28,
        "Drug treatment raises max $H_1$ across cell line, PDX, and patient cohorts",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=8, color=NAVY, style="italic",
    )


def _draw_story_arrow(fig, x0, y0, x1, y1, text, text_offset=(0, 0)):
    """Draw a narrative connector arrow at figure level + text label."""
    arrow = mpatches.FancyArrowPatch(
        (x0, y0), (x1, y1),
        transform=fig.transFigure, figure=fig,
        arrowstyle="->", mutation_scale=12, lw=0.9, color=GREY,
        zorder=20,
    )
    fig.patches.append(arrow)
    tx, ty = (x0 + x1) / 2 + text_offset[0], (y0 + y1) / 2 + text_offset[1]
    fig.text(tx, ty, text, fontsize=7, color=GREY, style="italic",
             ha="center", va="center", zorder=21)


def main():
    set_publication_style()
    rng = np.random.default_rng(SEED)

    # Layout: 2 rows + a header strip at the top of the figure.
    # Top — A (16:9 aspect, ~1.78:1) on left, B (1:1 aspect) on right.
    # Bottom — C (filtration process) full width.
    # Figure: 180mm × 145mm (extra 10mm reserved for the story banner).
    fig = plt.figure(figsize=(180 / 25.4, 145 / 25.4))

    # Story banner: bold headline + italic three-step subtitle.
    # Communicates the paper's finding at-a-glance.
    fig.text(
        0.5, 0.975,
        "Drug-tolerant cancer cells cycle between states, quantifiable "
        "as a rising max $H_1$",
        fontsize=10.5, fontweight="bold", ha="center", va="top", color=NAVY,
    )
    fig.text(
        0.5, 0.945,
        "(a) drug pressure induces cycling  $\\rightarrow$  "
        "(b) scRNA-seq sees it as a loop  $\\rightarrow$  "
        "(c) max $H_1$ rises across 4 cohorts",
        fontsize=7.5, style="italic", ha="center", va="top", color=GREY,
    )

    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        top=0.90,  # leave room above for the story banner
        height_ratios=[1.0, 1.5],
        width_ratios=[1.78, 1.0],
        hspace=0.32, wspace=0.10,
    )

    # Step-tagged panel labels announce each panel's role
    ax_a = fig.add_subplot(gs[0, 0])
    panel_mug_to_donut(ax_a)
    ax_a.text(-0.02, 1.02, r"a $\cdot$ biology",
              transform=ax_a.transAxes,
              fontsize=9.5, fontweight="bold", va="top", color=NAVY)

    ax_b = fig.add_subplot(gs[0, 1])
    panel_loop_cells(ax_b, rng)
    ax_b.text(-0.02, 1.02, r"b $\cdot$ scRNA-seq",
              transform=ax_b.transAxes,
              fontsize=9.5, fontweight="bold", va="top", color=NAVY)

    ax_c = fig.add_subplot(gs[1, :])
    panel_barcode(ax_c, rng)
    ax_c.text(-0.005, 1.02, r"c $\cdot$ result across cohorts",
              transform=ax_c.transAxes,
              fontsize=9.5, fontweight="bold", va="top", color=NAVY)

    # Narrative connector: vertical arrow from the top row down to panel (c).
    _draw_story_arrow(
        fig, x0=0.5, y0=0.605, x1=0.5, y1=0.555,
        text="quantify with persistent homology",
        text_offset=(0.19, 0),
    )

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
