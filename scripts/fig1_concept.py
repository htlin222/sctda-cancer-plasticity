#!/usr/bin/env python
"""
Fig 1 — Concept figure: drug-induced cell-state cycling, detected as a
topological loop, measured by persistent homology.

Three panels (concept only — no headline numbers; the cross-cohort
result lives in Fig 7 `fig:master`):
  (a) AI render: untreated cancer cells (uniform cluster) → drug
      arrow → cells cycling between distinct states. The biological
      motivation.
  (b) AI render: the same cycling appears as a closed loop in
      scRNA-seq state-space.
  (c) Matplotlib: persistent-homology filtration process. Three
      snapshots show cells with growing balls; one orange bar above
      gives the loop's lifespan = max H_1 persistence.

The mug-donut topological-equivalence analogy is a single line in the
Background prose, not Figure 1.

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
    """(c) Filtration process: cells + growing balls → max H_1 bar.

    Three matplotlib snapshots show the same cell loop with growing
    filtration balls. Stage 1: small balls, cells isolated. Stage 2:
    balls touch → loop appears. Stage 3: balls fill the interior →
    loop disappears. Orange bar above spans birth→death = max H_1
    persistence. Plain-language labels throughout; the technical
    ε / H_1 names sit in tiny grey parentheticals.
    """
    del rng
    s_y = 0.45
    s_x = [0.55, 1.50, 2.45]
    loop_r = 0.18
    ball_radii = [0.045, 0.115, 0.245]

    for cx, br in zip(s_x, ball_radii):
        _draw_filtration_step(ax, cx, s_y, loop_r, br, n_cells=12)
    stage_labels = [
        "small circles:\neach cell alone",
        "circles touch:\nloop appears",
        "circles fill in:\nloop disappears",
    ]
    for cx, label in zip(s_x, stage_labels):
        ax.text(
            cx, s_y - loop_r - ball_radii[-1] - 0.05, label,
            fontsize=7.5, color=NAVY, ha="center", va="top",
        )

    eps_axis_y = -0.42
    ax.annotate(
        "", xy=(2.90, eps_axis_y), xytext=(0.10, eps_axis_y),
        arrowprops=dict(arrowstyle="->", lw=0.7, color=GREY),
    )
    ax.text(
        1.50, eps_axis_y - 0.05,
        "circle size grows $\\rightarrow$",
        fontsize=8, color=GREY, ha="center", va="top",
    )
    ax.text(
        1.50, eps_axis_y - 0.16,
        r"(filtration scale $\varepsilon$)",
        fontsize=6, color=GREY, style="italic", ha="center", va="top",
    )

    bar_y = s_y + loop_r + ball_radii[-1] + 0.10
    bar_x0, bar_x1 = s_x[1], s_x[2]
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color=ORANGE, lw=8, solid_capstyle="round", zorder=8)
    for x in (bar_x0, bar_x1):
        ax.plot([x, x], [bar_y - 0.04, eps_axis_y],
                color=GREY, lw=0.6, linestyle=":", alpha=0.7, zorder=1)
    endpoint_labels = ["loop appears", "loop disappears"]
    for x, txt in zip([bar_x0, bar_x1], endpoint_labels):
        ax.plot([x, x], [eps_axis_y - 0.04, eps_axis_y + 0.04],
                color=GREY, lw=0.8, zorder=2)
        ax.text(x, eps_axis_y + 0.06, txt,
                fontsize=6.5, color=GREY, ha="center", va="bottom")
    ax.text(
        (bar_x0 + bar_x1) / 2, bar_y + 0.22,
        "How long the loop survives",
        fontsize=9, color=ORANGE, fontweight="bold",
        ha="center", va="bottom",
    )
    ax.text(
        (bar_x0 + bar_x1) / 2, bar_y + 0.13,
        r"(max $H_1$ persistence)",
        fontsize=7, color=ORANGE, style="italic",
        ha="center", va="bottom",
    )

    ax.set_xlim(0, 3.0)
    ax.set_ylim(eps_axis_y - 0.25, bar_y + 0.40)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


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
    # Frames the figure as a concept introduction (the cross-cohort
    # result itself is in Fig 7, fig:master).
    fig.text(
        0.5, 0.975,
        "Drug-induced cell-state cycling, detected as a topological "
        "loop and measured by persistent homology",
        fontsize=10.5, fontweight="bold", ha="center", va="top", color=NAVY,
    )
    fig.text(
        0.5, 0.945,
        "(a) drug pressure induces cycling  $\\rightarrow$  "
        "(b) scRNA-seq sees it as a loop  $\\rightarrow$  "
        "(c) persistent homology measures the loop's lifespan",
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
    ax_c.text(-0.005, 1.02, r"c $\cdot$ how we measure the loop",
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
