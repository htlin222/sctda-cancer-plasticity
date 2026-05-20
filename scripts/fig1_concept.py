#!/usr/bin/env python
"""
Fig 1 concept: Topology distinguishes shapes that geometry cannot.

Three panels:
  (a) Photorealistic 3D render of a coffee mug, a donut-with-handle
      intermediate, and a donut connected by congruence symbols ---
      illustrating that homeomorphism preserves the single hole even as
      geometry changes. The image is AI-generated (OpenAI gpt-image-1)
      with chroma-key background removal; provenance and disclosure are
      in manuscript/submission/declarations.md.
  (b) Single-cell transcriptomes lying on a noisy closed loop in
      cell-state coordinates (matplotlib scatter).
  (c) A persistence barcode: short H_0 bars (noise) and one long
      H_1 bar (the robust loop) (matplotlib).

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

PANEL_A_PNG = FIGURES_DIR / "fig1_panel_a_3d.png"
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
    """(a) AI 3D render: mug ≅ pinched torus ≅ donut. All genus = 1."""
    _imshow_alpha_panel(ax, PANEL_A_PNG, "Same hole, different geometry")


def panel_loop_cells(ax, rng):
    """(b) AI 3D render of cells on a loop + explanatory overlay labels."""
    del rng
    _imshow_alpha_panel(
        ax, PANEL_B_PNG,
        "Cells trace a state-space loop",
        overlay_fn=_overlay_panel_b,
    )


def panel_barcode(ax, rng):
    """(c) Filtration process → max H_1 persistence bar.

    Three snapshots show the same cell loop with growing filtration
    balls. Stage 1: small balls, cells are isolated components (H_0).
    Stage 2: balls touch their neighbours' balls and the loop closes
    --- the H_1 feature is *born*. Stage 3: balls fill the loop's
    interior --- the H_1 feature *dies*. The orange bar above spans
    birth-to-death = max H_1 persistence.
    """
    del rng
    # Wide layout: data range 0..3.0 in x to match the wide bottom-row
    # cell. Snapshots spaced across the full width.
    s_y = 0.45
    s_x = [0.55, 1.50, 2.45]
    loop_r = 0.18
    ball_radii = [0.045, 0.115, 0.245]  # small / just-touching / interior-filling

    # Three filtration snapshots
    for cx, br in zip(s_x, ball_radii):
        _draw_filtration_step(ax, cx, s_y, loop_r, br, n_cells=12)
    # Stage labels under each snapshot — plain language first, technical
    # ε-notation in small parens for the specialist reader only.
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

    # Plain-language axis at the bottom (technical name kept in tiny grey)
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
        r"(``filtration scale $\varepsilon$'' in TDA)",
        fontsize=6, color=GREY, style="italic", ha="center", va="top",
    )

    # The orange bar: spans from snapshot 2 (loop appears) to snapshot 3
    # (loop disappears) — its length is the loop's lifespan.
    bar_y = s_y + loop_r + ball_radii[-1] + 0.10
    bar_x0, bar_x1 = s_x[1], s_x[2]
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color=ORANGE, lw=8, solid_capstyle="round", zorder=8)
    # Vertical drop-lines tying the bar to the axis below
    for x in (bar_x0, bar_x1):
        ax.plot([x, x], [bar_y - 0.04, eps_axis_y],
                color=GREY, lw=0.6, linestyle=":", alpha=0.7, zorder=1)
    # Endpoint markers labelled in plain language
    endpoint_labels = ["loop appears", "loop disappears"]
    for x, txt in zip([bar_x0, bar_x1], endpoint_labels):
        ax.plot([x, x], [eps_axis_y - 0.04, eps_axis_y + 0.04],
                color=GREY, lw=0.8, zorder=2)
        ax.text(x, eps_axis_y + 0.06, txt,
                fontsize=6.5, color=GREY, ha="center", va="bottom")
    # Plain-language bar label, with technical name parenthesised
    ax.text(
        (bar_x0 + bar_x1) / 2, bar_y + 0.22,
        "How long the loop survives",
        fontsize=9, color=ORANGE, fontweight="bold",
        ha="center", va="bottom",
    )
    ax.text(
        (bar_x0 + bar_x1) / 2, bar_y + 0.13,
        r"$=$ this paper's ``max $H_1$ persistence''",
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


def main():
    set_publication_style()
    rng = np.random.default_rng(SEED)

    # Layout: 2 rows. Top — A (16:9 aspect, ~1.78:1) on left, B (1:1
    # aspect) on right. Bottom — C (filtration process) full width.
    # Width ratios match the image aspects so both A and B fill their
    # cells with minimal whitespace. Figure: 180mm × 135mm.
    fig = plt.figure(figsize=(180 / 25.4, 135 / 25.4))
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        height_ratios=[1.0, 1.5],
        width_ratios=[1.78, 1.0],
        hspace=0.18, wspace=0.05,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    panel_mug_to_donut(ax_a)
    ax_a.text(-0.02, 1.02, "a", transform=ax_a.transAxes,
              fontsize=10, fontweight="bold", va="top")

    ax_b = fig.add_subplot(gs[0, 1])
    panel_loop_cells(ax_b, rng)
    ax_b.text(-0.02, 1.02, "b", transform=ax_b.transAxes,
              fontsize=10, fontweight="bold", va="top")

    ax_c = fig.add_subplot(gs[1, :])
    panel_barcode(ax_c, rng)
    ax_c.text(-0.005, 1.02, "c", transform=ax_c.transAxes,
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
