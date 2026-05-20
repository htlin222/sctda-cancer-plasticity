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
import matplotlib.pyplot as plt
import numpy as np

from sctda_plasticity.config import FIGURES_DIR, PROJECT_ROOT, SEED
from sctda_plasticity.visualize import save_figure, set_publication_style

PANEL_A_PNG = FIGURES_DIR / "fig1_panel_a_3d.png"
PANEL_B_PNG = FIGURES_DIR / "fig1_panel_b_3d.png"
PANEL_C_PNG = FIGURES_DIR / "fig1_panel_c_3d.png"


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


def _overlay_panel_c(ax, w, h):
    """Add explanatory callouts on top of the persistence-barcode image."""
    # Long orange bar — call it out as max H_1
    ax.annotate(
        "max $H_1$ persistence\n(this paper's statistic)",
        xy=(w * 0.50, h * 0.18), xytext=(w * 1.04, h * 0.08),
        fontsize=6.5, color=ORANGE, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", lw=0.6, color=ORANGE,
                        shrinkA=0, shrinkB=2),
    )
    # Short navy stack — call out as noise
    ax.annotate(
        "noise\n($H_0$ components)",
        xy=(w * 0.55, h * 0.65), xytext=(w * 1.04, h * 0.70),
        fontsize=6.5, color=NAVY, ha="left", va="center",
        arrowprops=dict(arrowstyle="->", lw=0.6, color=NAVY,
                        shrinkA=0, shrinkB=2),
    )
    # Filtration-scale axis hint below the bars
    ax.annotate(
        "", xy=(w * 0.95, h * 0.96), xytext=(w * 0.05, h * 0.96),
        arrowprops=dict(arrowstyle="->", lw=0.5, color=GREY),
    )
    ax.text(
        w * 0.5, h * 1.01, r"filtration scale $\varepsilon$",
        fontsize=6.5, color=GREY, ha="center", va="top",
    )

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
    """(c) AI 3D barcode + explanatory overlay labels (max H_1, noise, ε)."""
    del rng
    _imshow_alpha_panel(
        ax, PANEL_C_PNG,
        "Loop $\\Rightarrow$ long $H_1$ bar",
        overlay_fn=_overlay_panel_c,
    )


def main():
    set_publication_style()
    rng = np.random.default_rng(SEED)

    # Layout: panel (a) is 3:1 (wide), (b) is 1:1 (square), (c) is 3:2.
    # Aspect-matched widths so each panel reads at similar visual height.
    # Total target: 180mm × 75mm.
    fig = plt.figure(figsize=(180 / 25.4, 75 / 25.4))
    gs = gridspec.GridSpec(
        1, 3, figure=fig,
        width_ratios=[3.0, 1.8, 2.5],
        wspace=0.05,
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
