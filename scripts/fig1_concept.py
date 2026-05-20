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

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Editorial palette
NAVY = "#1f3a68"
ORANGE = "#e98a3a"
GREY = "#7a7a7a"
LIGHT_GREY = "#b9b9b9"


def panel_mug_to_donut(ax):
    """(a) Inset the AI 3D render: mug ≅ donut-with-handle ≅ donut.

    The image was generated via codex-imagegen (OpenAI gpt-image-1) on a
    chroma-key background and matted to alpha. See declarations.md for
    AI-tool disclosure. All three shapes carry exactly one topological
    hole; the congruence symbols (≅) between them indicate homeomorphism.
    """
    if not PANEL_A_PNG.exists():
        ax.text(
            0.5, 0.5,
            f"Missing AI panel-a PNG:\n{PANEL_A_PNG}\n"
            "(generate via scripts/imagegen.sh)",
            ha="center", va="center", color="red",
            transform=ax.transAxes,
        )
    else:
        img = mpimg.imread(str(PANEL_A_PNG))
        ax.imshow(img, interpolation="bilinear", aspect="equal")
        # subtitle below the rendered shapes
        h = img.shape[0]
        w = img.shape[1]
        ax.text(
            w / 2, h * 1.06, "Same hole, different geometry",
            fontsize=7, color=NAVY, style="italic",
            ha="center", va="top",
        )
        ax.set_xlim(0, w)
        # extra y headroom for the subtitle
        ax.set_ylim(h * 1.18, -h * 0.04)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_loop_cells(ax, rng):
    """(b) 60 cells arranged on a noisy unit circle in state-space."""
    n = 60
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r = 1.0 + rng.normal(0, 0.07, n)
    t = theta + rng.normal(0, 0.07, n)
    x = r * np.cos(t)
    y = r * np.sin(t)
    cmap = plt.get_cmap("twilight_shifted")
    colors = cmap(theta / (2 * np.pi))
    ax.scatter(x, y, c=colors, s=22, edgecolors="white", linewidths=0.4, zorder=3)
    # axes with arrowheads (visual hint that these are coordinate axes)
    ax.annotate("", xy=(1.55, 0), xytext=(-1.55, 0),
                arrowprops=dict(arrowstyle="->", lw=0.5, color=LIGHT_GREY),
                zorder=1)
    ax.annotate("", xy=(0, 1.55), xytext=(0, -1.55),
                arrowprops=dict(arrowstyle="->", lw=0.5, color=LIGHT_GREY),
                zorder=1)
    # axis labels (data-coord)
    ax.text(1.55, -0.18, "cell-state\ncoord. 1", fontsize=6, color=GREY,
            ha="right", va="top")
    ax.text(-0.05, 1.55, "cell-state\ncoord. 2", fontsize=6, color=GREY,
            ha="right", va="top", rotation=0)
    # subtitle
    ax.text(0, -1.95, "A topological loop in data",
            fontsize=7, color=NAVY, style="italic",
            ha="center", va="top")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-2.1, 1.9)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_barcode(ax, rng):
    """(c) Persistence barcode: short H_0 bars + one long H_1 bar."""
    n_short = 8
    starts = rng.uniform(0.04, 0.18, n_short)
    lengths = rng.uniform(0.10, 0.30, n_short)
    order = np.argsort(lengths)[::-1]
    bar_dy = 0.075
    for i, idx in enumerate(order):
        y = 0.18 + i * bar_dy
        ax.plot(
            [starts[idx], starts[idx] + lengths[idx]], [y, y],
            color=NAVY, lw=2.4, solid_capstyle="round",
        )
    y_long = 0.18 + n_short * bar_dy + 0.06
    ax.plot([0.06, 0.96], [y_long, y_long],
            color=ORANGE, lw=3.6, solid_capstyle="round")
    # inline labels pointing at the two bar populations
    ax.annotate(
        "robust loop\n($H_1$)",
        xy=(0.55, y_long), xytext=(1.05, y_long - 0.02),
        fontsize=6.5, color=ORANGE, ha="left", va="center",
        arrowprops=dict(arrowstyle="-", lw=0.5, color=ORANGE),
    )
    ax.annotate(
        "noise / components\n($H_0$)",
        xy=(0.30, 0.18 + 3 * bar_dy), xytext=(1.05, 0.18 + 3 * bar_dy),
        fontsize=6.5, color=NAVY, ha="left", va="center",
        arrowprops=dict(arrowstyle="-", lw=0.5, color=NAVY),
    )
    # x-axis arrow + label (filtration scale)
    ax.annotate("", xy=(1.0, 0.08), xytext=(0.0, 0.08),
                arrowprops=dict(arrowstyle="->", lw=0.5, color=GREY))
    ax.text(0.5, 0.03, r"filtration scale $\varepsilon$",
            fontsize=7, color=GREY, ha="center", va="top")
    # subtitle below
    ax.text(0.5, -0.07, "Persistence barcode",
            fontsize=7, color=NAVY, style="italic",
            ha="center", va="top")
    ax.set_xlim(0, 2.0)  # extra room on right for the inline labels
    ax.set_ylim(-0.12, y_long + 0.12)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    set_publication_style()
    rng = np.random.default_rng(SEED)

    # Taller layout (180 × 75 mm) gives room for axis labels and subtitles
    fig = plt.figure(figsize=(180 / 25.4, 75 / 25.4))
    gs = gridspec.GridSpec(
        1, 3, figure=fig,
        width_ratios=[1.7, 1.0, 1.4],
        wspace=0.12,
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
