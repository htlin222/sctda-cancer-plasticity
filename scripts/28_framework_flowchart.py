"""
28 — Framework decision flowchart (comprehension aid, Fig 5).
A top-down decision tree: input -> dispersion test -> traversal test + cell-cycle
test -> verdict, all wrapped in bootstrap. Pure matplotlib, no data needed.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.size": 8,
                     "pdf.fonttype": 42, "font.family": "sans-serif"})
OUT = Path("figures/cautionary"); OUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(7.2, 5.4)); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, ec="#333", fs=8, bold=False):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec=ec, lw=1.1))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", wrap=True)

def arrow(x1, y1, x2, y2, label="", lx=0, ly=0, color="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                 lw=1.1, color=color))
    if label:
        ax.text((x1+x2)/2 + lx, (y1+y2)/2 + ly, label, ha="center", va="center",
                fontsize=7, style="italic", color=color)

STEP = "#dbe8f2"; STOP = "#e0e0e0"; REAL = "#cdebd2"; HEAD = "#f5d9c4"

box(5, 9.4, 6.2, 0.8, "Input: cell × PC embedding, per condition\n(wrapped in bootstrap subsampling)", HEAD, fs=8, bold=True)
box(5, 8.0, 4.6, 0.9, "① Dispersion test\nobserved max-$H_1$ / Gaussian-null ratio", STEP, bold=True)
arrow(5, 9.0, 5, 8.5)

box(1.7, 6.6, 2.8, 0.85, "No structure\nbeyond dispersion", STOP)
arrow(5, 7.55, 2.2, 7.0, "ratio ≈ 1", lx=0.1, ly=0.25)

box(5, 6.4, 4.6, 0.9, "② Traversal test\nRayleigh R of circular coordinates", STEP, bold=True)
arrow(5, 7.55, 5, 6.9, "ratio > 1\n(robustly)", lx=1.05, ly=0)

box(2.4, 4.9, 3.0, 0.85, "Structure,\nbut NOT a traversed loop", STEP)
arrow(5, 5.95, 2.7, 5.35, "R ≥ 0.6", lx=0.15, ly=0.25)

box(6.6, 5.0, 4.0, 0.85, "Traversed loop\n(the only 'cyclic' verdict)", REAL)
arrow(5, 5.95, 6.5, 5.45, "R < 0.6", lx=0.55, ly=0.15)

box(2.4, 3.2, 4.4, 0.95, "③ Cell-cycle test\nexcess survives S/G2M regression?", STEP, bold=True)
arrow(2.4, 4.45, 2.4, 3.7)

box(1.5, 1.6, 2.6, 0.85, "Excess is\nthe cell cycle", STOP)
arrow(2.4, 2.7, 1.6, 2.05, "no", lx=-0.25, ly=0.15)
box(5.6, 1.6, 3.4, 0.85, "Genuine non-cell-cycle\nstructure", REAL)
arrow(2.4, 2.7, 5.4, 2.05, "yes", lx=0.3, ly=0.2)

ax.text(5, 0.5, "Decision rule: 'cyclic plasticity' requires a traversed loop (green, top-right). "
        "Everything else is dispersion, cell cycle, or non-traversed structure.",
        ha="center", va="center", fontsize=7.5, style="italic", color="#444")

fig.tight_layout()
fig.savefig(OUT / "fig5_framework_flowchart.pdf", bbox_inches="tight")
print("wrote", OUT / "fig5_framework_flowchart.pdf")
