"""
24 — Publication figures for the cautionary/methods paper.
Generates the control-framework figures from locked results (scripts 16-23).
Outputs 300-dpi vector PDFs to figures/cautionary/.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "font.size": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
})
OUT = Path("figures/cautionary"); OUT.mkdir(parents=True, exist_ok=True)
RES = Path("results/foundation")

# ---- Locked numbers (scripts/17 dispersion null; scripts/21 cc-regression) ----
COH = ["osi D0","osi D3","osi D7","osi D14","erl D0","erl D9","erl D11",
       "PDX unt","PDX res","May TN","May PER","May PD"]
REAL = [1.41,1.82,2.67,3.78,1.49,1.77,1.08,2.79,5.46,5.49,6.04,5.58]
GAUSS= [1.49,2.00,2.04,2.66,1.60,1.80,1.55,2.30,2.80,3.62,3.57,3.68]
PCT  = [31,18,100,98,33,42,0,86,99,99,98,95]  # % bootstrap real>gauss

# circular coordinate R (scripts/16) across cohorts; control loop R=0.881
CC = {"erl D9":0.992,"osi D14":0.997,"PDX res":0.965,"Kim naive":0.987,"May PD":1.000}
CONTROL_R = 0.881

# pooled cell-cycle regression (scripts/21): ratio real/gauss, std vs cc-regressed
CC_STD = {"D0":1.01,"D14":1.18}; CC_REG = {"D0":1.07,"D14":0.99}

# clonal memory decay robustness grid (scripts/22): decay D0->D14 per setting
GRID = {"res0.5 k2":-0.505,"res0.5 k3":0.671,"res1.0 k2":-0.045,
        "res1.0 k3":0.218,"res1.5 k2":-0.166,"res1.5 k3":0.085}


def fig1():
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    x = np.arange(len(COH)); w = 0.38
    ax.bar(x - w/2, REAL, w, label="observed max-$H_1$", color="#c0392b")
    ax.bar(x + w/2, GAUSS, w, label="covariance-matched Gaussian null", color="#7f8c8d")
    for i, p in enumerate(PCT):
        ax.text(i, max(REAL[i], GAUSS[i]) + 0.1, f"{p}%", ha="center", va="bottom", fontsize=6,
                color="#c0392b" if p >= 95 else "#7f8c8d")
    ax.set_xticks(x); ax.set_xticklabels(COH, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("max $H_1$ persistence")
    ax.set_title("Apparent drug-associated topology vs a dispersion-matched null\n"
                 "(% = bootstrap fraction observed > null; excess only beyond baseline)", fontsize=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.axhline(0, color="k", lw=0.5)
    fig.tight_layout(); fig.savefig(OUT / "fig1_dispersion_null.pdf"); plt.close(fig)


def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    # (a) circular-coordinate negative control
    ax = axes[0]
    labels = list(CC.keys()); vals = [CC[k] for k in labels]
    ax.bar(range(len(labels)), vals, color="#2c3e50")
    ax.axhline(CONTROL_R, color="#27ae60", ls="--", lw=1.2,
               label=f"12%-occupied control loop (R={CONTROL_R})")
    ax.axhline(1.0, color="#bbb", lw=0.6)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Rayleigh R  (1 = one angle, 0 = spread)")
    ax.set_ylim(0.8, 1.02)
    ax.set_title("Negative control: cells do NOT\ntraverse the loop", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    # (b) cell-cycle confound
    ax = axes[1]
    tps = ["D0", "D14"]; x = np.arange(2); w = 0.38
    ax.bar(x - w/2, [CC_STD[t] for t in tps], w, label="standard", color="#c0392b")
    ax.bar(x + w/2, [CC_REG[t] for t in tps], w, label="S/G2M regressed", color="#2980b9")
    ax.axhline(1.0, color="#888", ls=":", lw=1, label="= dispersion null")
    ax.set_xticks(x); ax.set_xticklabels(tps)
    ax.set_ylabel("max $H_1$ / Gaussian-null ratio")
    ax.set_title("Drug excess is the cell-cycle loop:\nit vanishes after S/G2M regression", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    ax.set_ylim(0.8, 1.3)
    fig.tight_layout(); fig.savefig(OUT / "fig2_controls.pdf"); plt.close(fig)


def fig3():
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    labels = list(GRID.keys()); vals = [GRID[k] for k in labels]
    colors = ["#27ae60" if v > 0 else "#c0392b" for v in vals]
    ax.barh(range(len(labels)), vals, color=colors)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("clonal memory decay  D0$\\to$D14")
    ax.set_title("Clonal memory-decay is not robust:\nsign flips across clustering settings (3/6 > 0)", fontsize=8)
    ax.legend(handles=[Patch(color="#27ae60", label="erodes (>0)"),
                       Patch(color="#c0392b", label="no erosion ($\\leq$0)")],
              frameon=False, fontsize=6, loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / "fig3_memory_nonrobust.pdf"); plt.close(fig)


def fig4():
    """Synthetic ground-truth benchmark (reads results/foundation/synthetic_benchmark.json)."""
    j = json.loads((RES / "synthetic_benchmark.json").read_text())
    order = ["a_gaussian_blob", "b_traversed_loop", "c_sparse_arc", "d_two_clusters"]
    names = ["Gaussian blob\n(no structure)", "traversed loop\n(real cycle)",
             "sparse arc\n(structure, not cyclic)", "two clusters\n(non-Gaussian)"]
    gr = [j[k]["gaussian_ratio"] for k in order]
    R = [j[k]["rayleigh_R"] for k in order]
    verdict = [j[k]["toolkit_verdict"] for k in order]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    x = np.arange(4)
    ax = axes[0]
    ax.bar(x, gr, color=["#7f8c8d", "#c0392b", "#e67e22", "#7f8c8d"])
    ax.axhline(1.0, color="k", ls=":", lw=1, label="= dispersion null")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("max $H_1$ / Gaussian-null ratio")
    ax.set_title("Structure detection (ratio > 1)", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    ax = axes[1]
    ax.bar(x, R, color=["#7f8c8d", "#c0392b", "#e67e22", "#7f8c8d"])
    ax.axhline(0.6, color="#27ae60", ls="--", lw=1, label="traversal threshold")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("Rayleigh R")
    ax.set_title("Traversal test (R < 0.6 = traversed)", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    ax.set_ylim(0, 1.05)
    fig.suptitle("Synthetic ground truth: toolkit is sensitive (loop) and specific (blob, clusters)", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig4_synthetic_benchmark.pdf"); plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("wrote:", *[p.name for p in sorted(OUT.glob("*.pdf"))])
