#!/usr/bin/env python
"""
Harmony batch correction sensitivity analysis.

Re-computes max H1 on Harmony-corrected embeddings for:
1. PC9 osimertinib (batch = timepoint)  → test if monotonic trend survives
2. PDX ASCL1 (batch = pdx + condition as batch variables)
3. Kim 2020 LUAD atlas (batch = Sample = patient) → test if H1=4.77 survives

Purpose: directly address R3-R3 concern that cross-platform/inter-patient
variance confounds the monotonic H1 claim.

Run: python scripts/13_harmony_batch_correction.py --seed 42
"""

import argparse
import logging
import sys
import warnings

sys.path.insert(0, "src")
warnings.filterwarnings("ignore")

import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    FIGURES_DIR,
    PH_MAX_DIM,
    PH_N_PCS,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.topology import compute_persistent_homology
from sctda_plasticity.visualize import save_figure, set_publication_style

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def harmony_correct(adata: sc.AnnData, batch_key: str, seed: int = 42) -> np.ndarray:
    """Return Harmony-corrected PCA embedding."""
    meta = adata.obs[[batch_key]].copy()
    pca = adata.obsm["X_pca"]
    harmonized = hm.run_harmony(pca, meta, batch_key, random_state=seed, verbose=False)
    z = harmonized.Z_corr
    # harmonypy sometimes returns (n_cells, n_pcs), sometimes (n_pcs, n_cells); align
    if z.shape[0] == meta.shape[0]:
        return z
    return z.T


def main(seed: int = SEED):
    logger.info(f"=== Harmony batch correction sensitivity (seed={seed}) ===")

    results = []
    rng = np.random.RandomState(seed)

    # ── 1. Osimertinib: Harmony with batch=timepoint ────────
    logger.info("\n--- 1. PC9 osimertinib (batch=timepoint) ---")
    adata = sc.read_h5ad(DATA_PROCESSED / "pc9_osimertinib_processed.h5ad")
    adata.var_names_make_unique()
    X_harmony = harmony_correct(adata, "timepoint", seed=seed)
    adata.obsm["X_harmony"] = X_harmony

    for tp in ["D0", "D3", "D7", "D14"]:
        cells = adata[adata.obs["timepoint"] == tp]
        if cells.n_obs > 1000:
            cells = cells[rng.choice(cells.n_obs, 1000, replace=False)]

        # Raw PCA H1
        ph_raw = compute_persistent_homology(
            cells.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
        )
        h1_raw = ph_raw["max_persistence"].get(1, 0)

        # Harmony-corrected H1
        ph_har = compute_persistent_homology(
            cells.obsm["X_harmony"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
        )
        h1_har = ph_har["max_persistence"].get(1, 0)

        logger.info(f"  Osi {tp}: raw H1={h1_raw:.3f}, harmony H1={h1_har:.3f}")
        results.append({
            "dataset": "PC9+osimertinib",
            "group": tp,
            "n_cells": cells.n_obs,
            "max_h1_raw_pca": h1_raw,
            "max_h1_harmony": h1_har,
        })

    # ── 2. PDX ASCL1: batch=sample (encodes pdx × condition) ─
    logger.info("\n--- 2. ASCL1 PDX (batch=sample) ---")
    adata_pdx = sc.read_h5ad(DATA_PROCESSED / "pdx_ascl1_osimertinib_processed.h5ad")
    adata_pdx.var_names_make_unique()
    X_harmony_pdx = harmony_correct(adata_pdx, "sample", seed=seed)
    adata_pdx.obsm["X_harmony"] = X_harmony_pdx

    for pdx in sorted(adata_pdx.obs["pdx"].unique()):
        for cond in sorted(adata_pdx.obs["condition"].unique()):
            mask = (adata_pdx.obs["pdx"] == pdx) & (adata_pdx.obs["condition"] == cond)
            cells = adata_pdx[mask]
            if cells.n_obs < 30:
                continue
            if cells.n_obs > 1500:
                cells = cells[rng.choice(cells.n_obs, 1500, replace=False)]

            h1_raw = compute_persistent_homology(
                cells.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
            )["max_persistence"].get(1, 0)
            h1_har = compute_persistent_homology(
                cells.obsm["X_harmony"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
            )["max_persistence"].get(1, 0)

            logger.info(f"  PDX {pdx} {cond}: raw H1={h1_raw:.3f}, harmony H1={h1_har:.3f}")
            results.append({
                "dataset": f"PDX {pdx}",
                "group": cond,
                "n_cells": cells.n_obs,
                "max_h1_raw_pca": h1_raw,
                "max_h1_harmony": h1_har,
            })

    # ── 3. Kim atlas: batch=Sample (patient) ────────────────
    logger.info("\n--- 3. Kim LUAD atlas (batch=patient Sample) ---")
    adata_kim = sc.read_h5ad(DATA_PROCESSED / "kim_luad_atlas_processed.h5ad")
    adata_kim.var_names_make_unique()
    X_harmony_kim = harmony_correct(adata_kim, "Sample", seed=seed)
    adata_kim.obsm["X_harmony"] = X_harmony_kim

    # Subsample to 1500 for PH tractability
    if adata_kim.n_obs > 1500:
        idx = rng.choice(adata_kim.n_obs, 1500, replace=False)
        sub = adata_kim[idx]
    else:
        sub = adata_kim

    h1_raw = compute_persistent_homology(
        sub.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
    )["max_persistence"].get(1, 0)
    h1_har = compute_persistent_homology(
        sub.obsm["X_harmony"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
    )["max_persistence"].get(1, 0)

    logger.info(f"  Kim LUAD naive: raw H1={h1_raw:.3f}, harmony H1={h1_har:.3f}")
    results.append({
        "dataset": "Kim LUAD",
        "group": "naive",
        "n_cells": sub.n_obs,
        "max_h1_raw_pca": h1_raw,
        "max_h1_harmony": h1_har,
    })

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(TABLES_DIR / "harmony_batch_correction_comparison.csv", index=False)

    logger.info("\n=== Summary: raw vs Harmony-corrected max H1 ===")
    logger.info(df.to_string(index=False))

    # Plot: raw vs harmony comparison
    import matplotlib.pyplot as plt

    set_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 75 / 25.4))

    # Panel A: osimertinib time course (raw vs harmony)
    ax = axes[0]
    osi_df = df[df["dataset"] == "PC9+osimertinib"].copy()
    osi_df["day"] = osi_df["group"].str[1:].astype(int)
    osi_df = osi_df.sort_values("day")
    ax.plot(
        osi_df["day"], osi_df["max_h1_raw_pca"], "o-",
        color="#2196F3", label="Raw PCA", lw=1.5, markersize=6,
    )
    ax.plot(
        osi_df["day"], osi_df["max_h1_harmony"], "s--",
        color="#E91E63", label="Harmony-corrected", lw=1.5, markersize=6,
    )
    ax.set_xlabel("Days of osimertinib")
    ax.set_ylabel(r"Max $H_1$ persistence")
    ax.set_title("a  Osimertinib monotonic trend", loc="left", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.2)

    # Panel B: PDX raw vs harmony by group
    ax = axes[1]
    pdx_df = df[df["dataset"].str.startswith("PDX")].copy()
    x = np.arange(len(pdx_df))
    w = 0.35
    ax.bar(x - w / 2, pdx_df["max_h1_raw_pca"], w, label="Raw PCA", color="#2196F3", alpha=0.8)
    ax.bar(x + w / 2, pdx_df["max_h1_harmony"], w, label="Harmony", color="#E91E63", alpha=0.8)
    labels = [f"{r['dataset'].replace('PDX ','')}\n{r['group']}" for _, r in pdx_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5)
    ax.set_ylabel(r"Max $H_1$ persistence")
    ax.set_title("b  PDX raw vs Harmony", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    fig.tight_layout()
    save_figure(fig, "figS_harmony_batch_correction")
    plt.close(fig)

    logger.info("\nFigure saved: figS_harmony_batch_correction.pdf/.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
