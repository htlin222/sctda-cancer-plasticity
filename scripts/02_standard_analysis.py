#!/usr/bin/env python
"""
Step 3: Standard scRNA-seq baseline analysis.

UMAP, Leiden clustering, marker genes, EMT scoring, diffusion map.
This reproduces published results and serves as the baseline comparison for TDA.

Run: python scripts/02_standard_analysis.py --seed 42
"""

import argparse
import logging
import sys

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import DATA_PROCESSED, SEED, TABLES_DIR

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Hallmark EMT gene set (MSigDB core genes)
EMT_GENES = [
    "VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1",
    "MMP2", "MMP9", "CDH1", "OCLN", "TJP1", "DSP", "KRT19", "KRT8",
    "SERPINE1", "TGFBI", "SPARC", "COL1A1", "COL3A1", "COL5A2",
    "ACTA2", "TAGLN", "MYL9", "CALD1", "LGALS1", "TPM1", "TPM2",
    "FSTL1", "ITGB1", "LAMC2", "LOXL2",
]

# Epithelial markers (for E-score)
EPITHELIAL_GENES = ["CDH1", "OCLN", "TJP1", "DSP", "KRT19", "KRT8", "EPCAM", "CLDN4"]

# Mesenchymal markers (for M-score)
MESENCHYMAL_GENES = ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1"]


def main(seed: int = SEED):
    logger.info(f"=== Standard Analysis Pipeline (seed={seed}) ===")

    # ── Load preprocessed data ──────────────────────────────
    h5ad_path = DATA_PROCESSED / "pc9_erlotinib_processed.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Run preprocessing first: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── Neighbors graph ─────────────────────────────────────
    logger.info("Computing neighbors graph...")
    sc.pp.neighbors(adata, n_pcs=30, random_state=seed)

    # ── UMAP ────────────────────────────────────────────────
    logger.info("Computing UMAP...")
    sc.tl.umap(adata, random_state=seed)

    # ── Leiden clustering at multiple resolutions ───────────
    resolutions = [0.3, 0.5, 0.8, 1.0, 1.5]
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=seed)
        n_clusters = adata.obs[key].nunique()
        logger.info(f"Leiden res={res}: {n_clusters} clusters")

    # Default clustering
    adata.obs["leiden"] = adata.obs["leiden_0.8"]

    # ── Marker genes ────────────────────────────────────────
    logger.info("Finding marker genes per cluster...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)

    # Save top markers
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(TABLES_DIR / "marker_genes_leiden.csv", index=False)
    logger.info(f"Saved marker genes to {TABLES_DIR / 'marker_genes_leiden.csv'}")

    # ── EMT scoring ─────────────────────────────────────────
    logger.info("Computing EMT scores...")

    # Full EMT score
    emt_found = [g for g in EMT_GENES if g in adata.raw.var_names]
    if emt_found:
        sc.tl.score_genes(adata, gene_list=emt_found, score_name="emt_score", use_raw=True)
        logger.info(f"EMT score: {len(emt_found)}/{len(EMT_GENES)} genes found")

    # Epithelial score
    epi_found = [g for g in EPITHELIAL_GENES if g in adata.raw.var_names]
    if epi_found:
        sc.tl.score_genes(adata, gene_list=epi_found, score_name="epithelial_score", use_raw=True)

    # Mesenchymal score
    mes_found = [g for g in MESENCHYMAL_GENES if g in adata.raw.var_names]
    if mes_found:
        sc.tl.score_genes(adata, gene_list=mes_found, score_name="mesenchymal_score", use_raw=True)

    # ── Diffusion map ───────────────────────────────────────
    logger.info("Computing diffusion map...")
    sc.tl.diffmap(adata, random_state=seed)

    # Diffusion pseudotime (root = earliest timepoint)
    if "timepoint" in adata.obs.columns:
        # Pick a root cell from the untreated/D0 population
        d0_mask = adata.obs["timepoint"].isin(["D0", "untreated", "0"])
        if d0_mask.any():
            d0_cells = adata[d0_mask]
            # Pick the cell most central in diffmap space
            root_idx = np.argmin(np.sum(d0_cells.obsm["X_diffmap"] ** 2, axis=1))
            adata.uns["iroot"] = np.where(d0_mask)[0][root_idx]
            sc.tl.dpt(adata)
            logger.info("Diffusion pseudotime computed")

    # ── PAGA ────────────────────────────────────────────────
    logger.info("Computing PAGA...")
    sc.tl.paga(adata, groups="leiden")

    # ── Save ────────────────────────────────────────────────
    out_path = DATA_PROCESSED / "pc9_erlotinib_analyzed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved analyzed data: {out_path}")

    # Summary table
    summary = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_clusters_leiden_0.8": adata.obs["leiden"].nunique(),
    }
    if "timepoint" in adata.obs.columns:
        for tp in adata.obs["timepoint"].unique():
            summary[f"cells_{tp}"] = (adata.obs["timepoint"] == tp).sum()
    if "emt_score" in adata.obs.columns:
        summary["emt_score_mean"] = float(adata.obs["emt_score"].mean())
        summary["emt_score_std"] = float(adata.obs["emt_score"].std())

    pd.Series(summary).to_csv(TABLES_DIR / "analysis_summary.csv")
    logger.info("\n=== Analysis Summary ===")
    for k, v in summary.items():
        logger.info(f"  {k}: {v}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standard scRNA-seq baseline analysis")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
