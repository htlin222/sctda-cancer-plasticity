#!/usr/bin/env python
"""
Validation reference: Kim 2020 LUAD atlas (GSE131907).

208,506 cells from 44 patients, treatment-naive. Used as:
1. Baseline reference — H1 should be LOWER in treatment-naive patients
2. EGFR tumor subgroup — compare to treated datasets

Run: python scripts/09_validation_kim_atlas.py --seed 42
"""

import argparse
import gzip
import logging
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
import pyreadr
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    DATA_RAW,
    MAX_GENES,
    MAX_PCT_MITO,
    MIN_CELLS,
    MIN_GENES,
    N_PCS,
    N_TOP_GENES,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.data import preprocess, qc_filter
from sctda_plasticity.topology import compute_persistent_homology
from sctda_plasticity.visualize import (
    plot_persistence_barcodes,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_kim_atlas(data_dir: Path, subsample: int | None = None, seed: int = 42) -> sc.AnnData:
    """Load Kim 2020 LUAD atlas.

    Only loads tumor cells (Sample_Origin == 'tLung' or 'mLN' or 'mBrain')
    to keep memory manageable. Subsample if requested.
    """
    logger.info("Loading cell annotation...")
    annot = pd.read_csv(
        data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t"
    )
    annot = annot.set_index("Index")
    logger.info(f"Annotation: {annot.shape[0]} cells")
    logger.info(f"Sample origins: {annot['Sample_Origin'].value_counts().to_dict()}")
    logger.info(f"Cell types: {annot['Cell_type'].value_counts().to_dict()}")

    # Focus on primary lung tumor epithelial cells (where topology matters)
    # Use Sample_Origin = 'tLung' (primary tumor) and Cell_type containing 'Epithelial'
    tumor_mask = (annot["Sample_Origin"] == "tLung") & (
        annot["Cell_type"].str.contains("Epithelial", case=False, na=False)
    )
    tumor_cells = annot[tumor_mask].index.tolist()
    logger.info(f"Tumor epithelial cells: {len(tumor_cells)}")

    # Subsample
    rng = np.random.RandomState(seed)
    if subsample and len(tumor_cells) > subsample:
        tumor_cells = rng.choice(tumor_cells, subsample, replace=False).tolist()
        logger.info(f"Subsampled to: {len(tumor_cells)} cells")

    # Load matrix
    rds = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
    rds_gz = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"
    if not rds.exists():
        logger.info("Decompressing RDS...")
        with gzip.open(rds_gz, "rb") as fi, open(rds, "wb") as fo:
            shutil.copyfileobj(fi, fo)

    logger.info("Reading RDS (may take ~60s)...")
    result = pyreadr.read_r(str(rds))
    df = list(result.values())[0]
    logger.info(f"Matrix: {df.shape[0]} genes × {df.shape[1]} cells")

    # Subset columns (cells) to the ones we want
    keep_cells = [c for c in tumor_cells if c in df.columns]
    logger.info(f"Extracting {len(keep_cells)} cells from full matrix...")
    sub_df = df[keep_cells]

    # Build AnnData
    adata = sc.AnnData(sub_df.T)  # cells × genes
    adata.var_names_make_unique()

    # Attach metadata
    meta_sub = annot.loc[keep_cells]
    for col in ["Sample", "Sample_Origin", "Cell_type", "Cell_type.refined", "Cell_subtype"]:
        if col in meta_sub.columns:
            adata.obs[col] = meta_sub[col].values

    adata.obs["timepoint"] = "untreated"  # all Kim samples are treatment-naive
    logger.info(f"Built AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
    logger.info(f"Samples: {adata.obs['Sample'].value_counts().head(10).to_dict()}")
    return adata


def main(seed: int = SEED, subsample: int = 5000):
    logger.info(f"=== Kim 2020 LUAD Atlas Validation (seed={seed}) ===")

    data_dir = DATA_RAW / "GSE131907"
    adata = load_kim_atlas(data_dir, subsample=subsample, seed=seed)

    # ── QC + Preprocess ─────────────────────────────────────
    logger.info("\n--- QC + Preprocessing ---")
    adata = qc_filter(
        adata, min_genes=MIN_GENES, max_genes=MAX_GENES,
        max_pct_mito=MAX_PCT_MITO, min_cells=MIN_CELLS,
    )
    logger.info(f"After QC: {adata.n_obs} cells")

    adata = preprocess(
        adata, n_top_genes=N_TOP_GENES, n_pcs=N_PCS,
        batch_key="Sample", seed=seed,
    )

    out_path = DATA_PROCESSED / "kim_luad_atlas_processed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved: {out_path}")

    # ── PH on treatment-naive atlas ────────────────────────
    logger.info("\n--- Persistent Homology ---")

    ph_dir = RESULTS_DIR / "ph_kim_atlas"
    ph_dir.mkdir(parents=True, exist_ok=True)

    # Subsample for PH (too many cells for full PH)
    MAX_PH = 1500
    rng = np.random.RandomState(seed)
    if adata.n_obs > MAX_PH:
        idx = rng.choice(adata.n_obs, MAX_PH, replace=False)
        adata_ph = adata[idx].copy()
    else:
        adata_ph = adata

    ph = compute_persistent_homology(
        adata_ph.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
    )
    max_h1 = ph["max_persistence"].get(1, 0)
    logger.info(f"Treatment-naive LUAD atlas: {adata_ph.n_obs} cells, max H₁ = {max_h1:.4f}")

    for dim in range(PH_MAX_DIM + 1):
        np.save(ph_dir / f"dgm{dim}_kim_naive.npy", ph["dgms"][dim])

    pd.Series({
        "n_cells": adata_ph.n_obs,
        "max_h0": ph["max_persistence"].get(0, 0),
        "max_h1": max_h1,
        "n_h1": len(ph["dgms"][1]),
    }).to_csv(TABLES_DIR / "ph_kim_atlas.csv")

    # Plot barcodes
    fig = plot_persistence_barcodes(
        {"Kim treatment-naive LUAD": ph["dgms"]},
        max_dim=PH_MAX_DIM,
        title="Kim 2020 Atlas (Treatment-Naive)",
    )
    save_figure(fig, "fig7_d_kim_atlas_barcodes")

    logger.info("\n=== Kim Atlas Result ===")
    logger.info(f"Treatment-naive LUAD max H₁: {max_h1:.4f}")
    logger.info(
        "Compare to drug-treated: "
        "D9 erlotinib=1.77, D14 osimertinib=3.78, PDX residual=6.08"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--subsample", type=int, default=5000)
    args = parser.parse_args()
    main(seed=args.seed, subsample=args.subsample)
