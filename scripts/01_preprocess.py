#!/usr/bin/env python
"""
Step 2: Preprocessing pipeline.

Load raw count matrices → QC → normalize → HVG → PCA → cell cycle scoring.
Outputs: data/processed/pc9_erlotinib_processed.h5ad

Run: python scripts/01_preprocess.py --seed 42
"""

import argparse
import logging
import sys

sys.path.insert(0, "src")


from sctda_plasticity.config import (
    DATA_PROCESSED,
    DATA_RAW,
    DATASETS,
    MAX_GENES,
    MAX_PCT_MITO,
    MIN_CELLS,
    MIN_GENES,
    N_PCS,
    N_TOP_GENES,
    SEED,
)
from sctda_plasticity.data import load_geo_data, preprocess, qc_filter

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(seed: int = SEED):
    logger.info(f"=== Preprocessing Pipeline (seed={seed}) ===")

    # ── Load discovery dataset ──────────────────────────────
    accession = DATASETS["discovery"]["accession"]
    logger.info(f"Loading {accession} from {DATA_RAW}")
    adata = load_geo_data(DATA_RAW, accession)
    logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Ensure timepoint annotation exists
    if "timepoint" not in adata.obs.columns:
        logger.warning(
            "No 'timepoint' column found. Attempting to infer from sample/batch info."
        )
        # Try common column names
        for col in ["sample", "batch", "condition", "time", "group"]:
            if col in adata.obs.columns:
                adata.obs["timepoint"] = adata.obs[col]
                logger.info(f"Using '{col}' as timepoint column")
                break

    if "timepoint" in adata.obs.columns:
        logger.info(f"Timepoints: {adata.obs['timepoint'].value_counts().to_dict()}")

    # ── QC filtering ────────────────────────────────────────
    logger.info("Running QC filtering...")
    adata = qc_filter(
        adata,
        min_genes=MIN_GENES,
        max_genes=MAX_GENES,
        max_pct_mito=MAX_PCT_MITO,
        min_cells=MIN_CELLS,
    )

    # Save QC metrics for plotting
    qc_path = DATA_PROCESSED / "qc_metrics.csv"
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(qc_path)
    logger.info(f"QC metrics saved to {qc_path}")

    # ── Preprocess ──────────────────────────────────────────
    logger.info("Running preprocessing (normalize → HVG → regress → scale → PCA)...")
    batch_key = "timepoint" if "timepoint" in adata.obs.columns else None
    adata = preprocess(
        adata,
        n_top_genes=N_TOP_GENES,
        n_pcs=N_PCS,
        batch_key=batch_key,
        seed=seed,
    )

    # ── Save ────────────────────────────────────────────────
    out_path = DATA_PROCESSED / "pc9_erlotinib_processed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved processed data: {out_path}")
    logger.info(f"Final shape: {adata.n_obs} cells × {adata.n_vars} genes, {N_PCS} PCs")

    # Summary
    logger.info("\n=== Preprocessing Summary ===")
    logger.info(f"  Cells: {adata.n_obs}")
    logger.info(f"  Genes (HVG): {adata.n_vars}")
    logger.info(f"  PCs: {adata.obsm['X_pca'].shape[1]}")
    if "phase" in adata.obs.columns:
        logger.info(f"  Cell cycle: {adata.obs['phase'].value_counts().to_dict()}")
    if "timepoint" in adata.obs.columns:
        logger.info(
            f"  Timepoints: {adata.obs['timepoint'].value_counts().to_dict()}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess PC9 erlotinib scRNA-seq data")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
