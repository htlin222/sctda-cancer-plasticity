#!/usr/bin/env python
"""
Step 7: Validation on GSE150949 (PC9 osimertinib).

Repeat preprocessing + PH pipeline on validation dataset,
compare topological signatures between erlotinib and osimertinib.

Run: python scripts/06_validation.py --seed 42
"""

import argparse
import logging
import sys

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    DATA_RAW,
    MAX_GENES,
    MAX_PCT_MITO,
    MIN_CELLS,
    MIN_GENES,
    N_PCS,
    N_PERMUTATIONS,
    N_TOP_GENES,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.data import load_geo_data, preprocess, qc_filter
from sctda_plasticity.topology import (
    compare_persistence_diagrams,
    compute_persistent_homology,
    permutation_test_h1,
)
from sctda_plasticity.visualize import (
    plot_permutation_test,
    plot_persistence_comparison,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(seed: int = SEED):
    logger.info(f"=== Validation Pipeline — GSE150949 (seed={seed}) ===")

    # ── Load validation dataset ─────────────────────────────
    accession = "GSE150949"
    logger.info(f"Loading {accession} from {DATA_RAW}")

    data_dir = DATA_RAW / accession
    if not data_dir.exists():
        raise FileNotFoundError("Data not found. Run: make download")

    # Try to load the PC9 count matrix
    count_file = data_dir / "GSE150949_pc9_count_matrix.csv.gz"
    meta_file = data_dir / "GSE150949_metaData_with_lineage.txt.gz"

    if count_file.exists():
        # Load metadata first for subsetting
        meta = None
        if meta_file.exists():
            meta = pd.read_csv(meta_file, sep="\t", index_col=0)
            logger.info(f"Metadata: {meta.shape[0]} cells, timepoints: "
                        f"{meta['time_point'].value_counts().to_dict()}")

        # Subsample cells to keep memory manageable (56k cells → ~5k)
        # Stratified sample by timepoint
        if meta is not None and meta.shape[0] > 5000:
            rng = np.random.RandomState(seed)
            sample_idx = []
            for tp in meta["time_point"].unique():
                tp_cells = meta[meta["time_point"] == tp].index.tolist()
                n_sample = min(len(tp_cells), max(500, 5000 // meta["time_point"].nunique()))
                sample_idx.extend(rng.choice(tp_cells, size=n_sample, replace=False).tolist())
            keep_cells = sample_idx
            logger.info(f"Subsampled: {len(keep_cells)} cells (from {meta.shape[0]})")
        else:
            keep_cells = None

        logger.info(f"Loading count matrix: {count_file}")
        if keep_cells is not None:
            # Use lambda usecols filter to select subset of columns
            keep_set = set(keep_cells)
            logger.info(f"Reading {len(keep_set)} cells from CSV...")
            df = pd.read_csv(
                count_file, index_col=0,
                usecols=lambda col: col == "" or col in keep_set,
            )
        else:
            df = pd.read_csv(count_file, index_col=0)

        adata = sc.AnnData(df.T)  # genes x cells → cells x genes
        logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

        # Add metadata
        if meta is not None:
            common_cells = adata.obs_names.intersection(meta.index)
            if len(common_cells) > 0:
                adata = adata[common_cells].copy()
                for col in meta.columns:
                    adata.obs[col] = meta.loc[common_cells, col]
                # Map time_point to timepoint for consistency
                adata.obs["timepoint"] = "D" + adata.obs["time_point"].astype(str)
                logger.info(f"Metadata loaded: {adata.obs['timepoint'].value_counts().to_dict()}")
    else:
        # Fallback to generic loader
        adata = load_geo_data(DATA_RAW, accession)

    logger.info(f"Validation data: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── Preprocess ──────────────────────────────────────────
    logger.info("Preprocessing validation data...")
    adata = qc_filter(adata, min_genes=MIN_GENES, max_genes=MAX_GENES,
                       max_pct_mito=MAX_PCT_MITO, min_cells=MIN_CELLS)
    adata = preprocess(adata, n_top_genes=N_TOP_GENES, n_pcs=N_PCS, seed=seed)

    val_path = DATA_PROCESSED / "pc9_osimertinib_processed.h5ad"
    adata.write(val_path)
    logger.info(f"Saved: {val_path}")

    # ── Persistent homology ─────────────────────────────────
    logger.info("\n--- Persistent Homology (validation) ---")

    val_ph_dir = RESULTS_DIR / "ph_validation"
    val_ph_dir.mkdir(parents=True, exist_ok=True)

    X_pca = adata.obsm["X_pca"]
    ph = compute_persistent_homology(X_pca, maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS)

    for dim in range(PH_MAX_DIM + 1):
        np.save(val_ph_dir / f"dgm{dim}_osimertinib.npy", ph["dgms"][dim])

    max_h1 = ph["max_persistence"].get(1, 0.0)
    logger.info(f"Validation max H₁ persistence: {max_h1:.4f}")

    # ── Permutation test (skip with --skip-permtest for speed) ──
    skip_perm = "--skip-permtest" in sys.argv
    if max_h1 > 0 and not skip_perm:
        logger.info("\n--- Permutation Test (validation) ---")
        X_scaled = adata.X
        if hasattr(X_scaled, "toarray"):
            X_scaled = X_scaled.toarray()
        X_scaled = np.array(X_scaled, dtype=np.float64)

        perm = permutation_test_h1(
            X_scaled, observed_max_h1=max_h1,
            n_permutations=N_PERMUTATIONS, n_pcs=PH_N_PCS, seed=seed,
        )

        fig = plot_permutation_test(
            perm["null_distribution"], perm["observed"], perm["p_value"]
        )
        save_figure(fig, "fig6_permtest_osimertinib")
        logger.info(f"Validation H₁ p-value: {perm['p_value']:.4f}")

    # ── Cross-dataset comparison ────────────────────────────
    logger.info("\n--- Cross-Dataset Comparison ---")

    # Load discovery PH results
    discovery_ph_dir = RESULTS_DIR / "ph"
    discovery_dgms = {}
    for dgm_file in sorted(discovery_ph_dir.glob("dgm1_*.npy")):
        tp = dgm_file.stem.replace("dgm1_", "")
        discovery_dgms[tp] = np.load(dgm_file)

    if discovery_dgms and max_h1 > 0:
        val_dgm1 = ph["dgms"][1]

        # Compare validation vs each discovery timepoint
        distances = {}
        for tp, dgm in discovery_dgms.items():
            if len(dgm) > 0 and len(val_dgm1) > 0:
                d = compare_persistence_diagrams(dgm, val_dgm1, metric="wasserstein")
                distances[tp] = d
                logger.info(f"  Wasserstein(osimertinib, {tp}) = {d:.4f}")

        if distances:
            pd.Series(distances).to_csv(
                TABLES_DIR / "wasserstein_validation_vs_discovery.csv"
            )

        # Overlay persistence diagrams
        erl_dgms = {f"erl_{tp}": dgm for tp, dgm in discovery_dgms.items()}
        all_dgms = {**erl_dgms, "osimertinib": val_dgm1}
        fig = plot_persistence_comparison(all_dgms, dim=1)
        save_figure(fig, "fig6_a_persistence_comparison")

    logger.info("\n=== Validation Complete ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validation on GSE150949")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--skip-permtest", action="store_true", help="Skip permutation test")
    args = parser.parse_args()
    if args.skip_permtest:
        import sctda_plasticity.config as cfg
        cfg.N_PERMUTATIONS = 0
    main(seed=args.seed)
