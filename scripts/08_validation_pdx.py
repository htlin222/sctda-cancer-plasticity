#!/usr/bin/env python
"""
Validation on ASCL1 PDX (GSE243562) — osimertinib pre/post.

4 samples: YU-003 untreated/residual, YU-006 untreated/residual.
10x Cell Ranger triplet format.

Run: python scripts/08_validation_pdx.py --seed 42
"""

import argparse
import logging
import sys
import tarfile
from pathlib import Path

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
    N_TOP_GENES,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.data import preprocess, qc_filter
from sctda_plasticity.topology import (
    compare_persistence_diagrams,
    compute_persistent_homology,
)
from sctda_plasticity.visualize import (
    plot_persistence_barcodes,
    plot_persistence_comparison,
    plot_wasserstein_heatmap,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Sample mapping from GSE243562 series matrix
# GSM titles encode: GSMxxxxx_<number> = <PDX>_<condition>
SAMPLE_MAP = {
    "GSM7790998": ("YU-003", "untreated"),
    "GSM7790999": ("YU-003", "residual"),
    "GSM7791000": ("YU-006", "untreated"),
    "GSM7791001": ("YU-006", "residual"),
}


def extract_tar(tar_path: Path, extract_to: Path) -> list[Path]:
    """Extract GSE RAW.tar and return list of sample directories."""
    extract_to.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        tf.extractall(extract_to)

    # Files are named like GSM7790998_387_barcodes.tsv.gz
    # We need to organize them into per-sample subdirectories
    sample_dirs = []
    for gsm_id in SAMPLE_MAP:
        sample_dir = extract_to / gsm_id
        sample_dir.mkdir(exist_ok=True)
        # Find files starting with this GSM
        for f in extract_to.glob(f"{gsm_id}_*"):
            if f.is_file():
                # Rename to standard 10x Cell Ranger names
                if "barcodes" in f.name:
                    f.rename(sample_dir / "barcodes.tsv.gz")
                elif "features" in f.name or "genes" in f.name:
                    f.rename(sample_dir / "features.tsv.gz")
                elif "matrix" in f.name:
                    f.rename(sample_dir / "matrix.mtx.gz")
        if (sample_dir / "matrix.mtx.gz").exists():
            sample_dirs.append(sample_dir)

    return sample_dirs


def load_pdx_data(data_dir: Path) -> sc.AnnData:
    """Load all PDX samples and concatenate."""
    tar_path = data_dir / "GSE243562_RAW.tar"
    extract_dir = data_dir / "extracted"

    # Extract if needed
    if not extract_dir.exists() or not any(extract_dir.iterdir()):
        logger.info(f"Extracting {tar_path}...")
        sample_dirs = extract_tar(tar_path, extract_dir)
    else:
        sample_dirs = [
            d for d in extract_dir.iterdir()
            if d.is_dir() and (d / "matrix.mtx.gz").exists()
        ]

    logger.info(f"Found {len(sample_dirs)} sample directories")

    adatas = []
    for sample_dir in sorted(sample_dirs):
        gsm_id = sample_dir.name
        pdx, condition = SAMPLE_MAP.get(gsm_id, (gsm_id, "unknown"))
        logger.info(f"Loading {gsm_id} → {pdx} {condition}")

        adata_s = sc.read_10x_mtx(sample_dir)
        adata_s.obs["pdx"] = pdx
        adata_s.obs["condition"] = condition
        adata_s.obs["sample"] = gsm_id
        # Use condition as the "timepoint" for consistency with other scripts
        adata_s.obs["timepoint"] = condition
        # Make barcodes unique
        adata_s.obs_names = [f"{gsm_id}_{bc}" for bc in adata_s.obs_names]
        adata_s.var_names_make_unique()
        adatas.append(adata_s)
        logger.info(f"  {adata_s.n_obs} cells × {adata_s.n_vars} genes")

    adata = sc.concat(adatas, join="outer", fill_value=0)
    adata.obs_names_make_unique()
    logger.info(f"Combined: {adata.n_obs} cells × {adata.n_vars} genes")

    # Filter mouse cells (PDX contamination) — remove cells with > 10% mouse genes
    # Mouse gene symbols are lowercase. In 10x, genes are typically all human or all mouse.
    # Here we check gene naming convention — if features.tsv was human reference, genes are human.
    # Mouse genes would be present as mm10 ensembl prefixes or lowercase.
    mouse_genes = adata.var_names[adata.var_names.str.match(r"^mm10.*|^[a-z][a-z0-9]+$", na=False)]
    if len(mouse_genes) > 0 and len(mouse_genes) < adata.n_vars * 0.5:
        logger.info(f"Detected {len(mouse_genes)} possible mouse genes, filtering them out")
        adata = adata[:, ~adata.var_names.isin(mouse_genes)].copy()

    return adata


def main(seed: int = SEED):
    logger.info(f"=== PDX Validation Pipeline — GSE243562 (seed={seed}) ===")

    data_dir = DATA_RAW / "GSE243562"
    if not data_dir.exists():
        raise FileNotFoundError(f"Data not found at {data_dir}")

    # ── Load ────────────────────────────────────────────────
    adata = load_pdx_data(data_dir)
    logger.info(f"Loaded PDX data: {adata.n_obs} cells × {adata.n_vars} genes")
    logger.info(f"Samples: {adata.obs['sample'].value_counts().to_dict()}")
    logger.info(f"Conditions: {adata.obs['condition'].value_counts().to_dict()}")

    # ── QC + Preprocess ─────────────────────────────────────
    logger.info("\n--- QC + Preprocessing ---")
    adata = qc_filter(
        adata, min_genes=MIN_GENES, max_genes=MAX_GENES,
        max_pct_mito=MAX_PCT_MITO, min_cells=MIN_CELLS,
    )
    logger.info(f"After QC: {adata.n_obs} cells")
    logger.info(f"Conditions: {adata.obs['condition'].value_counts().to_dict()}")

    adata = preprocess(
        adata, n_top_genes=N_TOP_GENES, n_pcs=N_PCS,
        batch_key="pdx", seed=seed,
    )

    out_path = DATA_PROCESSED / "pdx_ascl1_osimertinib_processed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved: {out_path}")

    # ── Persistent Homology per condition ──────────────────
    logger.info("\n--- Per-Condition Persistent Homology ---")

    ph_dir = RESULTS_DIR / "ph_validation_pdx"
    ph_dir.mkdir(parents=True, exist_ok=True)

    # Combined PH (all cells)
    ph_all = compute_persistent_homology(
        adata.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
    )
    logger.info(f"All PDX cells: max H₁ = {ph_all['max_persistence'].get(1, 0):.4f}")

    # Per PDX + condition
    results = []
    dgms_by_group = {}
    for pdx in sorted(adata.obs["pdx"].unique()):
        for condition in sorted(adata.obs["condition"].unique()):
            mask = (adata.obs["pdx"] == pdx) & (adata.obs["condition"] == condition)
            cells = adata[mask]
            if cells.n_obs < 30:
                logger.info(f"  {pdx} {condition}: only {cells.n_obs} cells, skip")
                continue

            ph = compute_persistent_homology(
                cells.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
            )
            max_h1 = ph["max_persistence"].get(1, 0)
            n_h1 = len(ph["dgms"][1]) if PH_MAX_DIM >= 1 else 0

            logger.info(
                f"  {pdx} {condition}: {cells.n_obs} cells, max H₁ = {max_h1:.4f}, {n_h1} features"
            )

            key = f"{pdx}_{condition}"
            dgms_by_group[key] = ph["dgms"]
            results.append({
                "pdx": pdx,
                "condition": condition,
                "n_cells": cells.n_obs,
                "max_h0": ph["max_persistence"].get(0, 0),
                "max_h1": max_h1,
                "n_h1": n_h1,
            })

            for dim in range(PH_MAX_DIM + 1):
                np.save(ph_dir / f"dgm{dim}_{key}.npy", ph["dgms"][dim])

    df = pd.DataFrame(results)
    df.to_csv(TABLES_DIR / "ph_pdx_by_group.csv", index=False)

    # ── Plot barcodes ───────────────────────────────────────
    if dgms_by_group:
        fig = plot_persistence_barcodes(
            dgms_by_group, max_dim=PH_MAX_DIM,
            title="ASCL1 PDX — Osimertinib Pre/Post",
        )
        save_figure(fig, "fig7_a_pdx_barcodes")

        # H1 overlay
        h1_dgms = {k: v[1] for k, v in dgms_by_group.items() if len(v[1]) > 0}
        if h1_dgms:
            fig = plot_persistence_comparison(h1_dgms, dim=1)
            save_figure(fig, "fig7_b_pdx_h1_diagrams")

        # Wasserstein heatmap
        keys = list(h1_dgms.keys())
        if len(keys) > 1:
            n = len(keys)
            W = np.zeros((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    W[i, j] = compare_persistence_diagrams(
                        h1_dgms[keys[i]], h1_dgms[keys[j]], metric="wasserstein",
                    )
                    W[j, i] = W[i, j]
            fig = plot_wasserstein_heatmap(W, keys)
            save_figure(fig, "fig7_c_pdx_wasserstein")
            pd.DataFrame(W, index=keys, columns=keys).to_csv(
                TABLES_DIR / "wasserstein_pdx.csv",
            )

    # ── Summary ─────────────────────────────────────────────
    logger.info("\n=== PDX Validation Summary ===")
    for r in results:
        logger.info(
            f"  {r['pdx']} {r['condition']}: {r['n_cells']} cells, "
            f"max H₁ = {r['max_h1']:.3f}"
        )

    # Critical comparison: untreated vs residual per PDX
    logger.info("\n=== Critical Comparison (untreated vs residual) ===")
    for pdx in sorted(df["pdx"].unique()):
        pdx_data = df[df["pdx"] == pdx]
        untr = pdx_data[pdx_data["condition"] == "untreated"]
        resd = pdx_data[pdx_data["condition"] == "residual"]
        if len(untr) > 0 and len(resd) > 0:
            ratio = resd["max_h1"].iloc[0] / untr["max_h1"].iloc[0]
            logger.info(
                f"  {pdx}: untreated={untr['max_h1'].iloc[0]:.3f}, "
                f"residual={resd['max_h1'].iloc[0]:.3f}, "
                f"ratio={ratio:.2f}x {'↑' if ratio > 1 else '↓'}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
