#!/usr/bin/env python
"""
Validation on Maynard 2020 patient longitudinal data (neo-osi cohort).

Source: czbiohub-sf/scell_lung_adenocarcinoma (Google Drive).
The neo-osi cohort contains 3 EGFR-mutant NSCLC patients with matched
treatment-naive and post-osimertinib biopsies, sequenced by Smart-seq2.

This directly validates our framework on real patient tissue with osimertinib
treatment — the specific request of Round 2 Reviewer 2 that was previously
deferred.

Run: python scripts/14_validation_maynard.py --seed 42
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

sys.path.insert(0, "src")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    DATA_RAW,
    FIGURES_DIR,
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
from sctda_plasticity.topology import (
    compare_persistence_diagrams,
    compute_persistent_homology,
)
from sctda_plasticity.visualize import (
    plot_persistence_barcodes,
    plot_persistence_comparison,
    save_figure,
    set_publication_style,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_maynard_neoosi(data_dir: Path) -> sc.AnnData:
    """Load the neo-osi cohort (3 EGFR patients, Smart-seq2)."""
    csv_dir = data_dir / "Data_input" / "csv_files"

    logger.info("Loading neo-osi metadata...")
    meta = pd.read_csv(csv_dir / "neo-osi_metadata.csv", index_col=0)
    logger.info(f"Metadata: {meta.shape[0]} plates, "
                f"patients: {sorted(meta['patient_id'].unique())}")
    logger.info(f"Treatments: {meta['treatment_status'].value_counts().to_dict()}")

    # Build plate -> metadata map
    plate_to_meta = meta.set_index("plate")[
        ["patient_id", "treatment_status", "sample_name",
         "driver_gene", "driver_mutation", "biopsy_type"]
    ].to_dict(orient="index")

    logger.info(f"Loading raw counts (this is a ~200MB CSV)...")
    raw = pd.read_csv(csv_dir / "neo-osi_rawdata.csv", index_col=0)
    logger.info(f"Raw matrix: {raw.shape[0]} genes x {raw.shape[1]} cells")

    # Transpose: we want cells x genes
    adata = sc.AnnData(raw.T)
    adata.var_names_make_unique()

    # Extract plate from cell name: e.g. "A10_B000561_S10.homo" -> "B000561"
    plates = adata.obs_names.str.split("_").str[1]
    adata.obs["plate"] = plates

    # Map plate to metadata
    for col in ["patient_id", "treatment_status", "sample_name",
                "driver_gene", "driver_mutation", "biopsy_type"]:
        adata.obs[col] = adata.obs["plate"].map(
            lambda p: plate_to_meta.get(p, {}).get(col, "unknown")
        )

    # Normalize treatment_status labels (latin-1 encoding artifact -> "naive")
    adata.obs["treatment_status"] = adata.obs["treatment_status"].str.replace(
        "na\x95ve", "naive", regex=False
    ).str.replace("na.ve", "naive", regex=True)

    # Create timepoint column matching other datasets
    tp_map = {
        "naive": "TN",
        "post-osimertinib": "RD",  # residual disease = post-treatment
    }
    adata.obs["timepoint"] = adata.obs["treatment_status"].map(tp_map)

    # Filter to cells with known plate/treatment
    mask = adata.obs["treatment_status"].isin(["naive", "post-osimertinib"])
    adata = adata[mask].copy()
    logger.info(f"After plate-metadata matching: {adata.n_obs} cells")
    logger.info(f"Treatment: {adata.obs['treatment_status'].value_counts().to_dict()}")
    logger.info(f"Patient: {adata.obs['patient_id'].value_counts().to_dict()}")

    return adata


def smartseq2_preprocess(adata: sc.AnnData, n_top_genes: int = 3000,
                         n_pcs: int = 50, seed: int = 42) -> sc.AnnData:
    """Smart-seq2-tailored preprocessing (higher gene counts than 10x)."""
    adata = adata.copy()

    # QC: Smart-seq2 typically has higher gene count — relax max_genes
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    n_before = adata.n_obs
    # Smart-seq2: min_genes=500 (higher), max_genes=10000 (Smart-seq2 captures more)
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_cells(adata, max_genes=10000)
    adata = adata[adata.obs.pct_counts_mt < MAX_PCT_MITO, :].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    logger.info(f"QC: {n_before} → {adata.n_obs} cells")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e6)  # Smart-seq2 uses TPM-like normalization
    sc.pp.log1p(adata)
    adata.raw = adata

    # HVG — use patient as batch (3 patients)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key="patient_id")
    adata = adata[:, adata.var.highly_variable].copy()

    # Scale and PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1, adata.n_obs - 1),
              random_state=seed)

    logger.info(f"Preprocessed: {adata.n_obs} cells x {adata.n_vars} genes, "
                f"{adata.obsm['X_pca'].shape[1]} PCs")
    return adata


def main(seed: int = SEED):
    logger.info(f"=== Maynard 2020 Patient Validation (seed={seed}) ===")

    data_dir = DATA_RAW / "maynard2020"
    adata = load_maynard_neoosi(data_dir)

    # Preprocess
    logger.info("\n--- Preprocessing (Smart-seq2 tailored) ---")
    adata = smartseq2_preprocess(adata, seed=seed)

    out_path = DATA_PROCESSED / "maynard2020_neoosi_processed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved: {out_path}")

    # Per-patient + per-timepoint PH
    ph_dir = RESULTS_DIR / "ph_maynard"
    ph_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.RandomState(seed)
    results = []
    dgms_by_group = {}

    # Overall TN vs RD (pooled across patients)
    logger.info("\n--- Pooled TN vs RD across 3 EGFR patients ---")
    for tp in ["TN", "RD"]:
        cells = adata[adata.obs["timepoint"] == tp]
        n = cells.n_obs
        if n < 30:
            logger.info(f"  {tp}: only {n} cells, skip")
            continue
        # Use all cells or subsample to 1000
        if n > 1000:
            cells = cells[rng.choice(n, 1000, replace=False)]

        X = cells.obsm["X_pca"]
        ph = compute_persistent_homology(X, maxdim=PH_MAX_DIM, n_pcs=min(PH_N_PCS, X.shape[1]))
        max_h1 = ph["max_persistence"].get(1, 0)
        n_h1 = len(ph["dgms"][1])

        logger.info(f"  {tp} (pooled): n={cells.n_obs}, max H1={max_h1:.3f}, n_h1={n_h1}")
        dgms_by_group[f"pooled_{tp}"] = ph["dgms"]
        results.append({
            "patient": "pooled",
            "timepoint": tp,
            "n_cells": cells.n_obs,
            "max_h1": max_h1,
            "n_h1_features": n_h1,
        })
        for dim in range(PH_MAX_DIM + 1):
            np.save(ph_dir / f"dgm{dim}_pooled_{tp}.npy", ph["dgms"][dim])

    # Per-patient TN vs RD
    logger.info("\n--- Per-patient PH ---")
    for pat in sorted(adata.obs["patient_id"].unique()):
        for tp in ["TN", "RD"]:
            cells = adata[(adata.obs["patient_id"] == pat) &
                         (adata.obs["timepoint"] == tp)]
            if cells.n_obs < 30:
                logger.info(f"  {pat} {tp}: only {cells.n_obs} cells, skip")
                continue
            X = cells.obsm["X_pca"]
            ph = compute_persistent_homology(X, maxdim=PH_MAX_DIM, n_pcs=min(PH_N_PCS, X.shape[1]))
            max_h1 = ph["max_persistence"].get(1, 0)
            n_h1 = len(ph["dgms"][1])
            logger.info(f"  {pat} {tp}: n={cells.n_obs}, max H1={max_h1:.3f}")
            dgms_by_group[f"{pat}_{tp}"] = ph["dgms"]
            results.append({
                "patient": pat, "timepoint": tp,
                "n_cells": cells.n_obs, "max_h1": max_h1,
                "n_h1_features": n_h1,
            })
            for dim in range(PH_MAX_DIM + 1):
                np.save(ph_dir / f"dgm{dim}_{pat}_{tp}.npy", ph["dgms"][dim])

    df = pd.DataFrame(results)
    df.to_csv(TABLES_DIR / "ph_maynard_by_patient.csv", index=False)

    logger.info("\n=== Maynard neo-osi Summary ===")
    logger.info(df.to_string(index=False))

    # Per-patient TN vs RD ratio
    logger.info("\n=== Per-patient H1 change (TN -> RD) ===")
    for pat in sorted(df["patient"].unique()):
        if pat == "pooled": continue
        tn = df[(df["patient"] == pat) & (df["timepoint"] == "TN")]
        rd = df[(df["patient"] == pat) & (df["timepoint"] == "RD")]
        if len(tn) and len(rd):
            t, r = float(tn["max_h1"].iloc[0]), float(rd["max_h1"].iloc[0])
            ratio = r / t if t > 0 else float("inf")
            arrow = "↑" if r > t else "↓"
            logger.info(f"  {pat}: {t:.3f} → {r:.3f}  (ratio={ratio:.2f}x) {arrow}")

    # Pooled comparison
    pooled_tn = df[(df["patient"] == "pooled") & (df["timepoint"] == "TN")]
    pooled_rd = df[(df["patient"] == "pooled") & (df["timepoint"] == "RD")]
    if len(pooled_tn) and len(pooled_rd):
        t = float(pooled_tn["max_h1"].iloc[0])
        r = float(pooled_rd["max_h1"].iloc[0])
        logger.info(f"\n=== POOLED PATIENTS (3 EGFR): TN={t:.3f} → RD={r:.3f}, "
                    f"ratio={r/t:.2f}x ===")

    # Plot barcodes
    if dgms_by_group:
        fig = plot_persistence_barcodes(
            dgms_by_group, max_dim=PH_MAX_DIM,
            title="Maynard 2020 neo-osi (patient EGFR + osimertinib)",
        )
        save_figure(fig, "figS_maynard_barcodes")

        # H1 diagram overlay
        h1_dgms = {k: v[1] for k, v in dgms_by_group.items() if len(v[1]) > 0}
        if h1_dgms:
            fig = plot_persistence_comparison(h1_dgms, dim=1,
                                              figsize=(110/25.4, 85/25.4))
            save_figure(fig, "figS_maynard_h1_diagrams")

    logger.info("\n=== Maynard Validation Complete ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
