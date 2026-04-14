#!/usr/bin/env python
"""
Full Maynard 2020 cohort analysis (S01_datafinal).

27,489 cells across 33 patients with EGFR/ALK/BRAF/KRAS/ROS1/MET drivers.
Filter to EGFR-mutant subset (13,244 cells). Three timepoints per patient
in 'analysis' column: naive (TN) / grouped_pr (PER) / grouped_pd (PD).

This is the definitive patient-level validation of the framework.

Run: python scripts/15_maynard_full_cohort.py --seed 42
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
    MAX_PCT_MITO,
    MIN_CELLS,
    N_PCS,
    N_TOP_GENES,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.topology import compute_persistent_homology
from sctda_plasticity.visualize import (
    plot_persistence_barcodes,
    plot_persistence_comparison,
    save_figure,
    set_publication_style,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_maynard_full(data_dir: Path, driver: str = "EGFR") -> sc.AnnData:
    """Load the full Maynard S01 data, filter to a driver-gene subset."""
    csv_dir = data_dir / "Data_input" / "csv_files"

    logger.info("Loading S01 metadata (27k cells)...")
    meta_raw = pd.read_csv(csv_dir / "S01_metacells.csv", index_col=0)
    # The real cell IDs are in cell_id column, not the index
    meta = meta_raw.set_index("cell_id")
    logger.info(f"Metadata: {meta.shape[0]} cells, {meta['patient_id'].nunique()} patients")
    logger.info(f"Driver breakdown: {meta['driver_gene'].value_counts().to_dict()}")

    # Filter to driver-gene subset
    mask = meta["driver_gene"] == driver
    meta_sub = meta[mask].copy()
    logger.info(f"After {driver} filter: {meta_sub.shape[0]} cells, "
                f"{meta_sub['patient_id'].nunique()} patients")

    # Standardize timepoint labels: naive/grouped_pr/grouped_pd -> TN/PER/PD
    tp_map = {"naive": "TN", "grouped_pr": "PER", "grouped_pd": "PD"}
    meta_sub["timepoint"] = meta_sub["analysis"].map(tp_map)
    logger.info(f"Timepoint distribution: "
                f"{meta_sub['timepoint'].value_counts().to_dict()}")

    # Cell IDs to keep
    keep_cells = set(meta_sub.index.astype(str))
    logger.info(f"Target cells for filtering: {len(keep_cells)}")

    # Load CSV with column filter for efficiency (read only keep_cells)
    # The S01_datafinal.csv has rows=genes, cols=cells (incl. first cell-id column)
    logger.info("Loading S01_datafinal.csv (1.5GB, this may take 3-5 min)...")

    # Read just the header to find which columns to keep
    with open(csv_dir / "S01_datafinal.csv", "r") as f:
        header = f.readline().strip().split(",")
    # header[0] is probably empty or gene col; header[1:] are cell IDs
    all_cells = header[1:]
    col_names_clean = [c.strip('"').strip() for c in all_cells]

    # Match cell IDs
    keep_idx = [i for i, c in enumerate(col_names_clean) if c in keep_cells]
    logger.info(f"Cells matched: {len(keep_idx)} / {len(all_cells)}")

    if len(keep_idx) < 100:
        # Fallback — use index-based
        logger.warning("Few matches; falling back to full load")
        df = pd.read_csv(csv_dir / "S01_datafinal.csv", index_col=0)
    else:
        # Use lambda filter — pandas handles quote-stripping automatically
        keep_set = set(keep_cells)
        logger.info(f"Reading {len(keep_idx)} matching columns...")
        df = pd.read_csv(
            csv_dir / "S01_datafinal.csv", index_col=0,
            usecols=lambda c: c == "" or c in keep_set,
        )

    logger.info(f"Loaded raw: {df.shape[0]} genes x {df.shape[1]} cells")

    # Transpose to cells x genes
    adata = sc.AnnData(df.T)
    adata.var_names_make_unique()

    # Clean cell ID index — remove quotes
    adata.obs_names = [str(x).strip('"').strip() for x in adata.obs_names]

    # Attach metadata
    common = adata.obs_names.intersection(meta_sub.index.astype(str))
    logger.info(f"Intersection: {len(common)} cells")
    adata = adata[list(common)].copy()

    for col in ["patient_id", "timepoint", "analysis", "driver_gene",
                "driver_mutation", "biopsy_type", "biopsy_site",
                "early_treatment_status", "best_response_status",
                "sample_name"]:
        if col in meta_sub.columns:
            adata.obs[col] = meta_sub.loc[adata.obs_names, col].values

    logger.info(f"Final AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    logger.info(f"Timepoints: {adata.obs['timepoint'].value_counts().to_dict()}")
    logger.info(f"Patients: {adata.obs['patient_id'].nunique()}")

    return adata


def smartseq2_preprocess(adata: sc.AnnData, n_top: int = N_TOP_GENES,
                         n_pcs: int = N_PCS, seed: int = 42,
                         batch_key: str = "patient_id") -> sc.AnnData:
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_cells(adata, max_genes=10000)
    adata = adata[adata.obs.pct_counts_mt < MAX_PCT_MITO, :].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    logger.info(f"QC: {n0} -> {adata.n_obs} cells")

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, batch_key=batch_key)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1, adata.n_obs - 1),
              random_state=seed)
    return adata


def main(seed: int = SEED):
    logger.info(f"=== Full Maynard 2020 EGFR cohort analysis (seed={seed}) ===")

    data_dir = DATA_RAW / "maynard2020"
    adata = load_maynard_full(data_dir, driver="EGFR")

    adata = smartseq2_preprocess(adata, seed=seed)
    out_path = DATA_PROCESSED / "maynard2020_egfr_full_processed.h5ad"
    adata.write(out_path)
    logger.info(f"Saved: {out_path}")

    ph_dir = RESULTS_DIR / "ph_maynard_full"
    ph_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.RandomState(seed)

    # ── 1. Pooled EGFR across patients ──────────────────────
    logger.info("\n--- Pooled EGFR across all patients ---")
    results = []
    for tp in ["TN", "PER", "PD"]:
        cells = adata[adata.obs["timepoint"] == tp]
        if cells.n_obs < 30: continue
        if cells.n_obs > 1500:
            cells = cells[rng.choice(cells.n_obs, 1500, replace=False)]
        ph = compute_persistent_homology(cells.obsm["X_pca"], maxdim=PH_MAX_DIM,
                                         n_pcs=min(PH_N_PCS, cells.obsm["X_pca"].shape[1]))
        h1 = ph["max_persistence"].get(1, 0)
        logger.info(f"  pooled {tp}: n={cells.n_obs}, max H1={h1:.3f}")
        results.append({"patient": "pooled_EGFR", "timepoint": tp,
                       "n_cells": cells.n_obs, "max_h1": h1,
                       "n_h1_features": len(ph["dgms"][1])})
        for d in range(PH_MAX_DIM + 1):
            np.save(ph_dir / f"dgm{d}_pooled_{tp}.npy", ph["dgms"][d])

    # ── 2. Per-patient matched trajectories (only patients with >= 2 timepoints) ──
    logger.info("\n--- Per-patient matched trajectories ---")
    pat_tp_counts = adata.obs.groupby(["patient_id", "timepoint"], observed=True).size().reset_index(name="n")
    # Find patients with at least 2 timepoints and at least 30 cells per
    eligible_pivot = pat_tp_counts.pivot(index="patient_id", columns="timepoint",
                                          values="n").fillna(0)
    logger.info(f"\nCell counts per patient × timepoint:")
    logger.info(f"\n{eligible_pivot.to_string()}")

    for pat in eligible_pivot.index:
        row = eligible_pivot.loc[pat]
        has_tps = [tp for tp in ["TN", "PER", "PD"]
                   if tp in row.index and row[tp] >= 30]
        if len(has_tps) < 2:
            continue
        for tp in has_tps:
            cells = adata[(adata.obs["patient_id"] == pat) &
                         (adata.obs["timepoint"] == tp)]
            if cells.n_obs < 30: continue
            if cells.n_obs > 1000:
                cells = cells[rng.choice(cells.n_obs, 1000, replace=False)]
            ph = compute_persistent_homology(cells.obsm["X_pca"], maxdim=PH_MAX_DIM,
                                             n_pcs=min(PH_N_PCS, cells.obsm["X_pca"].shape[1]))
            h1 = ph["max_persistence"].get(1, 0)
            logger.info(f"  {pat} {tp}: n={cells.n_obs}, max H1={h1:.3f}")
            results.append({"patient": pat, "timepoint": tp,
                           "n_cells": cells.n_obs, "max_h1": h1,
                           "n_h1_features": len(ph["dgms"][1])})
            for d in range(PH_MAX_DIM + 1):
                np.save(ph_dir / f"dgm{d}_{pat}_{tp}.npy", ph["dgms"][d])

    df = pd.DataFrame(results)
    df.to_csv(TABLES_DIR / "ph_maynard_full_by_patient.csv", index=False)

    logger.info("\n=== Full summary ===")
    logger.info(df.to_string(index=False))

    # Per-patient trajectory stats
    logger.info("\n=== Per-patient trajectories (where matched) ===")
    for pat in df["patient"].unique():
        if pat == "pooled_EGFR": continue
        pat_df = df[df["patient"] == pat]
        tps = pat_df["timepoint"].tolist()
        vals = pat_df["max_h1"].tolist()
        pairs = dict(zip(tps, vals))
        logger.info(f"  {pat}: {pairs}")

    # Summary stats: TN vs post-treatment (PER or PD) across patients
    logger.info("\n=== Within-patient H1 change (TN -> post-treatment) ===")
    trajectories = []
    for pat in df["patient"].unique():
        if pat == "pooled_EGFR": continue
        sub = df[df["patient"] == pat].set_index("timepoint")["max_h1"].to_dict()
        if "TN" not in sub: continue
        tn = sub["TN"]
        post_vals = [sub[t] for t in ["PER", "PD"] if t in sub]
        if not post_vals: continue
        for post_tp, post_val in [(t, sub[t]) for t in ["PER", "PD"] if t in sub]:
            ratio = post_val / tn if tn > 0 else float("inf")
            logger.info(f"  {pat} TN={tn:.2f} -> {post_tp}={post_val:.2f} (ratio={ratio:.2f}x)")
            trajectories.append({"patient": pat, "TN": tn,
                                "post_tp": post_tp, "post_h1": post_val,
                                "ratio": ratio})

    if trajectories:
        traj_df = pd.DataFrame(trajectories)
        n_up = (traj_df["ratio"] > 1).sum()
        n_down = (traj_df["ratio"] < 1).sum()
        median_ratio = traj_df["ratio"].median()
        mean_ratio = traj_df["ratio"].mean()
        logger.info(f"\nOut of {len(traj_df)} matched TN->post pairs: "
                    f"{n_up} up, {n_down} down, median ratio={median_ratio:.2f}x, "
                    f"mean ratio={mean_ratio:.2f}x")
        traj_df.to_csv(TABLES_DIR / "maynard_full_trajectories.csv", index=False)

    logger.info("\n=== Maynard full EGFR analysis complete ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
