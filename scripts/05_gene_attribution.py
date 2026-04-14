#!/usr/bin/env python
"""
Step 6: Gene attribution on topological features.

Gene connectivity on Mapper graph, DE around loops,
GO/Hallmark enrichment, TF activity analysis.

Run: python scripts/05_gene_attribution.py --seed 42
"""

import argparse
import logging
import pickle
import sys

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.topology import gene_connectivity

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def run_enrichment(gene_list: list, output_dir: str, label: str) -> pd.DataFrame | None:
    """Run gene set enrichment using gseapy."""
    try:
        import gseapy as gp

        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=["MSigDB_Hallmark_2020", "GO_Biological_Process_2021"],
            organism="human",
            outdir=output_dir,
            no_plot=True,
        )
        results = enr.results
        results.to_csv(f"{output_dir}/{label}_enrichment.csv", index=False)
        top_term = results.iloc[0]["Term"] if len(results) > 0 else "none"
        logger.info(f"Enrichment ({label}): {len(results)} terms, top: {top_term}")
        return results
    except Exception as e:
        logger.warning(f"Enrichment failed: {e}. Skipping — may need internet access.")
        return None


def main(seed: int = SEED):
    logger.info(f"=== Gene Attribution Pipeline (seed={seed}) ===")

    # ── Load data ───────────────────────────────────────────
    h5ad_path = DATA_PROCESSED / "pc9_erlotinib_analyzed.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Run standard analysis first: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    enrichment_dir = RESULTS_DIR / "enrichment"
    enrichment_dir.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load Mapper graphs ──────────────────────────────────
    mapper_dir = RESULTS_DIR / "mapper"
    mapper_files = list(mapper_dir.glob("mapper_*.pkl"))
    if not mapper_files:
        raise FileNotFoundError("No Mapper graphs found. Run topology script first.")

    # Use the best Mapper (prefer EMT score filter)
    best_mapper_path = mapper_dir / "mapper_emt_score.pkl"
    if not best_mapper_path.exists():
        best_mapper_path = mapper_files[0]

    logger.info(f"Using Mapper graph: {best_mapper_path.name}")
    with open(best_mapper_path, "rb") as f:
        mapper_graph = pickle.load(f)

    # ── Gene connectivity ───────────────────────────────────
    logger.info("\n--- Gene Connectivity Analysis ---")

    # Use raw expression for gene connectivity
    if adata.raw is not None:
        expr_matrix = adata.raw.X
        gene_names = adata.raw.var_names.values
    else:
        expr_matrix = adata.X
        gene_names = adata.var_names.values

    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    gc = gene_connectivity(
        mapper_graph, np.array(expr_matrix), gene_names, n_top=100
    )

    gc["scores"].to_csv(TABLES_DIR / "gene_connectivity_scores.csv", index=False)
    logger.info(f"Top 20 connected genes:\n{gc['scores'].head(20).to_string()}")

    # ── Enrichment on top connected genes ───────────────────
    logger.info("\n--- Pathway Enrichment (top connected genes) ---")
    top_genes = gc["top_genes"]

    run_enrichment(
        top_genes, str(enrichment_dir), "mapper_top_connected"
    )

    # ── Per-timepoint gene connectivity ─────────────────────
    if "timepoint" in adata.obs.columns:
        logger.info("\n--- Per-Timepoint Gene Connectivity ---")

        timepoints = sorted(adata.obs["timepoint"].unique())
        tp_top_genes = {}

        for tp in timepoints:
            tp_mask = adata.obs["timepoint"] == tp
            tp_indices = np.where(tp_mask)[0]

            # Filter Mapper nodes to this timepoint
            tp_expr = np.array(expr_matrix)[tp_indices]
            tp_gene_names = gene_names

            # Simple approach: compute connectivity using full Mapper but only
            # cells from this timepoint
            tp_connectivity = np.zeros(len(tp_gene_names))
            for node_id, cell_indices in mapper_graph["nodes"].items():
                # Cells in this node AND this timepoint
                overlap = [i for i in cell_indices if i in tp_indices]
                if len(overlap) < 2:
                    continue
                # Remap indices
                local_idx = [np.searchsorted(tp_indices, i) for i in overlap if i in tp_indices]
                if len(local_idx) < 2:
                    continue
                node_expr = tp_expr[local_idx]
                node_mean = np.mean(node_expr, axis=0)
                global_mean = np.mean(tp_expr, axis=0)
                tp_connectivity += np.abs(node_mean - global_mean) * len(local_idx)

            tp_connectivity /= len(tp_indices)
            top_idx = np.argsort(tp_connectivity)[::-1][:50]
            tp_top = [tp_gene_names[i] for i in top_idx]
            tp_top_genes[tp] = tp_top

            logger.info(f"  {tp}: top gene = {tp_top[0]}")

        # Compare gene overlap across timepoints
        logger.info("\n--- Gene Overlap Across Timepoints ---")
        overlap_matrix = {}
        for tp1 in timepoints:
            for tp2 in timepoints:
                s1 = set(tp_top_genes[tp1][:50])
                s2 = set(tp_top_genes[tp2][:50])
                overlap_matrix[(tp1, tp2)] = len(s1 & s2)

        overlap_df = pd.DataFrame(
            {tp: {tp2: overlap_matrix[(tp, tp2)] for tp2 in timepoints} for tp in timepoints}
        )
        overlap_df.to_csv(TABLES_DIR / "gene_overlap_timepoints.csv")
        logger.info(f"\n{overlap_df.to_string()}")

    # ── DE analysis for drug-tolerant vs untreated ──────────
    if "timepoint" in adata.obs.columns:
        logger.info("\n--- Differential Expression (late vs. untreated) ---")

        adata_de = adata.copy()
        # Compare late timepoints (D9, D11) vs untreated (D0)
        late_mask = adata_de.obs["timepoint"].isin(["D9", "D11"])
        early_mask = adata_de.obs["timepoint"] == "D0"

        if late_mask.any() and early_mask.any():
            adata_de.obs["dtp_vs_untreated"] = "other"
            adata_de.obs.loc[late_mask, "dtp_vs_untreated"] = "DTP"
            adata_de.obs.loc[early_mask, "dtp_vs_untreated"] = "untreated"

            adata_sub = adata_de[adata_de.obs["dtp_vs_untreated"] != "other"]
            sc.tl.rank_genes_groups(
                adata_sub, groupby="dtp_vs_untreated", method="wilcoxon",
                reference="untreated", use_raw=True,
            )

            de_results = sc.get.rank_genes_groups_df(adata_sub, group="DTP")
            de_results.to_csv(TABLES_DIR / "de_dtp_vs_untreated.csv", index=False)

            # Enrichment on upregulated DE genes
            sig_up = de_results[
                (de_results["pvals_adj"] < 0.05) & (de_results["logfoldchanges"] > 0.5)
            ]["names"].tolist()

            if sig_up:
                logger.info(f"Significantly upregulated in DTP: {len(sig_up)} genes")
                run_enrichment(sig_up[:200], str(enrichment_dir), "dtp_upregulated")

            sig_down = de_results[
                (de_results["pvals_adj"] < 0.05) & (de_results["logfoldchanges"] < -0.5)
            ]["names"].tolist()

            if sig_down:
                logger.info(f"Significantly downregulated in DTP: {len(sig_down)} genes")
                run_enrichment(sig_down[:200], str(enrichment_dir), "dtp_downregulated")

    logger.info("\n=== Gene Attribution Complete ===")
    logger.info(f"Results in {TABLES_DIR} and {enrichment_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gene attribution on topological features")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
