#!/usr/bin/env python
"""
Gene attribution on PDX data.

Run Mapper + gene connectivity on ASCL1 PDX to validate EMT dominance.

Run: python scripts/11_pdx_gene_attribution.py --seed 42
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
    MAPPER_N_CUBES,
    MAPPER_OVERLAP,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.topology import gene_connectivity, run_mapper

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Same EMT gene set as standard analysis
EMT_GENES = [
    "VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1",
    "MMP2", "MMP9", "CDH1", "OCLN", "TJP1", "DSP", "KRT19", "KRT8",
    "SERPINE1", "TGFBI", "SPARC", "COL1A1", "COL3A1", "COL5A2",
    "ACTA2", "TAGLN", "MYL9", "CALD1", "LGALS1", "TPM1", "TPM2",
    "FSTL1", "ITGB1", "LAMC2", "LOXL2",
]


def main(seed: int = SEED):
    logger.info(f"=== PDX Gene Attribution (seed={seed}) ===")

    adata = sc.read_h5ad(DATA_PROCESSED / "pdx_ascl1_osimertinib_processed.h5ad")
    logger.info(f"Loaded PDX: {adata.n_obs} cells × {adata.n_vars} genes")

    # Compute EMT score on PDX (may need cleanup due to mouse gene prefix)
    # Strip hg38_ prefix from gene names if present
    logger.info(f"First few gene names: {list(adata.var_names[:5])}")
    logger.info(f"First few raw gene names: {list(adata.raw.var_names[:5])}")

    # Check if hg38_ prefix present — need to reconstruct AnnData since raw is immutable
    has_prefix = (
        adata.raw is not None and str(adata.raw.var_names[0]).startswith("hg38_")
    )
    if has_prefix:
        logger.info("Stripping hg38_ prefix from gene names")
        # Create a fresh AnnData with stripped names, keep raw
        import anndata as ad

        new_var_names = adata.var_names.str.replace("^hg38_", "", regex=True)
        raw_new_var_names = adata.raw.var_names.str.replace("^hg38_", "", regex=True)

        # Rebuild raw
        raw_new = ad.AnnData(X=adata.raw.X, var=pd.DataFrame(index=raw_new_var_names))
        # Rebuild main
        adata_new = adata.copy()
        adata_new.var_names = new_var_names
        adata_new.raw = raw_new
        adata = adata_new

    # Build EMT score — check both var and raw.var_names
    raw_genes = set(adata.raw.var_names) if adata.raw is not None else set()
    main_genes = set(adata.var_names)
    emt_found = [g for g in EMT_GENES if g in raw_genes or g in main_genes]
    logger.info(f"EMT genes found: {len(emt_found)}/{len(EMT_GENES)}")
    use_raw_for_score = adata.raw is not None and any(g in raw_genes for g in emt_found)
    if emt_found:
        sc.tl.score_genes(
            adata, gene_list=emt_found, score_name="emt_score",
            use_raw=use_raw_for_score,
        )
        logger.info(
            f"EMT score: mean={adata.obs['emt_score'].mean():.3f}, "
            f"by condition: {adata.obs.groupby('condition')['emt_score'].mean().to_dict()}"
        )

    # Subsample for Mapper (Mapper on 18k cells would be too slow)
    rng = np.random.RandomState(seed)
    MAX_MAPPER = 3000
    if adata.n_obs > MAX_MAPPER:
        idx = rng.choice(adata.n_obs, MAX_MAPPER, replace=False)
        adata_map = adata[idx].copy()
    else:
        adata_map = adata

    # ── Run Mapper with EMT filter ──────────────────────────
    mapper_dir = RESULTS_DIR / "mapper_pdx"
    mapper_dir.mkdir(parents=True, exist_ok=True)

    X_pca = adata_map.obsm["X_pca"][:, :PH_N_PCS]

    if "emt_score" not in adata_map.obs.columns:
        # Fallback: use PC1
        filter_vals = adata_map.obsm["X_pca"][:, 0]
        filter_name = "pc1"
    else:
        filter_vals = adata_map.obs["emt_score"].values
        filter_name = "emt_score"

    logger.info(f"Running Mapper with filter={filter_name}")
    result = run_mapper(
        X_pca, filter_vals, n_cubes=MAPPER_N_CUBES, overlap=MAPPER_OVERLAP,
    )
    n_nodes = len(result["graph"]["nodes"])
    n_edges = sum(len(v) for v in result["graph"]["links"].values())
    logger.info(f"Mapper: {n_nodes} nodes, {n_edges} edges")

    if n_nodes == 0:
        logger.error("Mapper produced 0 nodes! Cannot compute gene connectivity.")
        return

    import pickle

    with open(mapper_dir / f"mapper_pdx_{filter_name}.pkl", "wb") as f:
        pickle.dump(result["graph"], f)
    result["mapper"].visualize(
        result["graph"],
        path_html=str(mapper_dir / f"mapper_pdx_{filter_name}.html"),
        title=f"PDX Mapper ({filter_name})",
    )

    # ── Gene connectivity ──────────────────────────────────
    logger.info("\n--- Gene Connectivity ---")
    if adata_map.raw is not None:
        expr_matrix = adata_map.raw.X
        gene_names = adata_map.raw.var_names.values
    else:
        expr_matrix = adata_map.X
        gene_names = adata_map.var_names.values

    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()
    expr_matrix = np.array(expr_matrix)

    gc = gene_connectivity(result["graph"], expr_matrix, gene_names, n_top=100)
    gc["scores"].to_csv(TABLES_DIR / "pdx_gene_connectivity_scores.csv", index=False)

    logger.info("\nTop 20 Mapper-connected genes (PDX):")
    top20 = gc["scores"].head(20)
    logger.info(f"\n{top20.to_string()}")

    # Check overlap with EMT genes
    top100 = set(gc["scores"].head(100)["gene"].tolist())
    emt_in_top = top100.intersection(EMT_GENES)
    logger.info(f"\nEMT genes in PDX top-100 connected: {len(emt_in_top)}/{len(EMT_GENES)}")
    logger.info(f"Overlap: {sorted(emt_in_top)}")

    # Compare with cell line top genes
    try:
        cell_line_scores = pd.read_csv(TABLES_DIR / "gene_connectivity_scores.csv")
        cl_top100 = set(cell_line_scores.head(100)["gene"].tolist())
        overlap = top100.intersection(cl_top100)
        logger.info(f"\nTop-100 overlap PDX vs cell line: {len(overlap)}/100 genes")
        logger.info(f"Shared key genes: {sorted(overlap)[:30]}")

        # Hypergeometric test for overlap
        from scipy.stats import hypergeom

        # Total genes in both analyses
        all_genes_pdx = set(gc["scores"]["gene"])
        all_genes_cl = set(cell_line_scores["gene"])
        universe = all_genes_pdx & all_genes_cl

        k = len(overlap)  # observed
        M = len(universe)  # population
        n = len(top100 & universe)  # success in pop
        N = len(cl_top100 & universe)  # draws

        pval = hypergeom.sf(k - 1, M, n, N)
        logger.info(f"Hypergeometric test: k={k}, M={M}, n={n}, N={N}, p={pval:.3e}")
    except FileNotFoundError:
        logger.warning("Cell line gene connectivity not found for comparison")

    # ── DE: residual vs untreated ──────────────────────────
    logger.info("\n--- DE: residual vs untreated (PDX) ---")
    sc.tl.rank_genes_groups(
        adata, groupby="condition", method="wilcoxon",
        reference="untreated", use_raw=True,
    )
    de_df = sc.get.rank_genes_groups_df(adata, group="residual")
    de_df.to_csv(TABLES_DIR / "pdx_de_residual_vs_untreated.csv", index=False)

    sig_up = de_df[(de_df["pvals_adj"] < 0.05) & (de_df["logfoldchanges"] > 0.5)]
    sig_down = de_df[(de_df["pvals_adj"] < 0.05) & (de_df["logfoldchanges"] < -0.5)]
    logger.info(f"Significantly up in residual: {len(sig_up)}")
    logger.info(f"Significantly down in residual: {len(sig_down)}")

    # Enrichment
    try:
        import gseapy as gp

        enr_dir = RESULTS_DIR / "enrichment_pdx"
        enr_dir.mkdir(parents=True, exist_ok=True)
        if len(sig_up) > 10:
            enr = gp.enrichr(
                gene_list=sig_up["names"].head(200).tolist(),
                gene_sets=["MSigDB_Hallmark_2020"],
                organism="human",
                outdir=str(enr_dir),
                no_plot=True,
            )
            enr.results.to_csv(enr_dir / "pdx_residual_up_hallmark.csv", index=False)
            logger.info(
                f"Top Hallmark pathway (up): {enr.results.iloc[0]['Term']}"
            )
    except Exception as e:
        logger.warning(f"Enrichment skipped: {e}")

    logger.info("\n=== PDX Gene Attribution Complete ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
