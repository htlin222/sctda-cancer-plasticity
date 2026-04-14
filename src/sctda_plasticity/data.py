"""
Data loading, QC, and preprocessing for scRNA-seq data.
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

# Tirosh et al. 2016 cell cycle genes
# fmt: off
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]
# fmt: on


def load_geo_data(geo_dir: Path, dataset_id: str) -> sc.AnnData:
    """Load scRNA-seq data from GEO download directory.

    Handles GSE134839 (per-timepoint Drop-seq DGE files) and generic formats.

    Parameters
    ----------
    geo_dir : Path
        Directory containing downloaded GEO files
    dataset_id : str
        GEO accession (e.g., 'GSE134839')

    Returns
    -------
    AnnData with raw counts and timepoint annotation
    """
    data_path = geo_dir / dataset_id
    if not data_path.exists():
        raise FileNotFoundError(f"Data not found at {data_path}. Run `make download` first.")

    # ── GSE134839: per-timepoint DGE files ──────────────────
    dge_files = sorted(data_path.glob("*.dge.txt.gz")) + sorted(data_path.glob("*.dge.txt"))
    if dge_files:
        logger.info(f"Found {len(dge_files)} Drop-seq DGE files")
        adatas = []
        for f in dge_files:
            # Extract timepoint from filename (e.g., GSM3972657_D0.dge.txt.gz → D0)
            stem = f.name.split(".")[0]  # GSM3972657_D0
            tp = stem.split("_", 1)[1] if "_" in stem else stem

            logger.info(f"  Loading {f.name} (timepoint={tp})")
            df = pd.read_csv(f, sep="\t", index_col=0)
            adata_tp = sc.AnnData(df.T)  # genes × cells → cells × genes
            adata_tp.obs["timepoint"] = tp
            adata_tp.obs["sample"] = stem
            # Make cell barcodes unique across timepoints
            adata_tp.obs_names = [f"{tp}_{bc}" for bc in adata_tp.obs_names]
            adatas.append(adata_tp)
            logger.info(f"    {adata_tp.n_obs} cells × {adata_tp.n_vars} genes")

        adata = sc.concat(adatas, join="outer", fill_value=0)
        adata.obs_names_make_unique()
        logger.info(
            f"Combined: {adata.n_obs} cells × {adata.n_vars} genes, "
            f"timepoints: {adata.obs['timepoint'].value_counts().to_dict()}"
        )
        return adata

    # ── h5ad files ──────────────────────────────────────────
    h5ad_files = list(data_path.glob("*.h5ad"))
    if h5ad_files:
        logger.info(f"Loading h5ad: {h5ad_files[0]}")
        return sc.read_h5ad(h5ad_files[0])

    # ── 10x MTX format ──────────────────────────────────────
    mtx_files = list(data_path.glob("**/matrix.mtx*"))
    if mtx_files:
        logger.info(f"Loading 10x-style from: {mtx_files[0].parent}")
        return sc.read_10x_mtx(mtx_files[0].parent)

    # ── CSV/TSV matrices ────────────────────────────────────
    csv_files = sorted(data_path.glob("*.csv*")) + sorted(data_path.glob("*.tsv*"))
    if csv_files:
        logger.info(f"Loading tabular: {csv_files[0]}")
        sep = "\t" if "tsv" in str(csv_files[0]) else ","
        df = pd.read_csv(csv_files[0], index_col=0, sep=sep)
        return sc.AnnData(df.T)

    raise FileNotFoundError(f"No recognized data format in {data_path}")


def qc_filter(
    adata: sc.AnnData,
    min_genes: int = 200,
    max_genes: int = 5000,
    max_pct_mito: float = 20.0,
    min_cells: int = 3,
) -> sc.AnnData:
    """Standard QC filtering.

    Parameters
    ----------
    adata : AnnData
    min_genes, max_genes : int
        Cell filter thresholds
    max_pct_mito : float
        Maximum mitochondrial percentage
    min_cells : int
        Gene filter threshold

    Returns
    -------
    Filtered AnnData (copy)
    """
    adata = adata.copy()
    n_before = adata.n_obs

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mito, :].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_after = adata.n_obs
    logger.info(f"QC: {n_before} → {n_after} cells ({n_before - n_after} removed)")

    return adata


def preprocess(
    adata: sc.AnnData,
    n_top_genes: int = 3000,
    n_pcs: int = 50,
    batch_key: Optional[str] = None,
    seed: int = 42,
) -> sc.AnnData:
    """Full preprocessing: normalize → HVG → scale → PCA → cell cycle score.

    Parameters
    ----------
    adata : AnnData
        QC-filtered data with raw counts
    n_top_genes : int
        Number of highly variable genes
    n_pcs : int
        Number of principal components
    batch_key : str, optional
        Batch key for HVG selection
    seed : int
        Random seed

    Returns
    -------
    Preprocessed AnnData with PCA in .obsm['X_pca']
    """
    adata = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store raw for DE later
    adata.raw = adata

    # HVG
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, batch_key=batch_key
    )
    adata = adata[:, adata.var.highly_variable].copy()

    # Regress and scale
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=n_pcs, random_state=seed)

    # Cell cycle scoring
    s_genes_found = [g for g in S_GENES if g in adata.raw.var_names]
    g2m_genes_found = [g for g in G2M_GENES if g in adata.raw.var_names]
    if s_genes_found and g2m_genes_found:
        sc.tl.score_genes_cell_cycle(
            adata, s_genes=s_genes_found, g2m_genes=g2m_genes_found
        )
        logger.info(
            f"Cell cycle scored: {len(s_genes_found)} S genes, "
            f"{len(g2m_genes_found)} G2M genes"
        )

    logger.info(
        f"Preprocessed: {adata.n_obs} cells × {adata.n_vars} genes, "
        f"{n_pcs} PCs"
    )
    return adata


def remove_cell_cycle_genes(adata: sc.AnnData) -> sc.AnnData:
    """Remove cell cycle genes and return a new AnnData.

    Used for cell cycle ablation experiment.
    """
    cc_genes = set(S_GENES + G2M_GENES)
    mask = ~adata.var_names.isin(cc_genes)
    n_removed = (~mask).sum()
    adata_nocc = adata[:, mask].copy()
    logger.info(f"Removed {n_removed} cell cycle genes, {adata_nocc.n_vars} remaining")
    return adata_nocc
