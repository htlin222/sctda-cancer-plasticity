"""
27 — Are the drug-emergent ciliated cells in PDX YU-006 residual tumour-derived
(transdifferentiated) or a non-malignant population?

CNV inference (infercnvpy). Tumour cells carry copy-number aberrations; normal
cells are flat. If the ciliated cells share the bulk tumour's CNV profile, they are
tumour-lineage (transdifferentiation); if CNV-flat, they are normal.

Decision: compare (i) per-cell CNV burden and (ii) correlation of each cell's CNV
profile to the bulk-tumour consensus, between ciliated-high and the tumour bulk.
Run in conda `sctda` env.
"""
from __future__ import annotations
import logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
import infercnvpy as cnv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("cnv")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
CILIA = ["PIFO","RSPH1","CFAP126","CFAP45","CFAP157","CAPS","CAPSL","CETN2","TPPP3",
         "C20orf85","C1orf194","C9orf116","WDR38","LRRC23","FAM81B","MORN2","TUBA1A","C5orf49"]


def gene_positions(symbols):
    from pybiomart import Server
    ds = Server(host="http://www.ensembl.org").marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
    df = ds.query(attributes=["external_gene_name", "chromosome_name", "start_position", "end_position"])
    df.columns = ["gene", "chromosome", "start", "end"]
    df = df[df["chromosome"].isin([str(i) for i in range(1, 23)] + ["X", "Y"])]
    df["chromosome"] = "chr" + df["chromosome"].astype(str)
    df = df.dropna().drop_duplicates("gene").set_index("gene")
    return df


def main():
    P = ad.read_h5ad("data/processed/pdx_ascl1_osimertinib_processed.h5ad")
    sub = P[(P.obs.pdx == "YU-006")].copy()           # untreated + residual
    raw = sub.raw.to_adata()
    raw.var_names = [v.replace("hg38_", "") for v in raw.var_names]
    raw.var_names_make_unique()
    raw.obs["condition"] = sub.obs["condition"].values
    sc.pp.normalize_total(raw, target_sum=1e4); sc.pp.log1p(raw)
    sc.tl.score_genes(raw, [g for g in CILIA if g in raw.var_names], score_name="cilia")
    raw.obs["group"] = np.where(raw.obs["cilia"] > 0.5, "ciliated", "bulk")
    log.info("ciliated-high cells: %d / %d", (raw.obs.group == "ciliated").sum(), raw.n_obs)

    log.info("fetching gene positions from Ensembl ...")
    pos = gene_positions(list(raw.var_names))
    keep = [g for g in raw.var_names if g in pos.index]
    raw = raw[:, keep].copy()
    raw.var["chromosome"] = pos.loc[keep, "chromosome"].values
    raw.var["start"] = pos.loc[keep, "start"].astype(int).values
    raw.var["end"] = pos.loc[keep, "end"].astype(int).values
    log.info("genes with positions: %d", raw.n_vars)

    # infercnv with the whole population mean as baseline (no labelled normal)
    cnv.tl.infercnv(raw, window_size=100, step=10)
    Xc = raw.obsm["X_cnv"]
    Xc = Xc.toarray() if hasattr(Xc, "toarray") else np.asarray(Xc)

    # per-cell CNV burden = mean |signal|
    burden = np.abs(Xc).mean(1)
    raw.obs["cnv_burden"] = burden
    # bulk-tumour consensus (residual bulk, non-ciliated) and correlation of each cell to it
    bulk_mask = (raw.obs.group == "bulk").to_numpy()
    consensus = Xc[bulk_mask].mean(0)
    corr = np.array([np.corrcoef(Xc[i], consensus)[0, 1] for i in range(Xc.shape[0])])
    raw.obs["cnv_corr_to_tumour"] = corr

    g = raw.obs.groupby("group")
    print("\n==== CNV origin test: ciliated vs tumour bulk (YU-006) ====")
    print(g[["cnv_burden", "cnv_corr_to_tumour"]].median())
    print("\nby condition x group (residual is where ciliated emerge):")
    print(raw.obs.groupby(["condition", "group"])[["cnv_burden", "cnv_corr_to_tumour"]].median())
    cil = raw.obs[raw.obs.group == "ciliated"]
    print(f"\nciliated cells: n={len(cil)}  median CNV burden={cil.cnv_burden.median():.4f}  "
          f"median corr-to-tumour={cil.cnv_corr_to_tumour.median():.3f}")
    print(f"bulk tumour:   median CNV burden={raw.obs[bulk_mask].cnv_burden.median():.4f}")
    verdict = ("TUMOUR-DERIVED (transdifferentiation): ciliated cells share tumour CNVs"
               if cil.cnv_corr_to_tumour.median() > 0.3 and cil.cnv_burden.median() > 0.5 * raw.obs[bulk_mask].cnv_burden.median()
               else "AMBIGUOUS / possibly normal: ciliated CNV profile differs from tumour")
    print(f"\nVERDICT: {verdict}")
    raw.obs[["condition","group","cilia","cnv_burden","cnv_corr_to_tumour"]].to_csv(OUT/"pdx_cnv_origin.csv")


if __name__ == "__main__":
    main()
