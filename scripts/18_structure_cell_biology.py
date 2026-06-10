"""
18 — What are the cells that create the drug-induced non-Gaussian topological
structure? (biology of the rescued signal)

Strategy: at osi D14 (where real max-H1 exceeds the dispersion null), Leiden-
cluster the cells, then leave-one-cluster-out: the cluster whose removal most
reduces max-H1 is the structure-driving population. Recover gene symbols from the
raw count CSV (the processed h5ad lost them) and run differential expression of
the structure-driving cells vs the rest, controlling for count depth, to test
whether they are a coherent drug-tolerant/EMT/quiescence state.

Run in conda `sctda` env.
"""
from __future__ import annotations
import logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from ripser import ripser

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("bio")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
CSV = "data/raw/GSE150949/GSE150949_pc9_count_matrix.csv.gz"

MARKERS = {
 "EMT": ["VIM","FN1","SERPINE1","TGFBI","ZEB1","ZEB2","SNAI1","SNAI2","CDH2","SPARC",
         "TAGLN","COL1A1","LGALS1","S100A4","TPM1","KRT8","KRT19","ITGB1","CALD1"],
 "DTP_persister": ["NGFR","AXL","CDKN1A","CDKN2A","SOX2","ALDH1A1","NT5E","ABCB1",
                   "DDIT4","BEX2","FTH1","FTL","SLC7A11","GPX4"],
 "quiescence_lowcycle": ["CDKN1B","CDKN1A","FOXO3","BTG1","BTG2"],
 "stress_AP1": ["FOS","FOSB","JUN","JUNB","EGR1","ATF3","DDIT3","GADD45A","HSPA1A"],
 "proliferation": ["MKI67","TOP2A","PCNA","CCNB1","CDK1","BIRC5"],
}

def maxh1(X):
    d = ripser(X[:, :30], maxdim=1)["dgms"][1]
    return float((d[:,1]-d[:,0]).max()) if len(d) else 0.0

def main():
    O = ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad")
    D = O[O.obs.timepoint == "D14"].copy()
    X = np.asarray(D.obsm["X_pca"][:, :30], float)
    # neighbors + leiden on the same PCA used for topology
    D.obsm["X_pca"] = X
    sc.pp.neighbors(D, n_pcs=30, random_state=42)
    sc.tl.leiden(D, resolution=1.0, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    clusters = D.obs["leiden"].to_numpy()
    full = maxh1(X)
    log.info("D14 full max-H1=%.2f, %d Leiden clusters", full, len(set(clusters)))

    # leave-one-cluster-out
    rows = []
    for c in sorted(set(clusters), key=int):
        keep = clusters != c
        h = maxh1(X[keep])
        rows.append((c, int((clusters==c).sum()), h, full - h))
        log.info("  drop cluster %2s (n=%4d): max-H1 %.2f  (drop %.2f)", c, (clusters==c).sum(), h, full-h)
    rows.sort(key=lambda r: -r[3])
    driver = rows[0][0]
    log.info("structure-driving cluster = %s (removal drops max-H1 by %.2f)", driver, rows[0][3])
    D.obs["driver"] = np.where(clusters == driver, "structure", "rest")

    # recover gene symbols from raw CSV (genes=rows, cells=cols) for D14 cells
    bc = list(D.obs_names)
    log.info("reading raw counts for %d D14 cells from CSV ...", len(bc))
    df = pd.read_csv(CSV, index_col=0, usecols=[0] + list(range(1, 56420)))  # all cells
    df = df[bc]                      # genes x D14cells, real symbols as index
    counts = df.T                    # cells x genes
    R = ad.AnnData(counts.values.astype(np.float32),
                   obs=D.obs.loc[counts.index, ["driver"]].copy(),
                   var=pd.DataFrame(index=counts.columns))
    R.var_names_make_unique()
    sc.pp.normalize_total(R, target_sum=1e4); sc.pp.log1p(R)
    R.obs["total_counts"] = counts.values.sum(1)

    # QC sanity: is the structure cluster just low-quality cells?
    g = R.obs.groupby("driver")["total_counts"].median()
    log.info("median total counts: structure=%.0f  rest=%.0f", g.get("structure",0), g.get("rest",0))

    # DE structure vs rest
    sc.tl.rank_genes_groups(R, "driver", groups=["structure"], reference="rest", method="wilcoxon")
    res = sc.get.rank_genes_groups_df(R, group="structure")
    res = res.sort_values("scores", ascending=False)
    top = res.head(30)
    top.to_csv(OUT / "d14_structure_DE_top.csv", index=False)
    log.info("top up in structure-driving cells:\n%s", top[["names","logfoldchanges","pvals_adj"]].head(20).to_string(index=False))

    # marker-program enrichment among the up genes
    upset = set(res.head(200)["names"])
    print("\n==== program markers UP in structure-driving cells (top-200 DE) ====")
    for prog, gs in MARKERS.items():
        present = [g for g in gs if g in R.var_names]
        hit = [g for g in present if g in upset]
        print(f"  {prog:22s} {len(hit)}/{len(present)} present-up: {hit}")

    pd.DataFrame(rows, columns=["cluster","n","maxH1_drop_to","drop"]).to_csv(OUT/"d14_leave_one_cluster_out.csv", index=False)
    log.info("wrote DE + cluster tables to %s", OUT)

if __name__ == "__main__":
    main()
