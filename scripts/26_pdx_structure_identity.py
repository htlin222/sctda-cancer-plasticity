"""
26 — What IS the genuine, non-cell-cycle topological structure in the PDX residual?

The PDX YU-006 residual carries excess max-H1 over the dispersion null that survives
removing cycling cells (script: G1-only ratio 2.07, 100% bootstrap). Here we identify
the cells that create it (leave-one-cluster-out on max-H1) and characterise them by
differential expression + marker-program scoring, using the .raw gene layer.

If they are a coherent, recognisable program (EMT / mesenchymal / neuroendocrine /
alveolar / stress), the cautionary paper gains a positive finding.
Run in conda `sctda` env.
"""
from __future__ import annotations
import json, logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from ripser import ripser

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("pdxid")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)

PROGRAMS = {
 "EMT_mesenchymal": ["VIM","FN1","SERPINE1","TGFBI","ZEB1","ZEB2","SNAI2","CDH2","SPARC",
                     "TAGLN","COL1A1","COL3A1","LGALS1","S100A4","LGALS3","ACTA2","FBN1","TIMP1"],
 "epithelial": ["EPCAM","CDH1","KRT8","KRT18","KRT19","KRT7","CLDN4","CLDN7"],
 "alveolar_AT1_AT2": ["SFTPC","SFTPB","SFTPA1","NAPSA","AGER","PDPN","HOPX","LAMP3"],
 "neuroendocrine": ["ASCL1","CHGA","SYP","NCAM1","INSM1","CALCA"],
 "drug_tolerant_persister": ["NGFR","AXL","CDKN1A","CDKN2A","SOX2","ALDH1A1","NT5E","ABCB1",
                             "DDIT4","SLC7A11","GPX4","FTH1","FTL"],
 "stress_AP1_hypoxia": ["FOS","FOSB","JUN","JUNB","EGR1","ATF3","DDIT3","VEGFA","CA9","NDRG1","HILPDA"],
 "interferon_inflam": ["ISG15","IFI6","IFI27","STAT1","IRF7","CXCL10","B2M","HLA-A"],
}


def maxh1(X):
    d = ripser(X[:, :30], maxdim=1)["dgms"][1]
    return float((d[:, 1] - d[:, 0]).max()) if len(d) else 0.0


def main():
    P = ad.read_h5ad("data/processed/pdx_ascl1_osimertinib_processed.h5ad")
    D = P[(P.obs.pdx == "YU-006") & (P.obs.condition == "residual")].copy()
    rng = np.random.default_rng(42)
    if D.n_obs > 1500:
        D = D[rng.choice(D.n_obs, 1500, replace=False)].copy()
    X = np.asarray(D.obsm["X_pca"][:, :30], float)
    sc.pp.neighbors(D, n_pcs=30, random_state=42)
    sc.tl.leiden(D, resolution=1.0, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    clusters = D.obs["leiden"].to_numpy()
    full = maxh1(X)
    log.info("PDX YU-006 residual: %d cells, max-H1=%.2f, %d clusters", D.n_obs, full, len(set(clusters)))
    rows = []
    for c in sorted(set(clusters), key=int):
        h = maxh1(X[clusters != c])
        rows.append((c, int((clusters == c).sum()), full - h))
        log.info("  drop cluster %2s (n=%4d): drop max-H1 by %.2f", c, (clusters == c).sum(), full - h)
    rows.sort(key=lambda r: -r[2])
    driver = rows[0][0]
    log.info("structure-driving cluster = %s (drop %.2f)", driver, rows[0][2])
    D.obs["driver"] = np.where(clusters == driver, "structure", "rest")

    # raw expression for DE + scoring
    raw = D.raw.to_adata() if D.raw is not None else None
    if raw is None:
        full_raw = P.raw.to_adata()[D.obs_names]
        raw = full_raw
    raw.var_names = [v.replace("hg38_", "") for v in raw.var_names]
    raw.var_names_make_unique()
    sc.pp.normalize_total(raw, target_sum=1e4); sc.pp.log1p(raw)
    raw.obs["driver"] = D.obs["driver"].values

    sc.tl.rank_genes_groups(raw, "driver", groups=["structure"], reference="rest", method="wilcoxon")
    de = sc.get.rank_genes_groups_df(raw, group="structure").sort_values("scores", ascending=False)
    de.head(40).to_csv(OUT / "pdx_structure_DE_top.csv", index=False)
    log.info("top genes UP in structure cells:\n%s",
             de.head(25)[["names", "logfoldchanges", "pvals_adj"]].to_string(index=False))

    upset = set(de.head(150)["names"])
    prog_hits = {}
    print("\n==== marker-program enrichment in PDX structure-driving cells (top-150 DE) ====")
    for prog, genes in PROGRAMS.items():
        present = [g for g in genes if g in raw.var_names]
        hit = [g for g in present if g in upset]
        prog_hits[prog] = {"present": len(present), "up": hit}
        # also score-genes mean difference
        if len(present) >= 3:
            sc.tl.score_genes(raw, present, score_name=f"score_{prog}")
            ms = raw.obs.loc[raw.obs.driver == "structure", f"score_{prog}"].mean()
            mr = raw.obs.loc[raw.obs.driver == "rest", f"score_{prog}"].mean()
            prog_hits[prog]["score_struct_minus_rest"] = round(float(ms - mr), 3)
        print(f"  {prog:24s} up:{len(hit)}/{len(present)} {hit}"
              f"  Δscore={prog_hits[prog].get('score_struct_minus_rest')}")
    (OUT / "pdx_structure_identity.json").write_text(json.dumps(
        {"driver_cluster": driver, "maxh1_drop": rows[0][2],
         "top_DE": de.head(25)["names"].tolist(), "programs": prog_hits}, indent=2))
    log.info("wrote %s", OUT / "pdx_structure_identity.json")


if __name__ == "__main__":
    main()
