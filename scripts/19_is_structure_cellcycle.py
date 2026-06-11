"""
19 — Make-or-break: is the drug-induced non-Gaussian structure just the
CELL-CYCLE loop?

Script 18 showed the cells driving max-H1 at D14 are proliferating (cell-cycle)
cells. Cell cycle is intrinsically a loop, so the "structure beyond dispersion"
(real > Gaussian null) could be the cell-cycle loop, not drug-tolerant biology.

Test: rebuild D0/D7/D14 from raw counts WITH recovered gene symbols, then for
each compute real max-H1 and the covariance-matched Gaussian-null max-H1 under
two preprocessings:
  (i)  standard (paper pipeline; cell cycle NOT regressed)
  (ii) S/G2M score regressed out before PCA (proper cell-cycle removal)

If the excess (real > gauss) and the drug increase (D14 > D0) SURVIVE (ii),
the structure is not merely cell cycle -> rescue holds.
If they vanish under (ii) -> the signal is the cell-cycle loop -> rescue fails.

Run in conda `sctda` env.
"""
from __future__ import annotations
import logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from ripser import ripser

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("cc")
CSV = "data/raw/GSE150949/GSE150949_pc9_count_matrix.csv.gz"
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)

# Tirosh S and G2M gene lists (standard scanpy set)
S_GENES = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M_GENES = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()

def maxh1(X):
    d = ripser(X[:, :30], maxdim=1)["dgms"][1]
    return float((d[:,1]-d[:,0]).max()) if len(d) else 0.0

def gauss_null(Xs, ngauss, rng):
    mu = Xs.mean(0); cov = np.cov(Xs, rowvar=False)
    L = np.linalg.cholesky(cov + 1e-9*np.eye(cov.shape[0]))
    return np.median([maxh1(mu + rng.standard_normal(Xs.shape) @ L.T) for _ in range(ngauss)])

def build(counts, regress_cc):
    A = ad.AnnData(counts.values.astype(np.float32),
                   var=pd.DataFrame(index=counts.columns))
    A.var_names_make_unique()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.tl.score_genes_cell_cycle(A, s_genes=[g for g in S_GENES if g in A.var_names],
                                 g2m_genes=[g for g in G2M_GENES if g in A.var_names])
    sc.pp.highly_variable_genes(A, n_top_genes=3000); A = A[:, A.var.highly_variable].copy()
    if regress_cc:
        sc.pp.regress_out(A, ["S_score", "G2M_score"])
    sc.pp.scale(A, max_value=10)
    sc.tl.pca(A, n_comps=50, random_state=42)
    return A

def evaluate(name, counts, rng, ngauss=15):
    out = {}
    for tag, reg in [("std", False), ("ccreg", True)]:
        A = build(counts, reg)
        X = np.asarray(A.obsm["X_pca"][:, :30], float)
        real = maxh1(X); g = gauss_null(X, ngauss, rng)
        out[tag] = (real, g, real/g)
        log.info("%-6s %-5s  realH1=%.2f gauss=%.2f ratio=%.2f", name, tag, real, g, real/g)
    return out

def main():
    O = ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad")
    want = {tp: list(O.obs_names[O.obs.timepoint == tp]) for tp in ["D0", "D7", "D14"]}
    allbc = [b for v in want.values() for b in v]
    log.info("reading raw counts for %d cells (D0/D7/D14) ...", len(allbc))
    df = pd.read_csv(CSV, index_col=0, usecols=[0] + list(range(1, 56420)))
    rng = np.random.default_rng(42)
    res = {}
    for tp in ["D0", "D7", "D14"]:
        counts = df[want[tp]].T   # cells x genes
        res[tp] = evaluate(tp, counts, rng)

    print("\n==== IS THE STRUCTURE CELL CYCLE? real/gauss ratio ====")
    print(f"{'tp':5s} {'std ratio':>10s} {'ccreg ratio':>12s}  verdict")
    for tp in ["D0", "D7", "D14"]:
        s = res[tp]["std"][2]; c = res[tp]["ccreg"][2]
        print(f"{tp:5s} {s:10.2f} {c:12.2f}  {'excess survives CC removal' if c>1.15 else 'excess GONE after CC removal'}")
    d0c, d14c = res["D0"]["ccreg"][0], res["D14"]["ccreg"][0]
    print(f"\nAfter CC regression: D14 max-H1={d14c:.2f} vs D0={d0c:.2f}  "
          f"({'drug increase survives' if d14c>d0c else 'drug increase GONE'})")
    pd.DataFrame({tp: {f"{t}_{k}": v for t in ["std","ccreg"] for k,v in
                       zip(["real","gauss","ratio"], res[tp][t])} for tp in res}
                 ).to_csv(OUT/"cellcycle_structure_test.csv")

if __name__ == "__main__":
    main()
