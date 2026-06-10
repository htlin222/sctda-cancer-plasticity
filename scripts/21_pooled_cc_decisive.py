"""
21 — Decisive test: in the POOLED embedding (the paper's method, which is where
the excess-over-dispersion lives), does the drug-induced excess survive proper
cell-cycle regression? Bootstrapped.

Build the osimertinib pooled embedding from raw counts (recovered gene symbols)
two ways -- standard and S/G2M-regressed -- then bootstrap real max-H1 and the
covariance-matched Gaussian-null ratio for D0 and D14.

Verdict logic:
  - std D14 ratio should reproduce the stored-embedding excess (>~1.3).
  - if ccreg D14 ratio stays >~1.2 and > D0 -> excess is NOT just cell cycle (rescue holds).
  - if ccreg D14 ratio ~1.0 -> excess IS the cell-cycle loop (rescue fails).

Run in conda `sctda` env.
"""
from __future__ import annotations
import logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from ripser import ripser

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("pooled")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
CSV = "data/raw/GSE150949/GSE150949_pc9_count_matrix.csv.gz"
S="MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G="HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE NEK2 CBX5 CENPA".split()

def maxh1(X):
    d=ripser(X[:, :30],maxdim=1)["dgms"][1]; return float((d[:,1]-d[:,0]).max()) if len(d) else 0.0
def gnull(Xs,ng,rng):
    mu=Xs.mean(0);cov=np.cov(Xs,rowvar=False);L=np.linalg.cholesky(cov+1e-9*np.eye(cov.shape[0]))
    return np.median([maxh1(mu+rng.standard_normal(Xs.shape)@L.T) for _ in range(ng)])

def pooled_embed(counts, tp_labels, regress):
    A=ad.AnnData(counts.values.astype(np.float32),var=pd.DataFrame(index=counts.columns))
    A.var_names_make_unique(); A.obs["timepoint"]=list(tp_labels)
    sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    sc.tl.score_genes_cell_cycle(A,s_genes=[g for g in S if g in A.var_names],g2m_genes=[g for g in G if g in A.var_names])
    sc.pp.highly_variable_genes(A,n_top_genes=3000,batch_key="timepoint")
    A=A[:,A.var.highly_variable].copy()
    if regress: sc.pp.regress_out(A,["S_score","G2M_score"])
    sc.pp.scale(A,max_value=10); sc.tl.pca(A,n_comps=50,random_state=42)
    return A

def boot_ratio(X, nboot=15, ngauss=8, n=1000, seed=0):
    rng=np.random.default_rng(seed); nn=min(n,X.shape[0]); rs=[]
    for _ in range(nboot):
        Xs=X[rng.choice(X.shape[0],nn,replace=False)]
        rs.append(maxh1(Xs)/gnull(Xs,ngauss,rng))
    return np.array(rs)

def main():
    O=ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad")
    bc=list(O.obs_names); tp=O.obs["timepoint"].to_numpy()
    log.info("reading raw counts for %d osimertinib cells ...",len(bc))
    df=pd.read_csv(CSV,index_col=0,usecols=[0]+list(range(1,56420)))
    counts=df[bc].T
    for regress,tag in [(False,"std"),(True,"ccreg")]:
        A=pooled_embed(counts,tp,regress)
        res={}
        for t in ["D0","D14"]:
            X=np.asarray(A[A.obs.timepoint==t].obsm["X_pca"][:, :30],float)
            r=boot_ratio(X)
            res[t]=r
            log.info("%-5s %-5s ratio median=%.2f [%.2f,%.2f]  (realH1 med=%.2f)",
                     tag,t,np.median(r),np.percentile(r,10),np.percentile(r,90),maxh1(X))
        print(f"\n[{tag}] D14 ratio={np.median(res['D14']):.2f}  D0 ratio={np.median(res['D0']):.2f}  "
              f"D14>D0 in {100*np.mean(res['D14'][:,None]>res['D0'][None,:]):.0f}% ; "
              f"D14 excess(>1) in {100*np.mean(res['D14']>1):.0f}% of boots")

if __name__ == "__main__":
    main()
