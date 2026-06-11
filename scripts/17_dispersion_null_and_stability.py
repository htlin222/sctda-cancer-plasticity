"""
17 — The rescued core: is drug-induced max-H1 genuine topological structure,
or just increased transcriptional dispersion?

Two analyses that directly answer the peer-review killer objections:

  A. COVARIANCE-MATCHED GAUSSIAN NULL (reviewer #2's demanded control).
     For each cohort group, compare real max-H1 to max-H1 of a multivariate
     Gaussian with the SAME mean and covariance (identical dispersion, no
     real structure). ratio = real / gauss. ratio ~ 1 => statistic is a pure
     dispersion proxy; ratio >> 1 (real > null) => genuine non-Gaussian
     topological structure beyond spread.

  B. BOOTSTRAP STABILITY of the drug-vs-baseline gap (reviewer #3).
     Point estimates of max-H1 are subsample-unstable (CV 7-19%), so we test
     whether the ORDERING (drug > baseline; real > dispersion-null) is robust
     across many independent subsamples.

Conclusion from the run: under osimertinib and in EGFR-mutant PDX, max-H1
exceeds the dispersion-matched null under drug (not at baseline) and
drug > baseline in ~100% of bootstrap pairs. Erlotinib shows NO excess
structure (it dies). This is the evidence base for the reframed manuscript.

Run in conda `sctda` env.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
from ripser import ripser

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("dispnull")

N_PCS = 30
OUT = Path("results/foundation")
OUT.mkdir(parents=True, exist_ok=True)


def maxh1(X: np.ndarray) -> float:
    d = ripser(X, maxdim=1)["dgms"][1]
    return float((d[:, 1] - d[:, 0]).max()) if len(d) else 0.0


def gaussian_matched_h1(Xs: np.ndarray, ngauss: int, rng) -> np.ndarray:
    mu = Xs.mean(0)
    cov = np.cov(Xs, rowvar=False)
    L = np.linalg.cholesky(cov + 1e-9 * np.eye(cov.shape[0]))
    return np.array([maxh1(mu + rng.standard_normal(Xs.shape) @ L.T)
                     for _ in range(ngauss)])


def analyse(name, Xraw, nboot=30, ngauss=15, n=1200):
    X = np.asarray(Xraw, float)[:, :N_PCS]
    rng = np.random.default_rng(42)
    nn = min(n, X.shape[0])
    real, gauss = [], []
    for _ in range(nboot):
        idx = rng.choice(X.shape[0], nn, replace=False)
        Xs = X[idx]
        real.append(maxh1(Xs))
        gauss.extend(gaussian_matched_h1(Xs, max(1, ngauss // nboot * 1), rng).tolist())
    # one richer gaussian draw on a representative subsample
    idx = rng.choice(X.shape[0], nn, replace=False)
    gauss = np.array(gauss + gaussian_matched_h1(X[idx], ngauss, rng).tolist())
    real = np.array(real)
    rec = dict(
        n_cells=int(nn),
        real_median=float(np.median(real)),
        real_ci90=[float(np.percentile(real, 5)), float(np.percentile(real, 95))],
        real_cv_pct=float(100 * real.std() / real.mean()),
        gauss_median=float(np.median(gauss)),
        ratio=float(np.median(real) / np.median(gauss)),
        pct_real_gt_gauss=float(100 * np.mean(real[:, None] > gauss[None, :])),
    )
    log.info("%-22s real=%.2f[%.2f,%.2f] CV=%.0f%% gauss=%.2f ratio=%.2f real>gauss=%.0f%%",
             name, rec["real_median"], rec["real_ci90"][0], rec["real_ci90"][1],
             rec["real_cv_pct"], rec["gauss_median"], rec["ratio"], rec["pct_real_gt_gauss"])
    return rec, real


GROUPS = []  # (label, file, obs_filter)
def load():
    O = ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad"); O.var_names_make_unique()
    E = ad.read_h5ad("data/processed/pc9_erlotinib_processed.h5ad"); E.var_names_make_unique()
    P = ad.read_h5ad("data/processed/pdx_ascl1_osimertinib_processed.h5ad"); P.var_names_make_unique()
    M = ad.read_h5ad("data/processed/maynard2020_egfr_full_processed.h5ad"); M.var_names_make_unique()
    g = {}
    for tp in ["D0", "D3", "D7", "D14"]:
        g[f"osi {tp}"] = O[O.obs.timepoint == tp].obsm["X_pca"]
    for tp in ["D0", "D9", "D11"]:
        g[f"erl {tp}"] = E[E.obs.timepoint == tp].obsm["X_pca"]
    for c in ["untreated", "residual"]:
        g[f"PDX-YU006 {c}"] = P[(P.obs.pdx == "YU-006") & (P.obs.condition == c)].obsm["X_pca"]
    for tp in ["TN", "PER", "PD"]:
        g[f"Maynard {tp}"] = M[M.obs.timepoint == tp].obsm["X_pca"]
    return g


def main():
    groups = load()
    summary, reals = {}, {}
    for name, X in groups.items():
        rec, real = analyse(name, X)
        summary[name] = rec
        reals[name] = real

    # directional robustness for the key contrasts
    def pairgt(a, b):
        return float(100 * np.mean(reals[a][:, None] > reals[b][None, :]))
    contrasts = {
        "osi D14 > D0": pairgt("osi D14", "osi D0"),
        "osi D7 > D0": pairgt("osi D7", "osi D0"),
        "PDX residual > untreated": pairgt("PDX-YU006 residual", "PDX-YU006 untreated"),
        "erl D9 > D0": pairgt("erl D9", "erl D0"),
    }
    summary["_contrasts_pct_pairs"] = contrasts
    (OUT / "dispersion_null_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n==== RESCUED CORE: structure beyond dispersion ====")
    print(f"{'group':22s} {'realH1':>16s} {'gauss':>6s} {'ratio':>6s} {'real>gauss':>10s}")
    for name, r in summary.items():
        if name.startswith("_"):
            continue
        ci = r["real_ci90"]
        print(f"{name:22s} {r['real_median']:5.2f}[{ci[0]:4.2f},{ci[1]:4.2f}] "
              f"{r['gauss_median']:6.2f} {r['ratio']:6.2f} {r['pct_real_gt_gauss']:9.0f}%")
    print("\nDirectional robustness (% of bootstrap pairs):")
    for k, v in contrasts.items():
        print(f"   {k:28s} {v:5.0f}%")
    log.info("wrote %s", OUT / "dispersion_null_summary.json")


if __name__ == "__main__":
    main()
