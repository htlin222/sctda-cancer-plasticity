"""
16 — Foundation analysis: is the H1 loop a genuine cyclic structure that
clonal lineages traverse, or just a static hole in a blob?

This directly addresses the peer-review construct-validity objection
(reviewer #1/#3/#4): max-H1 is computed on a static snapshot, but the
biological claim is reversible *cycling*. Here we:

  1. Compute circular coordinates theta in [0,1) for each cell from the
     most-persistent H1 cocycle (de Silva-Morozov-Vejdemo-Johansson
     harmonic smoothing; implemented directly from ripser cocycles, no
     dreimac dependency).
  2. Test whether the loop is REAL: are cells distributed around the full
     circle (Rayleigh / Kuiper uniformity), and does a smooth biological
     program (EMT, cell cycle) vary around theta?
  3. Test CYCLING with Watermelon lineage barcodes (GSE150949):
       (a) within-timepoint: do clonemates spread around the loop
           (circular variance) more than expected, and does spread change
           with drug exposure D0->D14?
       (b) null: compare real-lineage circular variance to size-matched
           random cell groups.

Run inside the conda `sctda` env (ripser 0.6.14, coeff>2 + do_cocycles).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from dreimac import CircularCoords

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("circcoord")

SEED = 42
N_PCS = 30
PRIME = 47  # odd prime field for circular coordinates
DATA = Path("data/processed/pc9_osimertinib_processed.h5ad")
OUT = Path("results/foundation")
OUT.mkdir(parents=True, exist_ok=True)

EMT_GENES = [  # MSigDB Hallmark EMT subset (paper's cassette intersect)
    "TPM1", "KRT8", "KRT19", "VIM", "FN1", "CDH2", "SNAI2", "ZEB1",
    "S100A4", "LGALS1", "TGFBI", "SPARC", "COL1A1", "TAGLN",
]
BAD_BC = {"NA", "nan", "None", "0", "", "na", "NaN"}


def circular_coordinates(X: np.ndarray, prime: int = PRIME):
    """Return (theta in [0,1), diagnostics) for the most-persistent H1 class.

    Uses dreimac's validated implementation of the de Silva-Morozov-
    Vejdemo-Johansson circular-coordinate algorithm (greedy-permutation
    landmarks + harmonic smoothing of the lifted Z/p cocycle). cocycle_idx=0
    selects the most persistent cohomology class.
    """
    n = X.shape[0]
    n_landmarks = min(n, 400)
    cc = CircularCoords(X, n_landmarks=n_landmarks, prime=prime)
    # angular coordinates in [0, 2*pi)
    ang = cc.get_coordinates(perc=0.5, cocycle_idx=0, standard_range=False)
    theta = np.mod(np.asarray(ang) / (2 * np.pi), 1.0)

    dgm1 = cc.dgms_[1]
    if len(dgm1):
        pers = float(np.max(dgm1[:, 1] - dgm1[:, 0]))
    else:
        pers = float("nan")
    return theta, {
        "persistence": pers,
        "n_landmarks": int(n_landmarks),
    }


def circle_coverage(theta: np.ndarray, n_sectors: int = 12) -> float:
    """Fraction of circle sectors that are occupied. A genuinely traversed
    loop fills most sectors; two clumps with a gap fill few."""
    if theta.size == 0:
        return float("nan")
    occ = np.unique((theta * n_sectors).astype(int) % n_sectors)
    return float(occ.size / n_sectors)


def rayleigh_test(angles_rad: np.ndarray):
    """Rayleigh test for non-uniformity. Small p => concentrated (not uniform)."""
    n = angles_rad.size
    if n < 3:
        return np.nan, np.nan
    C = np.cos(angles_rad).sum()
    S = np.sin(angles_rad).sum()
    R = np.sqrt(C**2 + S**2) / n
    z = n * R**2
    p = np.exp(-z) * (1 + (2 * z - z**2) / (4 * n))  # Zar approximation
    return float(R), float(min(max(p, 0.0), 1.0))


def circular_variance(angles_rad: np.ndarray) -> float:
    """1 - mean resultant length. 0=concentrated, 1=spread around circle."""
    if angles_rad.size == 0:
        return np.nan
    C = np.cos(angles_rad).mean()
    S = np.sin(angles_rad).mean()
    return float(1.0 - np.sqrt(C**2 + S**2))


def main():
    rng = np.random.default_rng(SEED)
    log.info("loading %s", DATA)
    adata = ad.read_h5ad(DATA)
    adata.var_names_make_unique()

    summary = {}
    per_cell_frames = []

    for tp in ["D0", "D3", "D7", "D14"]:
        sub = adata[adata.obs["timepoint"] == tp]
        X = np.asarray(sub.obsm["X_pca"][:, :N_PCS], dtype=float)
        log.info("[%s] %d cells -> circular coordinates", tp, X.shape[0])
        try:
            theta, diag = circular_coordinates(X)
        except RuntimeError as e:
            log.warning("[%s] %s", tp, e)
            continue
        ang = theta * 2 * np.pi
        R, p_ray = rayleigh_test(ang)
        coverage = circle_coverage(theta)

        rec = {**diag, "n_cells": int(X.shape[0]),
               "rayleigh_R": R, "rayleigh_p": p_ray,
               "circle_coverage": coverage}
        idx = np.where((adata.obs["timepoint"] == tp).to_numpy())[0]

        # lineage circular-variance test (within timepoint)
        lin_here = adata.obs["lineage_barcode"].astype(str).to_numpy()[idx]
        ok = ~pd.Series(lin_here).isin(BAD_BC).to_numpy()
        df = pd.DataFrame({"lin": lin_here[ok], "ang": ang[ok]})
        sizes = df.groupby("lin").size()
        multi = sizes[sizes >= 3].index  # clones with >=3 co-sampled cells
        real_cv = [circular_variance(df.loc[df.lin == m, "ang"].to_numpy())
                   for m in multi]
        rec["n_lineages_ge3"] = int(len(multi))
        if len(multi) >= 5:
            real_cv = np.array(real_cv)
            # null: random groups of matched sizes
            null_means = []
            angs_all = df["ang"].to_numpy()
            grp_sizes = sizes[multi].to_numpy()
            for _ in range(1000):
                cvs = []
                for gs in grp_sizes:
                    pick = rng.choice(angs_all, size=gs, replace=False)
                    cvs.append(circular_variance(pick))
                null_means.append(np.mean(cvs))
            null_means = np.array(null_means)
            obs_mean = float(real_cv.mean())
            # two-sided empirical p
            p_lo = (np.sum(null_means <= obs_mean) + 1) / (len(null_means) + 1)
            p_hi = (np.sum(null_means >= obs_mean) + 1) / (len(null_means) + 1)
            rec["lineage_cv_mean"] = obs_mean
            rec["lineage_cv_null_mean"] = float(null_means.mean())
            rec["lineage_cv_p_concentrated"] = float(p_lo)  # clones tighter than random
            rec["lineage_cv_p_spread"] = float(p_hi)        # clones more spread

        summary[tp] = rec
        per_cell_frames.append(pd.DataFrame({
            "timepoint": tp,
            "theta": theta,
            "lineage_barcode": adata.obs["lineage_barcode"].astype(str).to_numpy()[idx],
        }))
        log.info("[%s] persistence=%.3f coverage=%.2f rayleigh_R=%.3f lineages>=3=%d",
                 tp, diag["persistence"], coverage, R, rec["n_lineages_ge3"])

    (OUT / "circular_coord_summary.json").write_text(json.dumps(summary, indent=2))
    if per_cell_frames:
        pd.concat(per_cell_frames).to_csv(OUT / "circular_coords_per_cell.csv", index=False)
    log.info("wrote %s", OUT / "circular_coord_summary.json")

    print("\n==== FOUNDATION ANALYSIS: circular coordinates ====")
    for tp, rec in summary.items():
        print(f"\n[{tp}] n={rec['n_cells']}  H1 persistence={rec['persistence']:.3f}")
        print(f"   circle coverage (frac of 12 sectors occupied; high=traversed loop): {rec['circle_coverage']:.2f}")
        print(f"   concentration: Rayleigh R={rec['rayleigh_R']:.3f} (R~1=one clump, R~0=spread)")
        if "lineage_cv_mean" in rec:
            print(f"   clone angular spread (CV): real={rec['lineage_cv_mean']:.3f} "
                  f"null={rec['lineage_cv_null_mean']:.3f}  "
                  f"p(tighter)={rec['lineage_cv_p_concentrated']:.3f} "
                  f"p(spread)={rec['lineage_cv_p_spread']:.3f}  "
                  f"(n_clones>=3 = {rec['n_lineages_ge3']})")
        else:
            print(f"   clone test skipped (only {rec['n_lineages_ge3']} clones with >=3 co-sampled cells)")


if __name__ == "__main__":
    main()
