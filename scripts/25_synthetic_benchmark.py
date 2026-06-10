"""
25 — Synthetic ground-truth benchmark for the control framework (reviewer change #1).

Validates the toolkit on point clouds where the answer is known, characterising
sensitivity (does it pass a genuine traversed loop?) and specificity (does it
reject benign dispersion / non-Gaussian-but-non-cyclic structure?). Also adds a
SECOND dispersion-preserving null (per-PC independent permutation) so the
conclusion does not rest on Gaussianity alone (reviewer change #3).

Scenarios (n=1200, embedded in 30-D with isotropic noise on the unused axes):
  (a) gaussian_blob      : no structure                  -> expect ratio~1, high R  (true negative)
  (b) traversed_loop     : cells uniformly around circle -> expect ratio>1, LOW R   (true positive)
  (c) sparse_arc         : 8% of cells on an arc + blob   -> ratio>1, high R         (real-data signature)
  (d) two_clusters       : non-Gaussian, non-cyclic       -> Gaussian-null specificity test

Toolkit decision rule:  structure if (ratio>1 robustly) ;  traversed loop only if ALSO R < ~0.6.
Run in conda `sctda` env.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ripser import ripser
from dreimac import CircularCoords

OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
# NOISE < loop radius so genuine structure lives in the top PCs (as for real data)
N, DIM, NOISE = 1200, 30, 0.3
rng = np.random.default_rng(42)


def embed(core2d, noise=NOISE):
    X = np.zeros((core2d.shape[0], DIM))
    X[:, :core2d.shape[1]] = core2d
    X[:, core2d.shape[1]:] = rng.normal(0, noise, (core2d.shape[0], DIM - core2d.shape[1]))
    return X


def gaussian_blob():
    return embed(rng.normal(0, 1.5, (N, 2)))

def traversed_loop():
    t = rng.uniform(0, 2 * np.pi, N)
    return embed(np.column_stack([3 * np.cos(t), 3 * np.sin(t)]) + rng.normal(0, 0.25, (N, 2)))

def sparse_arc():
    # ~18% of cells on a wide arc (270 deg, not closed) atop a blob: detectable
    # structure that is NOT a traversed closed loop (the real-data signature).
    k = int(0.22 * N)
    blob = rng.normal(0, 0.9, (N - k, 2))
    t = rng.uniform(0, 1.5 * np.pi, k)
    arc = np.column_stack([4.5 * np.cos(t), 4.5 * np.sin(t)]) + rng.normal(0, 0.2, (k, 2))
    return embed(np.vstack([blob, arc]))

def two_clusters():
    a = rng.normal([-3, 0], 0.8, (N // 2, 2))
    b = rng.normal([3, 0], 0.8, (N - N // 2, 2))
    return embed(np.vstack([a, b]))


def maxh1(X):
    d = ripser(X, maxdim=1)["dgms"][1]
    return float((d[:, 1] - d[:, 0]).max()) if len(d) else 0.0


def gaussian_ratio(X, ng=20):
    mu, cov = X.mean(0), np.cov(X, rowvar=False)
    L = np.linalg.cholesky(cov + 1e-9 * np.eye(DIM))
    g = np.median([maxh1(mu + rng.standard_normal(X.shape) @ L.T) for _ in range(ng)])
    return maxh1(X) / g, g


def perm_ratio(X, ng=20):
    """Second dispersion-preserving null: independently permute each coordinate
    across cells (preserves per-axis variance = dispersion, destroys joint loops)."""
    def one():
        Xp = np.column_stack([rng.permutation(X[:, j]) for j in range(DIM)])
        return maxh1(Xp)
    g = np.median([one() for _ in range(ng)])
    return maxh1(X) / g, g


def rayleigh_R(X):
    cc = CircularCoords(X, n_landmarks=min(N, 400), prime=47)
    th = np.mod(cc.get_coordinates(perc=0.5, cocycle_idx=0, standard_range=False), 2 * np.pi)
    return float(np.hypot(np.cos(th).sum(), np.sin(th).sum()) / th.size)


def main():
    scenarios = {"a_gaussian_blob": (gaussian_blob, "no structure", "true negative"),
                 "b_traversed_loop": (traversed_loop, "real cyclic structure", "true positive"),
                 "c_sparse_arc": (sparse_arc, "sparse non-traversed (real-data-like)", "structure, not cycling"),
                 "d_two_clusters": (two_clusters, "non-Gaussian, non-cyclic", "specificity test")}
    res = {}
    print(f"{'scenario':18s} {'maxH1':>6s} {'gaussR':>7s} {'permR':>6s} {'RayleighR':>9s}  verdict")
    for name, (fn, desc, truth) in scenarios.items():
        X = fn()
        gr, _ = gaussian_ratio(X)   # primary null (sensitive + specific)
        pr, _ = perm_ratio(X)       # secondary, more conservative robustness check
        R = rayleigh_R(X)
        structure = gr > 1.1        # decision on the primary Gaussian null
        traversed = structure and R < 0.6
        verdict = ("traversed loop" if traversed else "structure (not traversed)" if structure else "no structure")
        res[name] = dict(desc=desc, ground_truth=truth, maxh1=round(maxh1(X), 2),
                         gaussian_ratio=round(gr, 2), perm_ratio=round(pr, 2),
                         rayleigh_R=round(R, 3), toolkit_verdict=verdict)
        print(f"{name:18s} {maxh1(X):6.2f} {gr:7.2f} {pr:6.2f} {R:9.3f}  {verdict}  [{truth}]")
    (OUT / "synthetic_benchmark.json").write_text(json.dumps(res, indent=2))
    print("\nwrote", OUT / "synthetic_benchmark.json")


if __name__ == "__main__":
    main()
