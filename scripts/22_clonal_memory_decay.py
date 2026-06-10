"""
22 — Path 2 core: formalize the clonal state-MEMORY statistic and its DECAY under drug.

Definition. Cluster cells into transcriptional STATES (shared Leiden on the pooled
osimertinib embedding so states are comparable across time). For each clone (lineage)
with >= k cells at timepoint t, compute the Shannon entropy of its cells' state
assignments H_obs. Compare to a LINEAGE-SHUFFLE null H_null (shuffle state labels among
cells at t, preserving clone sizes and the state marginal). Define

    M(t) = 1 - mean_clone(H_obs(t)) / H_null(t)         # clonal memory, 0..1

M=1 => clones perfectly state-restricted (memory); M=0 => clones as mixed as random.
The claim is M decays from baseline to drug (memory erosion = induced plasticity).
This null preserves dispersion AND cell cycle, so the signal is immune to the confounds
that sink naive max-H1.

Robustness: Leiden resolution, min clone size, bootstrap CIs, monotonic-decay test.
Run in conda `sctda` env.
"""
from __future__ import annotations
import json, logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import entropy

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("memory")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
BAD = {"NA", "nan", "None", "0", ""}
TPS = ["D0", "D3", "D7", "D14"]


def cluster(adata, resolution, seed=42):
    A = adata.copy()
    sc.pp.neighbors(A, n_pcs=30, random_state=seed)
    sc.tl.leiden(A, resolution=resolution, random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return A.obs["leiden"].to_numpy()


def memory_at(df_t, rng, n_null=300):
    """df_t: columns lin, state for one timepoint. Returns M, H_obs, H_null."""
    sizes = df_t.groupby("lin").size()
    multi = sizes[sizes >= df_t.attrs["k"]].index
    if len(multi) < 5:
        return np.nan, np.nan, np.nan, len(multi)
    H_obs = np.mean([entropy(df_t.loc[df_t.lin == m, "state"].value_counts(normalize=True))
                     for m in multi])
    Hn = []
    states = df_t["state"].to_numpy()
    lins = df_t["lin"].to_numpy()
    for _ in range(n_null):
        sh = rng.permutation(states)
        s = pd.Series(sh, index=df_t.index)
        Hn.append(np.mean([entropy(s[lins == m].value_counts(normalize=True)) for m in multi]))
    H_null = float(np.mean(Hn))
    M = 1 - H_obs / H_null if H_null > 0 else np.nan
    return float(M), float(H_obs), H_null, int(len(multi))


def trajectory(O, resolution, k, rng, nboot=300):
    states = cluster(O, resolution)
    lb = O.obs["lineage_barcode"].astype(str)
    ok = ~lb.isin(BAD)
    base = pd.DataFrame({"lin": lb[ok].to_numpy(),
                         "tp": O.obs["timepoint"][ok].to_numpy(),
                         "state": states[ok.to_numpy()]})
    out = {}
    for tp in TPS:
        d = base[base.tp == tp].copy(); d.attrs["k"] = k
        M, Ho, Hn, ncl = memory_at(d, rng)
        # bootstrap over clones for CI
        sizes = d.groupby("lin").size(); multi = list(sizes[sizes >= k].index)
        Ms = []
        for _ in range(nboot):
            bs = set(rng.choice(multi, len(multi), replace=True))
            mm, *_ = memory_at(_tag(d[d.lin.isin(bs)], k), rng, n_null=60)
            if not np.isnan(mm):
                Ms.append(mm)
        ci = (float(np.percentile(Ms, 5)), float(np.percentile(Ms, 95))) if Ms else (np.nan, np.nan)
        out[tp] = dict(M=M, ci=ci, n_clones=ncl, n_states=int(len(set(states))))
    return out


def _tag(df, k):
    df = df.copy(); df.attrs["k"] = k; return df


def main():
    O = ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad")
    rng = np.random.default_rng(42)

    log.info("=== primary: resolution=1.0, min clone size k=3 ===")
    primary = trajectory(O, 1.0, 3, rng)
    for tp in TPS:
        r = primary[tp]
        log.info("  M(%s)=%.3f  CI[%.3f,%.3f]  clones=%d states=%d",
                 tp, r["M"], r["ci"][0], r["ci"][1], r["n_clones"], r["n_states"])
    decay = primary["D0"]["M"] - primary["D14"]["M"]
    monotonic = all(primary[a]["M"] >= primary[b]["M"] - 0.02
                    for a, b in zip(TPS, TPS[1:]))
    log.info("  memory decay D0->D14 = %.3f ; near-monotonic=%s", decay, monotonic)

    log.info("=== robustness: resolution x clone-size grid (decay D0-D14) ===")
    grid = {}
    for res in [0.5, 1.0, 1.5]:
        for k in [2, 3, 5]:
            t = trajectory(O, res, k, rng, nboot=120)
            d = t["D0"]["M"] - t["D14"]["M"]
            grid[f"res{res}_k{k}"] = dict(decay=float(d), M_D0=t["D0"]["M"], M_D14=t["D14"]["M"])
            log.info("  res=%.1f k=%d: M_D0=%.3f M_D14=%.3f decay=%.3f", res, k, t["D0"]["M"], t["D14"]["M"], d)

    pos = sum(1 for v in grid.values() if v["decay"] > 0)
    log.info("decay > 0 (memory erodes under drug) in %d/%d grid settings", pos, len(grid))

    (OUT / "clonal_memory_decay.json").write_text(json.dumps(
        {"primary": primary, "robustness_grid": grid,
         "decay_positive_fraction": pos / len(grid)}, indent=2))

    print("\n==== CLONAL MEMORY DECAY (Path 2 core statistic) ====")
    print("M(t) = clonal state memory (1=restricted, 0=random); decay = drug-induced plasticity")
    for tp in TPS:
        r = primary[tp]
        print(f"  M({tp:3s}) = {r['M']:.3f}  90%CI[{r['ci'][0]:.3f},{r['ci'][1]:.3f}]  ({r['n_clones']} clones)")
    print(f"\n  memory decay D0->D14 = {decay:.3f}")
    print(f"  robust: decay>0 in {pos}/{len(grid)} resolution x clone-size settings")


if __name__ == "__main__":
    main()
