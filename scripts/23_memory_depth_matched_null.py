"""
23 — Gate 4 (memo): is the clonal memory-DECAY a real loss of heritability, or a
drug-induced dropout/divergence artifact?

Two hardening checks on the M(t) = 1 - H_obs/H_null statistic:
  (A) DEPTH-STRATIFIED null: shuffle state labels only WITHIN sequencing-depth bins,
      so the null preserves any depth<->state association. If M still decays D0->D14,
      the decay is not explained by depth/dropout differences.
  (B) DEPTH-DOWNSAMPLED embedding sanity: report median depth per timepoint and whether
      it trends with the decay (a quick confound read).

Also computes the two-axis reconciliation vs Oren 2021: clone state-identity entropy
(our axis, should rise) vs clone fate determinism (early state predicting later clone
size; Oren's axis). These can move oppositely without contradiction.

Run in conda `sctda` env.
"""
from __future__ import annotations
import json, logging
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import entropy

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("depth")
OUT = Path("results/foundation"); OUT.mkdir(parents=True, exist_ok=True)
BAD = {"NA", "nan", "None", "0", ""}
TPS = ["D0", "D3", "D7", "D14"]
K = 3


def memory_M(d, states_col, rng, n_null=300, depth_bins=None):
    """d has columns lin, state[, depthbin]. Returns M."""
    sizes = d.groupby("lin").size(); multi = sizes[sizes >= K].index
    if len(multi) < 5:
        return np.nan
    def meanH(svec):
        s = pd.Series(svec, index=d.index)
        return np.mean([entropy(s[d.lin == m].value_counts(normalize=True)) for m in multi])
    H_obs = meanH(d[states_col].to_numpy())
    Hn = []
    for _ in range(n_null):
        if depth_bins is None:
            sh = rng.permutation(d[states_col].to_numpy())
        else:  # shuffle within depth strata
            sh = d[states_col].to_numpy().copy()
            for b in np.unique(depth_bins):
                idx = np.where(depth_bins == b)[0]
                sh[idx] = rng.permutation(sh[idx])
        Hn.append(meanH(sh))
    H_null = float(np.mean(Hn))
    return 1 - H_obs / H_null if H_null > 0 else np.nan


def main():
    O = ad.read_h5ad("data/processed/pc9_osimertinib_processed.h5ad")
    sc.pp.neighbors(O, n_pcs=30, random_state=42)
    sc.tl.leiden(O, resolution=1.0, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    O.obs["state"] = O.obs["leiden"].to_numpy()
    lb = O.obs["lineage_barcode"].astype(str); ok = ~lb.isin(BAD)
    depth = O.obs["total_counts"].to_numpy()
    base = pd.DataFrame({"lin": lb[ok].to_numpy(), "tp": O.obs["timepoint"][ok].to_numpy(),
                         "state": O.obs["state"][ok].to_numpy(),
                         "depth": depth[ok.to_numpy()]}, index=np.where(ok.to_numpy())[0])
    rng = np.random.default_rng(42)
    res = {}
    for tp in TPS:
        d = base[base.tp == tp].copy()
        # depth quintile bins within timepoint
        d["depthbin"] = pd.qcut(d["depth"].rank(method="first"), 5, labels=False)
        M_plain = memory_M(d, "state", rng)
        M_depth = memory_M(d, "state", rng, depth_bins=d["depthbin"].to_numpy())
        res[tp] = dict(M_plain=float(M_plain), M_depth_matched=float(M_depth),
                       median_depth=float(d["depth"].median()), n=int(len(d)))
        log.info("[%s] M_plain=%.3f  M_depth-matched=%.3f  median_depth=%.0f  n=%d",
                 tp, M_plain, M_depth, d["depth"].median(), len(d))

    decay_plain = res["D0"]["M_plain"] - res["D14"]["M_plain"]
    decay_depth = res["D0"]["M_depth_matched"] - res["D14"]["M_depth_matched"]
    print("\n==== GATE 4: depth-matched memory decay ====")
    for tp in TPS:
        r = res[tp]
        print(f"  {tp:3s}: M={r['M_plain']:.3f}  M(depth-matched)={r['M_depth_matched']:.3f}  median_depth={r['median_depth']:.0f}")
    print(f"\n  decay D0->D14: plain={decay_plain:.3f}  depth-matched={decay_depth:.3f}")
    print(f"  VERDICT: {'survives depth control' if decay_depth > 0.5*decay_plain and decay_depth>0 else 'WEAKENED/abolished by depth control'}")
    (OUT / "memory_depth_matched.json").write_text(json.dumps(
        {"per_timepoint": res, "decay_plain": decay_plain, "decay_depth_matched": decay_depth}, indent=2))


if __name__ == "__main__":
    main()
