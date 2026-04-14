# Analysis Protocol

Step-by-step protocol for reproducible analysis. Every step maps to a notebook.

---

## Step 1: Data Acquisition

```bash
make download
```

### GSE134839 (Primary)
- Source: GEO, supplementary files from Aissa et al. 2021
- Expected: count matrices for PC9 cells at Day 0, 1, 3, 11, drug holiday
- Format: likely 10x-style or Drop-seq DGE matrices
- Verify: cell counts per time point match publication

### GSE150949 (Validation)
- Source: GEO
- Expected: count matrix for PC9 cells treated with osimertinib
- ~8,139 cycling cells

---

## Step 2: Preprocessing

Notebook: `01-preprocessing.ipynb`

```python
import scanpy as sc

# Load and concatenate time points
adata = sc.read_...  # format-dependent

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='timepoint')

# Scale + PCA
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, random_state=42)

# Cell cycle scoring
s_genes = [...]   # Tirosh S-phase
g2m_genes = [...] # Tirosh G2/M-phase
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Save
adata.write('data/processed/pc9_erlotinib_processed.h5ad')
```

---

## Step 3: Standard Baseline Analysis

Notebook: `02-standard-analysis.ipynb`

Purpose: Reproduce published results. This is not the novel part — it's the control.

- Neighbors graph → UMAP → Leiden clustering
- DEG per cluster
- EMT score (scanpy gene set scoring with Hallmark EMT gene list)
- scVelo RNA velocity (if spliced/unspliced available)
- PAGA connectivity graph
- Diffusion map + pseudotime

Key output: UMAP colored by time point, cluster, EMT score, cell cycle phase.
This figure becomes Fig 1 of the manuscript.

---

## Step 4: Persistent Homology

Notebook: `03-topological-analysis.ipynb`

### 4a. Per-timepoint PH

```python
from ripser import ripser
import numpy as np

# For each time point
for tp in ['D0', 'D1', 'D3', 'D11', 'holiday']:
    cells = adata[adata.obs['timepoint'] == tp]
    X = cells.obsm['X_pca'][:, :30]  # top 30 PCs

    # Compute PH
    result = ripser(X, maxdim=1, thresh=np.inf)

    # Extract persistence diagrams
    dgm0 = result['dgms'][0]  # H₀: connected components
    dgm1 = result['dgms'][1]  # H₁: loops

    # Save diagrams
    np.save(f'results/ph/dgm0_{tp}.npy', dgm0)
    np.save(f'results/ph/dgm1_{tp}.npy', dgm1)
```

### 4b. Permutation test for H₁

```python
from persim import plot_diagrams
import numpy as np

observed_max_persistence_h1 = ...  # max(death - birth) for H₁ features

null_distribution = []
for i in range(500):
    # Permute gene labels independently per cell
    X_perm = np.copy(adata.X)
    for j in range(X_perm.shape[0]):
        np.random.shuffle(X_perm[j, :])

    # Re-run PCA on permuted data
    # ... (use sklearn PCA for speed)

    # Compute PH
    result_perm = ripser(X_pca_perm[:, :30], maxdim=1)
    if len(result_perm['dgms'][1]) > 0:
        max_pers = np.max(result_perm['dgms'][1][:, 1] - result_perm['dgms'][1][:, 0])
    else:
        max_pers = 0
    null_distribution.append(max_pers)

p_value = np.mean(np.array(null_distribution) >= observed_max_persistence_h1)
```

### 4c. Cross-timepoint comparison

- Plot persistence barcodes side by side for all time points
- Compute Wasserstein distance between persistence diagrams across time points
- Test: does H₁ persistence increase from D0 → D3 → D11?

### 4d. Mapper visualization

```python
import kmapper as km

mapper = km.KeplerMapper()

# Multiple filter functions
for filter_name, filter_func in [
    ('pc1', adata.obsm['X_pca'][:, 0]),
    ('emt_score', adata.obs['emt_score'].values),
    ('diffusion', adata.obsm['X_diffmap'][:, 0]),
]:
    projected = mapper.fit_transform(X, projection=filter_func.reshape(-1, 1))
    graph = mapper.map(projected, X, cover=km.Cover(n_cubes=15, perc_overlap=0.3))
    mapper.visualize(graph, path_html=f'results/mapper/mapper_{filter_name}.html')
```

---

## Step 5: Cell Cycle Ablation

Notebook: `04-cell-cycle-ablation.ipynb`

**This is the most critical control experiment.**

```python
# Remove cell cycle genes
cc_genes = s_genes + g2m_genes
non_cc_mask = ~adata.var_names.isin(cc_genes)
adata_nocc = adata[:, non_cc_mask].copy()

# Re-run PCA on cell-cycle-free data
sc.pp.scale(adata_nocc, max_value=10)
sc.tl.pca(adata_nocc, n_comps=50, random_state=42)

# Re-run PH
# Compare H₁ with and without cell cycle genes
```

**Interpretation:**
- H₁ disappears → loop was cell cycle. Report as: "TDA detects cell-cycle coupling
  to drug tolerance" (still publishable but less novel)
- H₁ persists → loop is independent of cell cycle. Report as: "Non-cell-cycle cyclic
  plasticity underlies drug tolerance" (HIGH IMPACT)
- H₁ changes but doesn't disappear → partial coupling. Most likely scenario.

---

## Step 6: Gene Attribution

Notebook: `05-gene-attribution.ipynb`

### 6a. Gene connectivity on Mapper graph

For each gene, compute connectivity score: does gene expression correlate with
position on the Mapper graph? High connectivity = gene is associated with a
specific topological feature.

### 6b. Loop-phase analysis

If significant H₁ loop found:
1. Identify cells participating in the loop (from representative cycles)
2. Order cells around the loop using persistent cohomology circular coordinate
3. Bin cells into 8-12 phases around the loop
4. DE analysis between adjacent phases → genes that change cyclically

### 6c. Pathway enrichment

```python
import gseapy as gp

# For loop-associated genes
enr = gp.enrichr(
    gene_list=loop_genes,
    gene_sets=['MSigDB_Hallmark_2020', 'GO_Biological_Process_2021'],
    organism='human',
    outdir='results/enrichment/'
)
```

Key gene sets to check:
- Hallmark EMT
- Hallmark cholesterol homeostasis
- Hallmark drug metabolism
- Hallmark MYC targets
- GO: epithelial cell differentiation

---

## Step 7: Validation

Notebook: `06-validation.ipynb`

Repeat Steps 2-6 on GSE150949 (osimertinib).

Key comparisons:
- Are H₁ features present in osimertinib-treated cells?
- Wasserstein distance between erlotinib and osimertinib persistence diagrams
- Do loop-associated genes overlap significantly? (hypergeometric test)
- Do enriched pathways replicate?

---

## Step 8: Figure Generation

Notebook: `07-figures.ipynb`

All figures generated programmatically with:
- matplotlib rcParams set for publication (font size 7-8pt, Arial/Helvetica)
- 300 dpi for raster, vector PDF for line art
- Consistent color scheme across all figures
- Exported to `results/figures/` and `manuscript/figures/`

See `plan/04-figure-plan.md` for detailed figure descriptions.
