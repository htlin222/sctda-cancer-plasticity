# Figure Plan

Target: 6 main figures + ~15 supplementary figures.
Style: Nature Communications format (180mm max width, 7-8pt font).

---

## Main Figures

### Figure 1: Study Design and Standard Analysis

**Message:** "Here's the system and what conventional methods show"

- (a) Schematic: PC9 cells → erlotinib treatment → time points → Drop-seq → analysis
- (b) UMAP colored by time point (D0, D1, D3, D11, holiday)
- (c) UMAP colored by Leiden cluster
- (d) UMAP colored by EMT score (continuous gradient)
- (e) UMAP colored by cell cycle phase (S, G2/M, G1)
- (f) RNA velocity stream plot on UMAP (or PAGA)
- (g) Key: "Note that velocity/PAGA suggest a branching tree structure"

### Figure 2: Persistent Homology Reveals Emerging Loops

**Message:** "TDA detects topological features that change with drug treatment"

- (a) Persistence barcodes (H₀ and H₁) for each time point, side by side
- (b) Persistence diagrams for H₁ — highlight features far from diagonal
- (c) Maximum H₁ persistence as a function of treatment time (line plot)
- (d) Permutation test: observed vs. null distribution of H₁ persistence
- (e) Wasserstein distance heatmap between time points' persistence diagrams

### Figure 3: Mapper Visualization of Loop Topology

**Message:** "The loop is a real topological structure, not a visualization artifact"

- (a) Mapper graph (filter=PC1) colored by time point
- (b) Mapper graph (filter=EMT_score) colored by EMT score
- (c) Mapper graph (filter=diffusion_component) colored by time point
- (d) Zoom on loop region: nodes colored by cluster identity
- (e) Comparison: UMAP shows branch, Mapper shows loop at same region

### Figure 4: Cell Cycle Ablation

**Message:** "The loop is NOT just cell cycle"

- (a) PH barcode with all genes vs. without cell cycle genes
- (b) Persistence diagram comparison (overlay)
- (c) Mapper graph after cell cycle gene removal — loop persists (or not)
- (d) If loop persists: show that cell cycle phase is randomly distributed around loop
- (e) If partially: show residual H₁ persistence after ablation

### Figure 5: Gene Attribution and Biological Mechanism

**Message:** "The loop corresponds to a reversible EMT/MET-like transition"

- (a) Gene connectivity heatmap on Mapper graph (top 50 genes)
- (b) Circular gene expression profiles around the loop (polar plot)
- (c) GO enrichment for loop-associated genes (dot plot)
- (d) Hallmark EMT and cholesterol metabolism scores around the loop
- (e) Key transcription factors identified (bar plot or network)
- (f) Model: schematic of cyclic state transition

### Figure 6: Cross-Dataset Validation

**Message:** "This is not dataset-specific — it generalizes"

- (a) PH barcode for osimertinib-treated PC9 (GSE150949)
- (b) Persistence diagram comparison: erlotinib vs. osimertinib
- (c) Mapper graph for osimertinib dataset
- (d) Venn diagram: overlap of loop-associated genes between datasets
- (e) Enriched pathways comparison (side-by-side dot plot)
- (f) Wasserstein distance: intra-drug vs. inter-drug PH similarity

---

## Supplementary Figures

- SFig 1: QC metrics (violin plots of genes, UMIs, mito%)
- SFig 2: HVG selection and PCA variance explained
- SFig 3: Leiden resolution sensitivity
- SFig 4: Marker gene dot plot per cluster
- SFig 5: PH sensitivity to number of PCs (20, 30, 40, 50)
- SFig 6: PH sensitivity to distance metric (Euclidean vs. correlation)
- SFig 7: Mapper parameter sensitivity (n_cubes, overlap)
- SFig 8: Full permutation null distributions for all time points
- SFig 9: Subsampling stability of H₁ features
- SFig 10: Cell cycle gene expression in Mapper graph
- SFig 11: Complete gene connectivity ranked list
- SFig 12: Full GO enrichment results
- SFig 13: Transcription factor activity (decoupler) across loop phases
- SFig 14: Drug holiday topology — does loop structure persist?
- SFig 15: Patient NSCLC validation (if Phase 2 completed)

---

## Supplementary Tables

- STable 1: Dataset metadata (cells per time point, genes detected)
- STable 2: QC filter statistics
- STable 3: Persistent homology summary (all Betti numbers, all time points)
- STable 4: Permutation test results (p-values)
- STable 5: Top 100 loop-associated genes with connectivity scores
- STable 6: Full GO/Hallmark enrichment results
- STable 7: Cross-dataset gene overlap statistics

---

## Color Scheme

```python
# Consistent across all figures
COLORS = {
    'D0': '#2196F3',       # blue — untreated
    'D1': '#4CAF50',       # green — early
    'D3': '#FF9800',       # orange — DTP emerging
    'D11': '#F44336',      # red — DTEP
    'holiday': '#9C27B0',  # purple — drug holiday
}

# EMT gradient
EMT_CMAP = 'RdYlBu_r'  # red = mesenchymal, blue = epithelial

# PH diagrams
PH_COLORS = {
    'H0': '#333333',
    'H1': '#E91E63',
    'H2': '#00BCD4',
}
```
