# Nature Communications portal-ready figures

Single, standalone files for direct upload to the NC submission portal. Vector
PDF is the primary deliverable; 300 dpi PNG is included as a portal-preview
fallback. All PDFs are 1-page vector; all PNGs carry 300 ppi metadata.

| File              | Pages | PDF page size                  | PNG resolution                 |
|-------------------|-------|--------------------------------|--------------------------------|
| `Fig1.pdf/.png`   | 1     | 18.18 cm x 12.79 cm            | 2148 x 1511 px @ 300 dpi       |
| `Fig2.pdf/.png`   | 1     | 15.56 cm x 11.58 cm            | 1838 x 1368 px @ 300 dpi       |
| `Fig3.pdf/.png`   | 1     | 18.52 cm x 10.87 cm            | 2188 x 1285 px @ 300 dpi       |
| `Fig4.pdf/.png`   | 1     | 17.52 cm x 14.02 cm            | 2070 x 1657 px @ 300 dpi       |
| `Fig5.pdf/.png`   | 1     | 17.54 cm x  8.85 cm            | 2073 x 1046 px @ 300 dpi       |
| `Fig6.pdf/.png`   | 1     | 17.69 cm x  7.20 cm            | 2090 x  851 px @ 300 dpi       |
| `Fig7.pdf/.png`   | 1     | 17.52 cm x 13.37 cm            | 2070 x 1579 px @ 300 dpi       |

## Per-figure provenance

### Fig 1 (`fig:concept`) — concept-only graphical abstract
Single panel.
- copy of `manuscript/figures/fig1_concept_topology.pdf`

### Fig 2 (`fig:study_overview`) — study design and dataset map
Single panel.
- copy of `manuscript/figures/fig1_study_overview.pdf`

### Fig 3 (`fig:discovery_ph`) — discovery persistent homology, 2x2 composite
Composed via `_build/Fig3.tex` (LaTeX `standalone` + `tikz`).
- (a) `manuscript/figures/fig2_a_persistence_barcodes.pdf` — top-left, H1 barcodes across PC9 erlotinib time points
- (b) `manuscript/figures/fig2_d_wasserstein_heatmap.pdf` — top-right, pairwise H1 Wasserstein distances between time points
- (c) `manuscript/figures/fig4_a_barcodes_with_cc.pdf` — bottom-left, barcodes with cell-cycle genes retained
- (d) `manuscript/figures/fig4_b_barcodes_without_cc.pdf` — bottom-right, barcodes after cell-cycle gene ablation

### Fig 4 (`fig:validation`) — validation, vertical 2-panel
Composed via `_build/Fig4.tex`.
- (a) `results/figures/fig6_d_h1_over_time.pdf` — top, H1 persistence over time across cohorts
- (b) `results/figures/fig7_a_pdx_barcodes.pdf` — bottom, PDX H1 barcodes

### Fig 5 (`fig:master`) — master H1 comparison panel
Single panel.
- copy of `results/figures/fig8_master_h1_comparison.pdf`

### Fig 6 (`fig:harmony`) — Harmony batch correction QC
Single panel.
- copy of `results/figures/figS_harmony_batch_correction.pdf`

### Fig 7 (`fig:patient`) — patient-cohort validation, vertical 2-panel
Composed via `_build/Fig7.tex`.
- (a) `results/figures/figS_maynard_barcodes.pdf` — top, Maynard et al. patient cohort H1 barcodes
- (b) `results/figures/fig7_d_kim_atlas_barcodes.pdf` — bottom, Kim et al. lung cancer atlas H1 barcodes

## Build provenance for composites

LaTeX `standalone` + `tikz` (TinyTeX, `pdflatex`), absolute-path `\includegraphics`,
explicit cm geometry per panel, `keepaspectratio`. Source `.tex` files retained in
`_build/` for reproducibility. PNG companions rendered with `pdftoppm -r 300 -png`
then post-processed with `magick -units PixelsPerInch -density 300` so the pHYs
chunk reports 300 dpi to the portal preview.

## Verification

- `pdfinfo Fig*.pdf` confirms one page each, vector producer (`pdfTeX` or
  `Matplotlib pdf backend`), page sizes in points as listed above.
- `identify -units PixelsPerInch -format '%wx%h @ %x dpi' Fig*.png` confirms
  300 dpi on all seven PNGs.

## Notes / aspect-ratio handling

- The four barcode source PDFs (`fig2_a`, `fig4_a`, `fig4_b`, `fig7_a`,
  `fig7_d_kim_atlas`, `figS_maynard`) all share aspect 2.91:1. The composite
  templates use `keepaspectratio` and let height adjust naturally — no
  cropping or stretching.
- Fig 3 (2x2): heatmap (b) is the only non-barcode panel; it sits in a narrower
  top-right cell (6.5 cm wide), letting the wider barcode (a) at 11 cm and the
  two bottom barcodes (c, d) split the bottom row equally at 8.75 cm each.
  Final page is 18.52 cm wide -- comfortably inside NC's 18 cm double-column
  guideline -- and 10.87 cm tall.
- Bold sans-serif panel labels (a, b, c, d) are placed at each panel's
  top-left corner via TikZ `anchor=north west` + `xshift/yshift`.
- DO NOT edit source PDFs. Re-running the build is idempotent: `pdflatex
  Fig{3,4,7}.tex` in `_build/`, copy `.pdf` to this dir, then `pdftoppm -r 300`
  and `magick -units PixelsPerInch -density 300`.
