# Path 2 Research Memo — Clonal Transcriptional-State Memory Decay Statistic

**Date:** 2026-06-11
**Decision sought:** Go / no-go on a methods paper proposing a statistic that quantifies *clonal transcriptional-state memory and its decay under drug exposure*, from lineage-traced scRNA-seq.
**Anchor result (GSE150949, PC9 / osimertinib, Watermelon):** within-clone state entropy < lineage-shuffle null (p = 0.002–0.014); memory erodes under drug (entropy 0.15 @ D0 → 0.51 @ D14). Robust to dispersion + cell-cycle confounds because the lineage-shuffle null preserves both.

**Bottom line:** **GO — conditional.** Novelty risk is **MEDIUM-HIGH** on the *static* memory statistic (well-precedented) but the *time-resolved decay-under-drug* framing is an unoccupied, defensible niche. Proceed only if novelty is reframed onto the decay dynamic, the closest competitors are cited and benchmarked head-on, and generalizability is demonstrated on ≥1 independent lineage-traced drug timecourse (recommended: **MeRLin melanoma, GSE299589**). This is a Cell Systems / Genome Biology-Methods-shaped contribution, **not** a Nature-biology discovery — calibrate the target accordingly.

---

## 1. NOVELTY — is the statistic genuinely novel, and what is the defensible niche?

**Verdict: the static heritability claim is NOT novel; the decay-under-drug trajectory statistic IS defensibly novel, but narrowly so.** Every individual ingredient is published 2020–2026; the contribution lives entirely in the *combination + framing*. Treat "clones share states more than chance" and "within-clone entropy vs lineage-shuffle null" as **established**, not novel.

### What is already taken (do not claim as novel)

| Component | Prior art | Status |
|---|---|---|
| Within-vs-between-clone heritability statistic (continuous) | **Lin et al. 2025/2026** (bioRxiv 10.1101/2025.08.21.671653 → Cell Systems 2026): per-gene ANOVA ω² = variance explained by clone identity, lineage-resolved scRNA in cancer (incl. A549 lung) + stem cells | **Largely pre-empts the static statistic** |
| Quantified cell-state **heritability vs plasticity** as a named statistic | **PATH** (Schiffman/Landau, Nat Genet 2024, DOI 10.1038/s41588-024-01920-6; PMID 39317739): Moran's-I phylogenetic-autocorrelation on lineage trees, applied to cancer (pancreatic EMT, GBM, B-ALL). 1 − heritability ≈ your within-clone entropy | **Owns "quantified heritability/plasticity in cancer"** |
| Lineage-shuffle / permutation null preserving marginals | Mathew et al. 2023 (Life Sci Alliance, PMC9840405): "clonal index" = within-clone UMI-SD vs 10,000-permutation clone-label null. Also CoSpar fate-bias; PATH tree-label permutation | **Field-standard null, not a novel device** |
| Memory decays across generations (under drug) | **Iyer/Granada/Chakrabarti, PLOS Comp Biol 2025** (DOI 10.1371/journal.pcbi.1013446, PMC12469175): cousin-correlation decay (memory ~2–3 generations) + Luria-Delbrück VMR, persister fate under cisplatin (HCT116/U2OS) | **"Memory decays under drug" is established (microscopy/fate, not transcriptional entropy)** |
| Heritable expression states predict therapy resistance | MemorySeq (Shaffer/Emert 2020, Cell, PMC7496637); Harmange/Shaffer 2023 (Nat Commun s41467-023-41811-8); Nat Rev Cancer 2024 review (s41568-024-00780-w) | **"Memory genes in drug tolerance" is a named subfield** |

### Same-dataset risk (must be neutralized in framing)

- **Oren et al. 2021** (Nature, DOI 10.1038/s41586-021-03796-6; PMC9209846) **is** the source of GSE150949. Verified by full-text term search: they have **no** occurrence of "heritability," "clonal memory," "entropy," "within-clone," or "lineage shuffle." Their permutation test is on a **discrete proliferative fate** label (uni- vs multi-fate, P≈1e-5), not on continuous transcriptional state. Their clonal-state analysis correlated expression with **D14 clone size**, and that coupling **INCREASES** D0→D14 (programs sharpen). Your "entropy 0.15→0.51 decay" is the **inverse temporal framing** of an effect they qualitatively touched — defensible as a distinct statistic, but a reviewer who knows this paper will flag apparent tension. Frame explicitly: you measure **state-identity entropy within clone**, a different axis from **fate/clone-size determinism**.

### The defensible niche (the quadruple no one ships pre-packaged)

A within-clone transcriptional-state **entropy** (not ω², not energy distance, not Moran's I) computed as an explicit **time-resolved decay trajectory** across a **drug** timecourse (yielding a memory half-life / decay-rate readout), benchmarked against a **lineage-shuffle null that preserves dispersion AND cell-cycle structure**, on lineage-traced cancer. The key scientific hook: this is the **OPPOSITE** of Lin et al.'s headline (verified from their abstract: memory is "maintained by robust epigenetic mechanisms **resistant to environmental perturbations**"). Demonstrating erosion under drug is a genuine, citable tension — **if it holds and generalizes**.

### Mandatory framing rules (else novelty collapses)
1. Do **not** claim "first to quantify clonal state memory" or "first to show memory decays under drug." Both are taken.
2. **Cite Lin et al., PATH, Oren, and Iyer head-on.** Benchmark your entropy statistic against PATH heritability and Iyer's cousin-correlation/VMR on the *same* Watermelon data.
3. Lead on (a) the **decay trajectory as a perturbation-response biomarker** and (b) the **cell-cycle-preserving null** (the property that sinks naive persistent-homology approaches — but note this is a null-model property, not unique to you).
4. **Rename away from "inheritance entropy"** — that exact phrase is already taken (bioRxiv 2025.10.21.683667, a different cell-cycle-exit metric).
5. Pre-register a **monotonic decay-trend test** across D0/D3/D7/D14 rather than leaning on per-timepoint p-values (0.002–0.014 are modest given multiplicity).

### Tool landscape (corroborates white space)
None of CoSpar, CellRank/CellRank2, LineageOT, Waddington-OT, or moslin/moscot computes a within-clone state-entropy statistic, a lineage-shuffle null, or a memory-decay trajectory. They answer "where do cells go / which clone biases to which fate," not "how much heritable state-memory remains and how fast it erodes." A PubMed query for the exact concept returned zero indexed articles.

---

## 2. GENERALIZABILITY — best external dataset to validate on

**Single best pick: MeRLin melanoma — GEO GSE299589.**

- **Paper:** "Clonal dynamics shaped by diverse drug-tolerant persister states in melanoma resistance" (PMC12458950). Code: github.com/Yeqing95/MeRLin.
- **System:** human melanoma PDX (WM4237-1, WM4007, WM4380-2, …) — in vivo, different cancer, different drug class.
- **Barcode:** MeRLin expressed heritable clonal barcode (2.89M-library; 265 bp in luciferase/mNeptune2.5 3'UTR) — directly analogous to Watermelon.
- **Timepoints:** 4 longitudinal — D0 (pre-tx), D21 (early MRD), D57 (pre-recurrence), Endpoint (resistant).
- **Drug:** BRAFi/MEKi (0.3 µM dabrafenib + 30 nM trametinib).
- **Clones:** 1,127 pre-treatment → ~261 in recurrent tumors (steep selection — ideal substrate for measuring memory and its decay).
- **Access:** counts matrices AND lineage barcode tables public; code public.
- **Feasibility: 9/10.** Statistic transfers nearly verbatim (entropy-by-clone-by-timepoint vs lineage shuffle). Different drug class + in vivo PDX = strong generalizability claim.
- **One caveat to check first:** per-timepoint clone counts thin by endpoint — verify enough multi-cell clones survive at late timepoints to estimate entropy.

**Backups, in priority order:**
- **#2 GALILEO TNBC** (Nat Commun 2024, s41467-024-51424-4; PMC11366763): SUM159PT, paclitaxel, dense in-vitro reverse-time-course (D5,7,9,11,13,15) + **Multiome (ATAC+GEX)** — lets you test whether memory is epigenetically encoded. Feasibility 8/10; confirm accession is open (Data Availability behind redirect) and note severe paclitaxel bottleneck.
- **#3 SPLINTR AML** (Fennell/Dawson 2022 Nature; GSE161676): mouse MLL-AF9, in vivo chemo, third orthogonal system. Feasibility 6/10; verify ≥2 well-spaced timepoints per clone (design reads baseline-vs-endpoint).
- **Avoid:** TraCe-seq (EGA controlled access). FateMap/ClonMapper are endpoint-heavy — verify intermediate timepoints before relying on them.
- **For further mining:** scLTdb (Nucleic Acids Res 2025, scltdb.com) — 109 curated lineage-tracing datasets, filterable.

**Recommended de-risking order:** pull MeRLin (GSE299589) first; confirm late-timepoint multi-cell clone counts support entropy estimation; then re-derive the decay curve. If MeRLin reproduces erosion, add GALILEO for the cross-modality/epigenetic-encoding angle.

---

## 3. NOVELTY-RISK RATING & GO/NO-GO

**Novelty-risk rating: MEDIUM-HIGH.**
- HIGH on the static "clones share state / within-clone heritability" claim (pre-empted by Lin 2025/2026, PATH, Mathew 2023, MemorySeq).
- LOW-MEDIUM on the specific time-resolved **decay-under-drug entropy trajectory** with a cell-cycle-preserving null (unoccupied by the four trajectory tools, by Oren, and not directly published).

**Recommendation: GO, conditional on all four gates:**
1. **Reframe** novelty exclusively onto the decay/erosion dynamic + drug-tolerance link. Drop any "memory exists / decays" first-claims.
2. **Cite + benchmark head-on** vs PATH (heritability ≈ 1 − entropy), Iyer (cousin-correlation/VMR decay), Lin et al. (stable-memory thesis you contradict), and Oren (same dataset, fate axis).
3. **Generalize** on MeRLin GSE299589 before committing months — this is the load-bearing de-risk. The entire thesis currently rests on one dataset.
4. **Tighten statistics:** per-timepoint null matched for clone-size distribution and sequencing depth (the entropy rise can be partly drug-induced divergence/dropout, not loss of heritability); pre-register a monotonic decay-trend test.

**Kill criteria** (declare no-go if): MeRLin shows no erosion (signal is GSE150949-specific), OR the per-timepoint depth/clone-size-matched null abolishes the decay (confound, not memory), OR Lin et al.'s published Cell Systems version already reports a drug-perturbation decay curve.

**Target journal:** Cell Systems or Genome Biology (Methods). Not Nature Communications as a biology discovery — the contribution is a statistic + null + decay-curve biomarker, not new biology.

---

### Source ledger (primary verification)
- Oren et al. 2021, Nature, 10.1038/s41586-021-03796-6 (PMC9209846) — full-text term search confirms no memory/entropy/shuffle statistic; clone-state coupling increases D0→D14.
- Lin et al. 2025/2026, bioRxiv 10.1101/2025.08.21.671653 (Cell Systems S2405-4712(26)00103-1) — abstract verified: memory "resistant to environmental perturbations" (opposite of erosion thesis); A549 lung, not PC9/osimertinib.
- PATH, Schiffman/Landau, Nat Genet 2024, 10.1038/s41588-024-01920-6 (PMID 39317739, verified) — phylogenetic heritability statistic in cancer.
- Iyer/Granada/Chakrabarti, PLOS Comp Biol 2025, 10.1371/journal.pcbi.1013446 (PMC12469175) — drug-driven generational memory decay (fate, microscopy).
- Mathew et al. 2023, Life Sci Alliance, PMC9840405 — lineage-shuffle null + clonal index (closest null-design precedent).
- MemorySeq, Shaffer/Emert 2020, Cell, PMC7496637; Harmange/Shaffer 2023, Nat Commun s41467-023-41811-8.
- Trajectory tools: WOT 10.1016/j.cell.2019.01.006; LineageOT 10.1038/s41467-021-25133-1; CoSpar 10.1038/s41587-022-01209-1; CellRank 10.1038/s41592-021-01346-6 / CellRank2 10.1038/s41592-024-02303-9; moslin 10.1186/s13059-024-03422-4.
- Validation datasets: MeRLin GSE299589 (PMC12458950); GALILEO s41467-024-51424-4 (PMC11366763); SPLINTR GSE161676; scLTdb (Nucleic Acids Res 2025, scltdb.com).
- Name collision to avoid: "inheritance entropy," bioRxiv 2025.10.21.683667.
