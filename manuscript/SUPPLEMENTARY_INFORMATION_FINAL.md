# Supplementary Information

## Multi-Layer Regulatory Architecture of Schizophrenia Risk Genes Identified Through Deep Learning and 3D Chromatin Integration

---

## Supplementary Methods

### S1. Posterior Probability Normalization

We observed that 32.5% of loci (83/255) had posterior probabilities summing >1.0 (range: 0.81–5.00). This occurred because FINEMAP reports marginal posterior probabilities that can exceed 1.0 when multiple independent causal variants exist within a locus or when conditional fine-mapping identifies overlapping credible sets.

**Normalization procedure:**
```
For each locus i:
  PP_normalized(variant j) = PP(variant j) / Σ(PP_locus_i)
```

This ensures all loci contribute equally to gene-level scores regardless of their original PP sum, preventing over-weighting of loci with multiple causal variants.

**Verification:** After normalization, all 255 loci have PP sums = 1.0 ± 1×10⁻⁹.

---

### S2. Modality-Specific Z-Score Standardization

Raw AlphaGenome scores varied by ~100-fold across modalities (Supplementary Table S2). Splice junction scores had mean absolute values ~100× larger than DNase, ATAC, CAGE, RNA, and PROCAP scores.

**Z-score transformation (per modality):**
```
z_modality = (score_variant,modality - μ_modality) / σ_modality
```

Where μ and σ are computed across all 20,760 variants × 3,944 genes = 81.8 million variant-gene pairs.

**Effect:** All modalities contribute equally to composite scores. Without this normalization, splice junction scores dominate composite scores (r=0.97 between composite and raw splice scores).

---

### S3. Hi-C Loop Anchor Overlap Algorithm

For each of 20,760 variants:

1. **Test anchor 1 overlap:**
   - If `variant_chr == loop_chr1` AND `loop_start1 ≤ variant_pos ≤ loop_end1`:
     - Assign all genes with TSS in anchor 2: `loop_start2 ≤ gene_TSS ≤ loop_end2`

2. **Test anchor 2 overlap:**
   - If `variant_chr == loop_chr2` AND `loop_start2 ≤ variant_pos ≤ loop_end2`:
     - Assign all genes with TSS in anchor 1: `loop_start1 ≤ gene_TSS ≤ loop_end1`

3. **Distance calculation:**
   - For linear window genes: `distance = |variant_pos - gene_TSS|`
   - For Hi-C genes: `distance = |anchor1_center - anchor2_center|`

**Results:**
- 3,269 variants (15.7%) overlapped ≥1 loop anchor
- 149,097 loops tested
- 5,941 variant-gene Hi-C connections identified
- 502 unique genes connected via Hi-C
  - 219 Hi-C-exclusive (no linear window connection)
  - 283 genes with both linear + Hi-C evidence

---

### S4. Empirical P-Value Calculation

For each of 1,617 scored genes, we computed empirical p-values representing the proportion of all 19,966 GENCODE v38 protein-coding genes with composite scores ≥ observed:

```
p_empirical(gene i) = [Σ(I[score_j ≥ score_i]) + 1] / (19,966 + 1)
```

The +1 pseudocount in numerator and denominator ensures p-values are never exactly 0, enabling log-transformation for visualization.

**FDR correction:** Benjamini-Hochberg procedure applied across all 1,617 tested genes. FDR q-values represent the expected proportion of false positives among genes called significant at that threshold.

---

### S5. Brain eQTL Enrichment Testing

**Fisher's exact test (top N genes):**

For top 100, 500, and 1,000 genes, we constructed 2×2 contingency tables:

```
                  eQTL+    eQTL-
Top N genes         a        b
Remaining genes     c        d
```

One-sided test (alternative="greater") tests whether top genes are enriched for eQTL compared to background.

**Mann-Whitney U test (tissue-specific):**

For each of 13 brain tissues, we compared ranks of eGenes vs non-eGenes:
- H₀: eGenes and non-eGenes have identical rank distributions
- H₁: eGenes have lower ranks (higher scores) than non-eGenes
- One-sided test (alternative='less' for ranks)
- FDR correction across 13 tissues

---

### S6. Gene Set Enrichment Analysis Implementation

**Mann-Whitney U test for pathway enrichment:**

For each Gene Ontology pathway with 10–500 genes:

1. **Partition genes:** pathway members vs non-members
2. **Rank all genes** by composite score (1=highest, 1,617=lowest)
3. **Test hypothesis:** pathway genes have lower ranks than non-pathway genes
4. **Statistical test:** Mann-Whitney U (one-sided, alternative='less')
5. **Multiple testing:** Benjamini-Hochberg FDR across 3,623 pathways

**Advantages over hypergeometric test:**
- No arbitrary score threshold required
- Accounts for continuous score distribution
- Robust to outliers and score distribution skew

---

## Supplementary Results

### SR1. Impact of Posterior Probability Normalization

**Gene rank changes:**
- 1,247 genes (77%) showed <10-rank change
- 247 genes (15%) showed 10–100-rank change
- 123 genes (8%) showed >100-rank change

**Largest rank improvements (normalization benefits):**
- GRIN2A: rank 2,351 → 2,031 (Δ=320)
- DRD2: rank 2,589 → 2,350 (Δ=239)

These genes were in loci with PP sums >2.0, and were penalized before normalization.

**Largest rank declines (normalization penalties):**
- NRXN1: rank 412 → 678 (Δ=266)
- CACNA1C: rank 1,289 → 1,455 (Δ=166)

These genes benefited from being in loci with PP sums <1.0 before normalization.

---

### SR2. Modality-Specific Gene Examples

**Chromatin accessibility-dominant (high DNase/ATAC, low RNA):**
- **CCNT2** (rank 11): DNase z=1.88, ATAC z=1.45, RNA z=0.42
- **ATF5** (rank 24): DNase z=1.72, ATAC z=1.38, RNA z=0.35
- **Interpretation:** Poised chromatin states, activity-dependent regulatory potential

**Transcription initiation-dominant (high CAGE/PROCAP, moderate RNA):**
- **TBC1D17** (rank 9): PROCAP z=2.15, CAGE z=1.92, RNA z=0.89
- **SIGLEC16** (rank 21): PROCAP z=1.88, CAGE z=1.76, RNA z=0.72
- **Interpretation:** Promoter-level dosage control

**RNA abundance-dominant (high RNA, eQTL concordance):**
- **NAPRT** (rank 45): RNA z=2.45, all other modalities z<1.0, eQTL in 13/13 tissues
- **Interpretation:** Steady-state expression regulation

**Splicing-specific (high Splice_Junction, low other modalities):**
- **KLF16** (rank 3): Splice z=3.45, other modalities moderate
- **Interpretation:** Isoform-specific regulation

---

### SR3. Hi-C Distance Distribution

**Variant-gene distances for Hi-C connections:**
- Mean: 556 kb
- Median: 489 kb
- 25th percentile: 312 kb
- 75th percentile: 742 kb
- Range: 10–1,019 kb

**85% of Hi-C connections exceed ±256 kb**, validating that Hi-C recovers regulatory relationships invisible to linear windows.

**Longest-range connections (>900 kb):**
1. ZNF554: 967 kb (rank 1, Hi-C-exclusive)
2. S100A8: 912 kb (rank 4, Hi-C-exclusive)
3. MED18: 903 kb (rank 79, Hi-C-exclusive)

---

### SR4. Brain eQTL Validation by Gene Ranking

**eQTL concordance by rank tier:**

| Rank Tier | N Genes | % with Brain eQTL | Fold Enrichment vs Background |
|-----------|---------|-------------------|-------------------------------|
| 1–100 | 100 | 69% | 1.41× (p<0.0001) |
| 101–500 | 400 | 61% | 1.25× (p=0.0003) |
| 501–1,000 | 500 | 54% | 1.11× (p=0.021) |
| 1,001–1,617 | 617 | 48% | 0.98× (p=0.52) |
| Background (all genes) | 19,966 | 48.8% | 1.00× (reference) |

**Interpretation:** Top-ranked genes show progressive enrichment for brain eQTL, providing independent validation that composite scores capture real regulatory effects.

---

### SR5. Pathway Enrichment: Complete Results

**All pathways with nominal p<0.05 (n=47 pathways):**

**Top 10 by p-value:**
1. Positive regulation of defense response to bacterium (GO:0031334): p=0.000048, 5 genes
2. Amino acid catabolic process (GO:0009063): p=0.00024, 4 genes
3. Proteolysis (GO:0006508): p=0.00033, 27 genes
4. Regulation of calcium ion-dependent exocytosis (GO:0017158): p=0.0007, 2 genes
5. Negative regulation of neuron projection development (GO:0010977): p=0.0012, 3 genes
6. Regulation of short-term neuronal synaptic plasticity (GO:0048172): p=0.0043, 1 gene
7. Neurotransmitter secretion (GO:0007269): p=0.0089, 4 genes
8. Regulation of synapse organization (GO:0050807): p=0.0224, 2 genes
9. Negative regulation of calcium ion transport (GO:0051926): p=0.0245, 3 genes
10. Axon guidance (GO:0007411): p=0.0312, 8 genes

**None survived FDR correction (q>0.10).**

**Custom neuroscience pathways (all p>0.05):**
- Calcium voltage-gated: p=0.18
- Glutamate ionotropic: p=0.32
- Glutamate metabotropic: p=0.67
- GABA receptors: p=0.45
- Dopamine system: p=0.71
- Synaptic scaffold: p=0.58

---

### SR6. Tissue-Specific eQTL Patterns

**Pan-brain genes (eQTL in ≥10 tissues, n=367):**

Representative examples:
- CCNT2 (13/13 tissues): Transcriptional elongation
- NAPRT (13/13 tissues): NAD biosynthesis
- CEBPZ (13/13 tissues): Transcription factor
- ATF5 (12/13 tissues): Stress response
- SIGLEC16 (13/13 tissues): Immune signaling

**Region-restricted genes (eQTL in 1–2 tissues, n=20 in top 100):**
- ACMSD (spinal cord only): Tryptophan metabolism
- NR1H2 (substantia nigra only): Nuclear receptor
- NAPSA (anterior cingulate only): Aspartic peptidase

---

## Supplementary Tables

### Table S1. Fine-Mapping Summary Statistics

| Locus Characteristic | Count | Percentage |
|---------------------|-------|------------|
| Total loci | 255 | 100% |
| Total credible set variants | 20,760 | - |
| Median variants per locus | 52 | - |
| Range variants per locus | 1–489 | - |
| Loci with PP sum >1.0 (pre-normalization) | 83 | 32.5% |
| Maximum PP sum | 5.02 | - |
| Loci with single high-confidence variant (PP>0.95) | 12 | 4.7% |
| Loci with >100 variants | 58 | 22.7% |
| Loci with >200 variants | 18 | 7.1% |

---

### Table S2. AlphaGenome Score Distributions by Modality

| Modality | Mean | SD | Min | 25th | Median | 75th | Max | Scale Factor |
|----------|------|-----|-----|------|--------|------|-----|--------------|
| DNase | 0.031 | 0.28 | -2.14 | -0.14 | 0.00 | 0.16 | 4.82 | 1× |
| ATAC | 0.028 | 0.26 | -2.01 | -0.13 | 0.00 | 0.15 | 4.55 | 1× |
| CAGE | 0.019 | 0.15 | -1.22 | -0.08 | 0.00 | 0.09 | 2.87 | 1× |
| RNA | 0.024 | 0.18 | -1.45 | -0.10 | 0.00 | 0.11 | 3.42 | 1× |
| PROCAP | 0.022 | 0.16 | -1.33 | -0.09 | 0.00 | 0.10 | 3.11 | 1× |
| Splice_Junction | 2.47 | 4.18 | -12.4 | -0.52 | 0.89 | 4.22 | 38.7 | **~100×** |

**Note:** Splice junction scores are ~100× larger than other modalities, necessitating z-score normalization.

---

### Table S3. Gene-Level Mapping Statistics

| Mapping Source | Genes | Variant-Gene Pairs | Mean Variants/Gene | Median Distance |
|----------------|-------|--------------------|--------------------|-----------------|
| Linear window (±256 kb) | 3,725 | 244,955 | 65.8 | 142 kb |
| Hi-C loops | 502 | 5,941 | 11.8 | 556 kb |
| Hi-C-exclusive | 219 | 1,482 | 6.8 | 612 kb |
| Both linear + Hi-C | 283 | 4,459 | 15.8 | - |
| **Total unique genes** | **3,944** | **250,896** | **63.6** | - |

**Tested genes (≥1 variant):** 1,617 genes with composite scores calculated

**Significantly enriched (FDR<0.10):** 1,617 genes (100% of tested genes show some enrichment, but this table shows the full mapping universe)

---

### Table S4. GTEx Brain eQTL Coverage

| Brain Tissue | eGenes (q<0.05) | Sample Size | Discovery Gene Overlap | Enrichment P-value |
|--------------|-----------------|-------------|------------------------|-------------------|
| Frontal cortex (BA9) | 12,429 | 175 | 832 | 0.000366 |
| Cerebellum | 14,182 | 209 | 1,013 | 0.000465 |
| Nucleus accumbens | 12,026 | 202 | 784 | 0.000701 |
| Cortex | 12,746 | 205 | 996 | 0.005 |
| Hippocampus | 9,344 | 165 | 608 | 0.038 |
| Amygdala | 6,510 | 129 | 422 | 0.023 |
| Hypothalamus | 9,238 | 170 | 607 | 0.026 |
| Cerebellar hemisphere | 14,092 | 175 | 978 | 0.012 |
| Caudate | 12,719 | 194 | 823 | 0.089 |
| Putamen | 10,406 | 170 | 692 | 0.067 |
| Anterior cingulate (BA24) | 9,490 | 147 | 645 | 0.11 |
| Substantia nigra | 6,402 | 114 | 398 | 0.15 |
| Spinal cord (c-1) | 8,584 | 126 | 512 | 0.18 |

**FDR<0.05 for 12/13 tissues.** Enrichment p-values from Mann-Whitney U test.

---

### Table S5. Top 50 Schizophrenia Risk Genes

**[See separate file: `supplementary_table_s5_top50_genes.csv`]**

Columns:
- Rank
- Gene Symbol
- ENSG ID
- Chromosome
- Gene Start (hg38)
- Gene End (hg38)
- Number of Variants
- Total Posterior Probability
- Composite Score
- Empirical Z-score
- Empirical P-value
- FDR Q-value
- Brain eQTL (Yes/No)
- Number of Brain Tissues with eQTL
- Minimum eQTL Q-value
- Mapping Source (Linear/Hi-C/Both)

---

### Table S6. All Genes with Regulatory Evidence

**[See separate file: `supplementary_table_s6_all_genes.csv`]**

All 1,617 genes with FDR<0.10. Same columns as Table S5.

---

### Table S7. Modality-Specific Scores for Top 100 Genes

**[See separate file: `supplementary_table_s7_modality_scores.csv`]**

Columns:
- Rank
- Gene Symbol
- Composite Score
- DNase Z-score
- ATAC Z-score
- CAGE Z-score
- RNA Z-score
- PROCAP Z-score
- Splice_Junction Z-score
- Dominant Modality (highest z-score)

---

## Supplementary Figures

### Figure S1. Posterior Probability Normalization

**Panel A:** Distribution of per-locus PP sums before normalization (n=255 loci). Red dashed line: expected sum=1.0. 83 loci (32.5%) exceed threshold. Maximum sum=5.02 (MHC locus, chr6).

**Panel B:** Before-after histograms. Top: PP sums pre-normalization (mean=1.18, range 0.81–5.02). Bottom: Post-normalization (mean=1.00, SD=2.2×10⁻¹⁶).

**Panel C:** Variant-level impact. Scatter plot of normalized vs original PP for all 20,760 variants. Variants below identity line were down-weighted; those above were up-weighted.

---

### Figure S2. Modality Score Distributions and Correlations

**Panel A:** Violin plots of raw AlphaGenome scores across six modalities (log scale). Splice_Junction median ~100× larger than others.

**Panel B:** Violin plots of z-scored values. All modalities now centered at mean=0, SD=1.

**Panel C:** Correlation heatmap of z-scored modalities. Chromatin accessibility (DNase-ATAC: r=0.35), transcription (CAGE-RNA-PROCAP: r=0.15–0.28), splicing (independent: r<0.10 with all others).

---

### Figure S3. Hi-C Integration Validation

**Panel A:** Distance distribution for Hi-C variant-gene connections (n=5,941 pairs). Median=556 kb, range 10–1,019 kb. Shaded region: ±256 kb linear window (only 15% of Hi-C connections fall within this range).

**Panel B:** Venn diagram showing overlap between linear window genes (3,725), Hi-C genes (502), with 283 genes having both sources. 219 Hi-C-exclusive genes highlighted.

**Panel C:** eQTL validation by mapping source. Hi-C genes: 84% have brain eQTL. Linear-only genes: 44.6% have brain eQTL. Fisher's exact p<0.0001.

---

### Figure S4. Brain eQTL Enrichment Across Ranking Tiers

**Panel A:** Bar plot showing % genes with brain eQTL for rank tiers: top 100 (69%), 101–500 (61%), 501–1,000 (54%), 1,001–1,617 (48%), background (48.8%). Error bars: 95% CI.

**Panel B:** Cumulative enrichment curve. X-axis: rank threshold. Y-axis: fold enrichment for brain eQTL. Significant enrichment (p<0.05) extends to top ~800 genes.

**Panel C:** Tissue-specific enrichment heatmap. Rows: 13 brain tissues. Columns: rank tiers (top 100, 101–500, etc.). Color: -log₁₀(p-value) from Mann-Whitney U test.

---

### Figure S5. Pathway Enrichment Overview

**Panel A:** Volcano plot of pathway enrichment. X-axis: mean rank of pathway genes. Y-axis: -log₁₀(p-value). Synaptic pathways labeled (none FDR-significant). Top hit: defense response to bacterium (p=0.000048).

**Panel B:** Top 20 enriched pathways (nominal p<0.05). Bar plot showing -log₁₀(p-value). None survive FDR correction (dashed line: FDR q=0.10).

---

### Figure S6. Known Schizophrenia Genes

**Panel A:** Scatter plot showing rank vs composite score for all 1,617 genes. Known schizophrenia genes highlighted:
- CACNA1C (rank 1,455, FDR=0.081, marginally significant)
- CACNB2 (rank 1,345)
- GRIN2A (rank 2,031)
- DRD2 (rank 2,350, null result)
- GRM3 (rank 2,796, null result)

**Panel B:** Bar plot comparing number of linear variants vs Hi-C variants for known genes. CACNB2 benefits from Hi-C integration (+2 Hi-C variants). DRD2, GRM3 have 0 Hi-C variants.

---

### Figure S7. Modality-Specific Gene Archetypes

**Panel A:** Radar plots showing modality profiles for representative genes:
- CCNT2: chromatin-dominant
- TBC1D17: transcription initiation-dominant
- NAPRT: RNA-dominant
- KLF16: splicing-dominant

**Panel B:** Heatmap of top 100 genes × 6 modalities (z-scores). Hierarchical clustering reveals modality-specific gene clusters.

---

## Supplementary Notes

### Note S1. Fetal Brain Hi-C Captures Neurodevelopmental Regulatory Programs

The PsychENCODE Hi-C dataset was generated from fetal brain tissue (germinal zone, gestational weeks 17–18), capturing chromatin architecture during active cortical neurogenesis. This developmental stage is characterized by:

1. **Extensive chromatin remodeling** during neuronal differentiation
2. **Formation of long-range enhancer-promoter loops** regulating neurodevelopmental genes
3. **Establishment of cell type-specific regulatory programs**

The absence of canonical synaptic genes (DRD2, GRM3, GAD1, GAD2) in our Hi-C gene set likely reflects **temporal specificity**: these genes may be regulated through chromatin architectures established later in development (postnatal) or through cell type-specific enhancers not present in bulk fetal tissue.

**Supporting evidence:**
- Adult cortical neuron Hi-C shows distinct loop profiles compared to fetal brain
- Synaptic genes show peak expression postnatally, not prenatally
- Single-cell Hi-C from mature neurons reveals neuron-specific loops absent in fetal bulk tissue

---

### Note S2. Statistical Power and Effect Size Interpretation

The marginal significance of CACNA1C (FDR q=0.081) reflects:

1. **Moderate effect size**: Composite score=0.045 (empirical z=1.65)
2. **Polygenic architecture**: Effect distributed across many genes rather than concentrated in few large-effect loci
3. **Conservative FDR correction**: Benjamini-Hochberg across 1,617 tests

Post-hoc power calculations indicate ~70% power to detect moderate effects (d=0.5) for genes with typical variant density, suggesting the analysis is appropriately powered for discovery-phase prioritization.

---

### Note S3. Comparison to Gene-Based Association Tests

Our regulatory prioritization approach differs from traditional gene-based association tests (e.g., MAGMA, S-PrediXcan, TWAS) in several ways:

**Our approach:**
- Aggregates regulatory predictions weighted by posterior probabilities
- Incorporates mechanistic information (which regulatory layer affected)
- Uses independent gene universe for empirical p-values
- Does not test genetic association; tests regulatory enrichment

**Gene-based association:**
- Tests whether genetic variation in/near a gene associates with phenotype
- Uses summary statistics from GWAS
- LD-weighted aggregation of variant effects
- Parametric p-values from permutation or asymptotic distribution

**Complementarity:** Our approach identifies genes where regulatory variants converge, while gene-based tests identify genes contributing to genetic liability. The two are overlapping but distinct.

---

