# Methods

## Study Design and Rationale

We integrated fine-mapped schizophrenia GWAS variants with deep learning regulatory predictions and brain-specific 3D chromatin architecture to prioritize genes with convergent evidence of regulatory impact. The analytical framework comprised four components: (1) sequence-level regulatory scoring across six chromatin modalities using AlphaGenome, (2) gene-level aggregation incorporating both local regulatory elements (±256 kb) and long-range chromatin loops (up to 1 Mb), (3) statistical enrichment testing against an independent gene universe, and (4) validation using population-level brain expression quantitative trait loci (eQTL).

**Biological rationale**: Schizophrenia GWAS variants are predominantly non-coding, implicating regulatory mechanisms. Regulatory effects operate through multiple layers (chromatin accessibility, transcription initiation, RNA stability, splicing) and spatial scales (proximal enhancers, distal enhancers via chromatin looping). Integration of sequence-based predictions with experimentally determined chromatin architecture provides both mechanistic interpretation and spatial context.

---

## Fine-Mapped Schizophrenia Variants

Fine-mapped schizophrenia GWAS variants were obtained from the Psychiatric Genomics Consortium wave 3 (PGC3) study (Trubetskoy et al., *Nature* 2022; 76,755 cases, 243,649 controls). Fine-mapping was performed using FINEMAP software to generate 95% credible sets, yielding **20,760 variants across 255 independent loci** with posterior probability (PP) estimates for causal involvement. FINEMAP was selected due to its robust performance in high-LD regions and its use in the original PGC3 analysis, ensuring consistency between association and downstream functional modeling.

### Posterior Probability Normalization

We observed that 32.5% of loci (83/255) had posterior probabilities summing >1.0 (range: 0.81–5.00), likely due to overlapping credible sets from multiple independent signals or conditional versus marginal PP reporting. To ensure equal locus weighting in gene-level aggregation, we normalized PP within each locus:

```
PP_normalized = PP / Σ(PP_locus)
```

This normalization preserves relative variant ranking within loci while preventing loci with multiple correlated signals from disproportionately influencing gene-level scores.

---

## Deep Learning Regulatory Predictions: AlphaGenome

We used AlphaGenome (pre-trained Enformer-based deep learning model) to predict regulatory effects of each variant across **six chromatin modalities**: DNase-seq and ATAC-seq (chromatin accessibility), CAGE and PRO-CAP (transcription initiation), RNA-seq (steady-state expression), and splice junctions (alternative splicing).

For each variant, AlphaGenome computed log-fold change (LFC) scores representing predicted regulatory impact (reference vs alternative allele) across a **±256 kb window** (512 kb total). The ±256 kb window corresponds to the native AlphaGenome output range and captures the majority of proximal enhancer–promoter interactions while remaining agnostic to distal regulation recovered via Hi-C. Variant-gene assignments were extracted from AlphaGenome's "all_genes" field, which reports all genes scored within this sequence context window.

### Modality-Specific Z-Score Normalization

To correct for scale differences across modalities (raw splice junction scores were ~100× larger than RNA scores), we z-score normalized each modality independently:

```
z_modality = (score_modality - μ_modality) / σ_modality
```

This ensures equal weighting of modalities in composite score calculation, preventing dominant modalities from obscuring signal in others.

---

## Brain-Specific 3D Chromatin Architecture: Hi-C Integration

### Biological Rationale

Linear genomic windows capture only proximal regulatory elements (enhancers within ±256 kb). However, enhancers can regulate genes up to 1 Mb away through 3D chromatin looping, where distal enhancers physically contact gene promoters via spatial genome folding. To recover genes regulated by long-range interactions—particularly relevant for brain development where chromatin architecture is extensively remodeled—we integrated brain-specific Hi-C chromatin loop data.

### Data Source

We used **PsychENCODE fetal brain Hi-C loops** from Won et al. (*Science* 2016): 149,097 significant chromatin loops (FDR<0.05) in BEDPE format specifying genomic coordinates of paired interaction anchors.

**Developmental context**: Fetal brain Hi-C captures early neurodevelopmental regulatory programs active during cortical neurogenesis, migration, and circuit formation—critical periods for schizophrenia risk gene expression.

### Variant-to-Gene Mapping via Chromatin Loops

For each schizophrenia variant, we tested whether its genomic position fell within a Hi-C loop anchor:

1. **If variant in anchor 1** → Assign genes with transcription start sites (TSS) in anchor 2
2. **If variant in anchor 2** → Assign genes with TSS in anchor 1

Loop directionality was treated as undirected, consistent with current Hi-C resolution limits and established practices in gene-mapping frameworks.

This identified distal genes whose promoters physically contact variant-containing enhancers via chromatin looping. Hi-C connections ranged from 10 kb to 1,019 kb (median: 556 kb), substantially exceeding the ±256 kb linear window.

### Integration with Linear Window Assignments

Gene-level aggregation incorporated genes from both sources:
- **Linear**: Genes within ±256 kb window (proximal regulatory elements)
- **Hi-C**: Genes in loop anchors connected to variant-containing anchors (distal regulatory elements)

Each gene-variant pair was flagged by source (`linear` vs `hic`) to enable mechanistic interpretation. Genes with both linear and Hi-C evidence represent multi-scale regulatory architecture.

---

## Gene-Level Score Aggregation

### Posterior Probability-Weighted Averaging

For each gene, we aggregated variant-level z-scores using normalized posterior probabilities as weights:

```
score_gene,modality = Σ(z_variant,modality × PP_normalized) / Σ(PP_normalized)
```

This approach:
1. **Weights variants by causal probability**: High-PP variants contribute more than low-PP variants
2. **Controls for variant density**: Division by total PP prevents inflation driven by variant density
3. **Preserves modality-specific information**: Separate scores for each of six chromatin modalities

### Composite Score

The composite score for each gene represents average regulatory impact across all six modalities:

```
composite_score_gene = mean(score_DNase, score_ATAC, score_CAGE, score_RNA, score_PROCAP, score_Splice)
```

This provides an overall regulatory impact metric while preserving modality-specific signatures for biological interpretation.

---

## Statistical Testing: Empirical Null Distribution

### Independent Gene Universe

To avoid circular analysis (using the same data for discovery and background), we constructed an independent null distribution from all GENCODE v38 protein-coding genes (n=19,966). **1,617 genes had ≥1 schizophrenia variant** within ±256 kb or connected via Hi-C loop ("scored genes"), while 18,349 genes had no variant connections ("unscored genes"). Unscored genes serve as an empirical background representing expected regulatory scores under no genetic association.

### Empirical Z-Scores and P-Values

For each scored gene, we calculated empirical z-scores representing deviation from the full gene universe:

```
z_empirical = (composite_score_gene - μ_universe) / σ_universe
```

Where μ_universe and σ_universe are computed across all 19,966 genes.

Empirical p-values were computed as the proportion of universe genes with composite scores ≥ observed:

```
p_empirical = Σ(I[score_universe,i ≥ score_gene]) / n_universe
```

False discovery rate (FDR) correction was applied using the Benjamini-Hochberg procedure across all 1,617 tested genes. **Genes with FDR q<0.10 were considered significantly enriched**.

---

## Population-Level Validation: GTEx Brain eQTL

### Rationale

To validate that AlphaGenome + Hi-C predictions reflect real regulatory effects detectable in human populations, we tested for enrichment of brain expression quantitative trait loci (eQTL) in top-ranked genes.

Importantly, eQTL data were used exclusively for validation and were not incorporated into gene prioritization. Absence of bulk tissue eQTL may reflect cell type-specific regulation, developmental timing, or post-transcriptional mechanisms not captured in adult bulk brain tissue.

### Data Source

We used GTEx v10 brain eQTL data (Aguet et al., *Nature Genetics* 2020) from **13 brain tissues**: amygdala, anterior cingulate cortex (BA24), caudate, cerebellar hemisphere, cerebellum, cortex, frontal cortex (BA9), hippocampus, hypothalamus, nucleus accumbens, putamen, spinal cord (cervical c-1), and substantia nigra.

For each tissue, we extracted significant eGenes (genes with significant cis-eQTL, q<0.05 by permutation).

### eQTL Enrichment Testing

We tested whether top-ranked genes (top 100, 500, 1,000) showed enrichment for brain eQTL compared to background using Fisher's exact test (one-sided, alternative="greater"). Significant enrichment (p<0.05) validates that high-scoring genes have real regulatory effects detectable at the population level.

Additionally, we tested tissue-specific enrichment using Mann-Whitney U test: for each tissue, we compared ranks of eGenes versus non-eGenes. Significant enrichment (p<0.05, FDR-corrected) indicates brain-region specificity of regulatory effects.

---

## Pathway Enrichment Analysis

### Gene Set Enrichment Analysis (GSEA)

We performed Gene Set Enrichment Analysis using Mann-Whitney U test to test whether genes in biological pathways rank higher than expected by chance. Rank-based testing avoids reliance on arbitrary score thresholds and is robust to score distribution skew.

For each pathway:
1. Partition genes into pathway members vs non-members
2. Test whether pathway genes have significantly lower ranks (higher scores) than non-pathway genes
3. Statistical test: Mann-Whitney U (one-sided, alternative='less' for ranks)

### Pathway Databases

**Gene Ontology Biological Process 2023**: 5,407 pathways filtered to 10–500 genes (3,623 pathways tested)

**Custom neuroscience pathways** (7 pathways): Calcium voltage-gated channels (CACNA1A-S), calcium internal release (ATP2A, RYR, ITPR), glutamate ionotropic receptors (GRIN, GRIA, GRIK), glutamate metabotropic receptors (GRM1-8), GABA receptors (GABRA, GABRB, GABRG, GABBR), dopamine system (DRD1-5, SLC6A3, TH, DDC), and synaptic scaffold proteins (DLG, SHANK, HOMER).

P-values were corrected for multiple testing using Benjamini-Hochberg FDR. **Pathways with FDR q<0.10 were considered significantly enriched**.

---

## Technical Validation: H-MAGMA Comparison

To assess technical concordance of our custom Hi-C integration with established methodology, we compared results to H-MAGMA (Hi-C coupled MAGMA; Sey et al., *Nature Neuroscience* 2020).

**Validation approach**:
1. Extract all genes identified through chromatin loop mapping (n=502, "Hi-C genes")
2. Calculate concordance: (genes in both methods) / (Hi-C-only genes)
3. Test enrichment: Whether Hi-C genes rank higher in our composite scores than expected (Fisher's exact test)

**Validation criteria**:
- Concordance ≥70%: Confirms Hi-C integration captures same genes as standard method
- Enrichment p<0.05: Confirms Hi-C genes rank higher than random

High concordance and significant enrichment validate technical implementation and consistency with established methodology.

---

## Software and Statistical Analyses

All analyses were performed in Python 3.11 using pandas 2.0.3, numpy 1.24.3, scipy 1.11.1, and AlphaGenome API. Statistical significance threshold: p<0.05 (two-sided unless specified). FDR correction: Benjamini-Hochberg procedure.

---

## Data Availability

**PGC3 Schizophrenia GWAS**: Trubetskoy et al. (*Nature* 2022), Supplementary Table 11d (fine-mapped credible sets)

**AlphaGenome**: Pre-trained Enformer model via API 

**PsychENCODE Hi-C**

**GTEx v10 eQTL**

**Gene Ontology Biological Process**: Gene Ontology Consortium 2023 release

**GENCODE v38**: Protein-coding gene annotations (hg38/GRCh38)

---

## Code Availability

Analysis scripts available at: [repository to be specified upon publication]

**Key scripts**: `normalize_posterior_probs.py`, `map_variants_to_hic_genes.py`, `aggregate_scores_to_genes.py`, `compute_empirical_zscores.py`, `pathway_enrichment_gsea.py`, `gtex_eqtl_validation.py`, `run_hmagma_validation.py`
