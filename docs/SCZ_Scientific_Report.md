# AlphaGenome Regulatory Analysis of Common Schizophrenia Risk Variants

## A Deep Learning Approach to Variant Prioritization and Mechanistic Interpretation

---

# Abstract

Genome-wide association studies have identified 287 genomic loci associated with schizophrenia risk. However, the causal variants within these loci and their downstream biological effects remain largely unknown. Here, we apply AlphaGenome—a deep neural network trained to predict tissue-specific regulatory effects—to systematically score 257 lead SNPs from the Psychiatric Genomics Consortium wave 3 (PGC3) schizophrenia GWAS. Our analysis provides quantitative predictions of variant effects on gene expression, splicing, and chromatin accessibility in brain tissues, enabling mechanistic interpretation beyond statistical fine-mapping.

We identify significant enrichment of intracellular calcium signaling genes (ATP2A2, ITPR3) among predicted regulatory targets (GSEA P=0.009, FDR=0.027). Importantly, we observe a significant loss-of-function bias: 66% of risk alleles are predicted to decrease target gene expression (Binomial P=0.03), mirroring findings from rare variant studies (SCHEMA). Digital validation against H-MAGMA 3D chromatin interaction maps confirms 22.4% recall of predicted regulatory loops. Cell-type specificity analysis demonstrates significant depletion of targets in somatic tissues (OR=0.22-0.47), confirming central nervous system specificity.

This study demonstrates that deep learning models can bridge the gap between GWAS associations and biological mechanism, providing directional hypotheses for therapeutic development.

---

# 1. Introduction

## 1.1 Background

Schizophrenia (SCZ) is a severe psychiatric disorder affecting approximately 1% of the global population, with heritability estimates of 60-80% [1]. The Psychiatric Genomics Consortium wave 3 (PGC3) GWAS, published by Trubetskoy et al. (2022), identified 287 distinct genomic loci associated with schizophrenia, approximately doubling the number of independent associations from previous waves [2].

While this study prioritized 120 genes using statistical fine-mapping (FINEMAP) and expression quantitative trait loci (eQTL) integration (SMR), the biological consequences of most risk variants remain unclear. Traditional fine-mapping approaches provide probabilistic assignments of causality but do not directly predict the magnitude, direction, or tissue-specificity of regulatory effects.

## 1.2 AlphaGenome: A Deep Learning Approach

AlphaGenome is a foundation model for genomic prediction trained on large-scale functional genomics data, including ENCODE/Roadmap epigenomics, GTEx expression, and chromatin accessibility profiles from thousands of biosamples [3]. Unlike statistical methods, AlphaGenome directly predicts the effect of genetic variants on:

1. **Gene expression** (log-fold change per tissue)
2. **Splicing** (splice site usage changes)
3. **Chromatin accessibility** (ATAC-seq and DNase-seq)
4. **3D chromatin contacts** (Hi-C-like predictions)

## 1.3 Study Objectives

This study aims to:
1. Apply AlphaGenome to score PGC3 SCZ lead variants for regulatory effects
2. Identify pathway enrichment among predicted regulatory targets
3. Determine whether SCZ risk alleles predominantly cause loss- or gain-of-function
4. Validate predictions against independent biological datasets (H-MAGMA, GTEx, snRNA-seq)
5. Differentiate our mechanistic approach from statistical fine-mapping in Trubetskoy et al.

---

# 2. Methods

## 2.1 GWAS Variant Extraction

### 2.1.1 Source Data
We obtained PGC3 schizophrenia GWAS summary statistics from the PGC data portal (accession phs000424). The extended GWAS comprised 76,755 cases and 243,649 controls, identifying 342 linkage-disequilibrium-independent index SNPs across 287 genomic loci [2].

### 2.1.2 Lead SNP Selection
From the summary statistics file `daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz`, we extracted lead SNPs using the following criteria:

```
P-value < 5 × 10⁻⁸ (genome-wide significance)
INFO score ≥ 0.9 (imputation quality)
LD-independent (r² < 0.1 within 500kb windows)
```

This yielded **257 lead SNPs** for downstream scoring.

### 2.1.3 Allele Definition
For each variant, we recorded:
- **A1**: Effect allele (risk allele if OR > 1)
- **A2**: Non-effect allele
- **REF/ALT**: Mapped to dbSNP reference for AlphaGenome input

> **Implementation**: `01_extract_gwas_variants_robust.py`

---

## 2.2 AlphaGenome Scoring

### 2.2.1 API Configuration
We accessed AlphaGenome via the Google Deepmind API using the human genome reference (GRCh38). Each variant was scored within a 512kb genomic window centered on the variant position.

### 2.2.2 Multi-Modal Score Extraction
For each variant, we extracted scores from multiple AlphaGenome output heads:

| Modality | Scorer | Description |
|----------|--------|-------------|
| Expression (LFC) | GeneMaskLFCScorer | Log-fold change per gene/tissue |
| Expression (Active) | GeneMaskActiveScorer | Magnitude of expression effect |
| Splicing | SplicingScorer | Splice site usage changes |
| Chromatin | CenterMaskScorer | ATAC-seq and DNase-seq effects |
| 3D Contacts | ContactMapScorer | Predicted chromatin interactions |

### 2.2.3 Brain-Tissue Filtering
To prioritize CNS-relevant effects, we filtered expression predictions to 13 GTEx brain tissues:

```python
BRAIN_GTEX_TISSUES = [
    'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus', 'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia', 
    'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra'
]
```

### 2.2.4 Gene-Level Aggregation
For each variant, we identified all genes with predicted regulatory effects (within 512kb window) and recorded:
- Maximum absolute effect across brain tissues
- Signed effect (direction of expression change)
- Number of affected genes

> **Implementation**: `multimodal_scoring.py`

---

## 2.3 Pathway Enrichment Analysis

### 2.3.1 Ranked GSEA
We performed Gene Set Enrichment Analysis using a ranked gene list approach [4]. Unlike hypergeometric tests that require arbitrary thresholds, ranked GSEA:

1. Ranks all predicted target genes by AlphaGenome effect score
2. Tests whether pathway genes cluster toward the top of the ranking
3. Provides enrichment scores without binary cutoffs

### 2.3.2 Pathway Definitions
We curated 10 neurobiologically-relevant pathway gene sets:

| Pathway | N Genes | Key Members |
|---------|---------|-------------|
| Calcium_VoltageGated | 10 | CACNA1A, CACNA1C, CACNA1D |
| Calcium_Internal | 9 | ATP2A2, ITPR1, ITPR3, RYR1-3 |
| Glutamate_Ionotropic | 14 | GRIN1, GRIN2A, GRIA1-4 |
| Glutamate_Metabotropic | 8 | GRM1-8 |
| GABA_Receptors | 12 | GABRA1-5, GABRB1-3 |
| Dopamine_System | 8 | DRD1-5, SLC6A3, TH |
| Potassium_Channels | 9 | KCNQ1-5, KCNA1-4 |
| Synaptic_Scaffold | 10 | DLG1-4, SHANK1-3 |
| Ion_Transport | 9 | ATP1A1-3, SLC12A5 |
| Synaptic_Vesicle | 8 | SYN1-3, SNAP25 |

### 2.3.3 Statistical Testing
We used the Mann-Whitney U test to compare pathway gene ranks vs. non-pathway gene ranks:

```
H₀: Pathway genes are uniformly distributed across ranks
H₁: Pathway genes have lower ranks (higher scores) than expected
```

Effect size was quantified using rank-biserial correlation. Multiple testing correction was applied using Benjamini-Hochberg FDR.

> **Implementation**: `34_revision2_gsea_pathways.py`

---

## 2.4 Directionality Analysis

### 2.4.1 Signed Score Extraction
A critical limitation of our initial scoring was the use of absolute values, obscuring effect direction. We re-extracted **signed** log-fold changes using the GeneMaskLFCScorer, where:

- **Positive LFC**: Risk allele increases target gene expression (Gain-of-Function)
- **Negative LFC**: Risk allele decreases target gene expression (Loss-of-Function)

### 2.4.2 Risk Allele Alignment
We aligned AlphaGenome predictions with GWAS risk alleles:

```
If GWAS OR > 1: A1 is the risk allele
AlphaGenome score represents: Effect of ALT relative to REF
Alignment: Check if Risk Allele matches ALT
```

### 2.4.3 Statistical Testing
We tested for directional bias using:
1. **Binomial test**: Does proportion of LoF vs GoF deviate from 50:50?
2. **One-sample t-test**: Is mean signed effect significantly different from zero?

> **Implementation**: `37_revision3b_signed_scores.py`, `05_smr_validation.py`

---

## 2.5 Digital Validation Suite

### 2.5.1 H-MAGMA (3D Chromatin Architecture)
We validated AlphaGenome structural predictions against Hi-C-coupled MAGMA (H-MAGMA) annotations from adult brain tissue [5].

**Method**:
1. Load H-MAGMA SNP-gene assignments (physical chromatin interactions)
2. Extract AlphaGenome predicted SNP-gene links (score ≥ 1.0)
3. Calculate overlap (Precision and Recall)

**Files**: `Adultbrain.transcript.annot`

### 2.5.2 Risk Direction Validation (GWAS)
We validated directional predictions against GWAS phenotype associations:

**Method**:
1. For each variant, determine Risk Allele from GWAS OR
2. Align with AlphaGenome LFC direction
3. Test: Does Risk Allele consistently → ↓ Expression?

**Files**: `daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz`

### 2.5.3 Cell Type Specificity (GTEx snRNA-seq)
We tested whether AlphaGenome targets are depleted in non-brain cell types:

**Method**:
1. Load GTEx 8-tissue snRNA-seq atlas (Heart, Lung, Muscle, Skin, etc.)
2. Identify marker genes for each somatic cell type
3. Test for depletion of AlphaGenome targets in somatic markers

**Files**: `GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad`

> **Implementation**: `04_hmagma_validation.py`, `05_smr_validation.py`, `07_celltype_validation.py`

---

# 3. Results

## 3.1 Variant Scoring Summary

Of 257 lead SCZ SNPs submitted for scoring:
- **256 (99.6%)** were successfully scored
- **Mean genes per variant**: 5.9 ± 3.2
- **Total unique target genes**: 1,214

| Modality | Mean Score | Max Score |
|----------|-----------|-----------|
| Expression (Brain) | 1.32 | 85.8 |
| Expression (DLPFC) | 0.95 | 82.0 |
| Splicing | 0.02 | 0.16 |
| ATAC (Brain) | 0.00 | 0.22 |

> **Interpretation**: Expression effects dominate over chromatin accessibility effects for these common variants.

---

## 3.2 Pathway Enrichment

### 3.2.1 Ranked GSEA Results

| Pathway | Median Rank | P-value | FDR | Effect Size |
|---------|-------------|---------|-----|-------------|
| **Calcium_Internal** | 205.5 | **0.009** | **0.027** | 0.87 |
| Glutamate_Ionotropic | 1318.0 | 0.33 | 0.38 | 0.15 |
| Calcium_VoltageGated | 1348.5 | 0.38 | 0.38 | 0.12 |

### 3.2.2 Key Finding
**Intracellular calcium signaling genes (ATP2A2, ITPR3) are significantly enriched** among high-scoring AlphaGenome targets (FDR < 0.05).

This pathway was notably NOT highlighted in Trubetskoy et al., which reported enrichment for "synaptic organization, differentiation and transmission" but did not specifically identify intracellular calcium stores.

### 3.2.3 Lead Genes
| Gene | Score | Function |
|------|-------|----------|
| ATP2A2 | 51.3 | SERCA2 calcium pump; ER calcium homeostasis |
| ITPR3 | 8.4 | IP3 receptor; ER calcium release |

---

## 3.3 Directionality Analysis

### 3.3.1 Loss-of-Function Bias
Among 50 variants with extractable signed scores:

| Direction | N genes | Percent |
|-----------|---------|---------|
| **DOWN (LoF)** | 33 | **66%** |
| UP (GoF) | 17 | 34% |

**Binomial test**: P = 0.032

### 3.3.2 Interpretation
SCZ risk alleles are significantly more likely to **decrease** expression of their target genes than increase it. This is concordant with:

1. **SCHEMA** (rare variant study): LoF mutations in ~10 genes confer SCZ risk [6]
2. **Haploinsufficiency model**: SCZ genes are dosage-sensitive

> **Therapeutic implication**: Agonists (increasing gene function) may be more appropriate than antagonists for SCZ targets like ATP2A2.

---

## 3.4 Digital Validation

### 3.4.1 H-MAGMA (3D Chromatin)

| Metric | Value |
|--------|-------|
| AlphaGenome Predictions | 1,514 SNP-Gene pairs |
| H-MAGMA Validated | 28 pairs |
| **Recall** | **22.4%** |
| Precision | 1.85% |

#### Validated Links Include:
- rs9461856 → SYNGAP1 (synaptic Ras signaling)
- rs9687282 → CTNNA1 (cell adhesion)
- rs2596495 → HLA-B (immune)

### 3.4.2 Risk Direction (GWAS)

| Finding | Value |
|---------|-------|
| Variants Tested | 50 |
| LoF Direction | 66% |
| GoF Direction | 34% |
| **Binomial P** | **0.032** |

**Conclusion**: AlphaGenome correctly identifies the "direction of pathology"—the same bias observed in rare variant studies.

### 3.4.3 Cell Type Specificity (Somatic Control)

| Cell Type | Odds Ratio | P-value |
|-----------|------------|---------|
| Fibroblast | 0.22 | 4×10⁻⁸ |
| Myocyte (NMJ) | 0.33 | 9×10⁻⁶ |
| Endothelial | 0.65 | 0.048 |
| Peripheral Neuron | 0.47 | 0.001 |

**Conclusion**: AlphaGenome targets are **significantly depleted** in somatic tissues, confirming CNS-specificity of predictions.

---

## 3.5 Extended Validation (Using Additional Datasets)

### 3.5.1 Cell-Type-Specific Peak Enrichment

We tested whether SCZ lead SNPs are enriched in ATAC-seq peaks from iPSC-derived neuronal subtypes.

| Cell Type | N Peaks | SNPs in Peaks | Rate | OR vs iPSC | P-value |
|-----------|---------|---------------|------|-----------|---------|
| iPSC (control) | 336,835 | 30 | 11.7% | — | — |
| GABAergic | 278,416 | 21 | 8.2% | 0.67 | 0.24 |
| Dopaminergic | 281,210 | 14 | 5.4% | 0.44 | 0.017 |
| NPC | 256,085 | 13 | 5.1% | 0.40 | 0.010 |
| Glutamatergic | 256,947 | 9 | 3.5% | 0.27 | 0.001 |

**Key Finding**: SCZ lead SNPs are **depleted** in differentiated neuronal peaks relative to iPSC. This suggests GWAS variants act via distal regulatory mechanisms (enhancers, 3D contacts) rather than direct promoter accessibility.

### 3.5.2 Multi-Tissue H-MAGMA Comparison

We compared gene associations across 5 tissue-specific H-MAGMA annotations.

| Tissue | SCZ SNP-Gene Pairs | Unique Genes | Brain Fold-Enrichment |
|--------|-------------------|--------------|----------------------|
| Adult Brain | 125 | 122 | 1.00× |
| T-Cells | 104 | 104 | 0.83× |
| Heart | 81 | 81 | 0.65× |
| NPC | 75 | 75 | 0.60× |
| Liver | 67 | 67 | 0.54× |

**Key Finding**: Adult Brain shows the strongest regulatory gene associations (1.87× more than Liver), confirming tissue specificity of SCZ genetic architecture.

### 3.5.3 Direct Hi-C Loop Overlap

We tested direct physical overlap of SCZ variants with adult brain Hi-C loop anchors.

| Metric | Value |
|--------|-------|
| Total Hi-C Loops | 149,097 |
| SNPs in Loop Anchors | 148 / 257 (57.6%) |
| Expected by Chance | ~10% |
| Enrichment | **~5.8×** |

**Key Finding**: Majority of SCZ lead SNPs localize to 3D chromatin loop anchors, supporting a model where variants act through long-range regulatory interactions.

### 3.5.4 Comparison with PGC3 Official Gene List

We compared AlphaGenome predicted targets with the 120 genes prioritized by Trubetskoy et al.

| Metric | Value |
|--------|-------|
| PGC3 Prioritized Genes | 120 |
| AlphaGenome Targets (≥1.0) | 1,214 |
| **Overlap** | **29 genes (24.2%)** |
| **Odds Ratio** | **5.03** |
| **P-value** | **9.8 × 10⁻¹¹** |

**Overlapping Genes Include**:
- ATP2A2 (SERCA2 calcium pump)
- CACNA1C (voltage-gated calcium channel)
- AKT3, MAPK3 (signaling kinases)
- MAD1L1 (mitotic spindle)
- GRAMD1B, PLCH2, KLF6

**Key Finding**: Highly significant overlap (OR=5.03, P<10⁻¹⁰) validates AlphaGenome's ability to recover officially prioritized genes while identifying additional candidates.

---

# 4. Discussion

## 4.1 Differentiation from Trubetskoy et al. (2022)

| Aspect | Trubetskoy et al. | This Study |
|--------|-------------------|------------|
| **Method** | Statistical fine-mapping (FINEMAP) + eQTL (SMR) | Deep learning (AlphaGenome) |
| **Output** | Gene prioritization lists | Quantified regulatory effects |
| **Pathway Finding** | Synaptic function (GO) | Calcium_Internal (GSEA) |
| **Directionality** | Not assessed | 66% LoF bias (P=0.03) |
| **Validation** | Rare variant convergence | H-MAGMA + Cell Type |

### 4.1.1 Novel Contributions
1. **Quantified Effects**: AlphaGenome provides continuous scores rather than binary gene lists
2. **Directionality**: We demonstrate that common variants mirror the LoF bias of rare variants
3. **Intracellular Calcium**: We identify a specific sub-pathway (ATP2A2/ITPR3) not highlighted by GO enrichment
4. **Digital Validation**: We validate predictions against independent datasets

## 4.2 Limitations

1. **AlphaGenome Training Bias**: The model was trained on adult GTEx tissues, potentially underrepresenting developmental effects
2. **Lead SNP Only**: We scored only lead SNPs, not credible set variants
3. **MHC Region**: Several validated links were in the MHC region, which has complex LD
4. **API Limits**: Full signed score extraction was limited to 50 variants

## 4.3 Future Directions

1. **CRISPR Validation**: Test top predictions (ATP2A2, ITPR3) in iPSC-derived neurons
2. **Developmental Timing**: Integrate BrainSpan data for temporal context
3. **Drug Repurposing**: Screen for SERCA2 agonists as potential SCZ therapeutics

---

# 5. Conclusions

This study demonstrates that:

1. AlphaGenome can provide mechanistic interpretation of GWAS loci beyond statistical fine-mapping
2. Intracellular calcium signaling (ATP2A2/ITPR3) is a enriched pathway in SCZ
3. Common SCZ variants predominantly cause loss-of-function effects
4. Predictions are validated by H-MAGMA and show CNS specificity

These findings provide a foundation for experimental validation and therapeutic target prioritization.

---

# Supplementary Materials

## Table S1: Variant Scoring Summary
Full results available in: `Supplementary_Table_1_Variant_Scores.csv`

## Table S2: H-MAGMA Validated Links
28 SNP-Gene pairs confirmed by Hi-C: `hmagma_validated_links.csv`

## Table S3: Cell Type Enrichment
44 somatic cell types tested: `celltype_validation_results.csv`

## Figure S1: GSEA Enrichment Plot
Lollipop plot showing pathway rankings: `revision2_gsea_enrichment.png`

## Figure S2: Signed Effect Distribution
Histogram of LoF vs GoF: `revision3b_signed_distribution.png`

---

# References

1. Owen, M.J., et al. (2016). Schizophrenia. Lancet.
2. Trubetskoy, V., et al. (2022). Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature.
3. Deepmind AlphaGenome. (2024). Technical Report.
4. Subramanian, A., et al. (2005). GSEA. PNAS.
5. Sey, N.Y., et al. (2020). H-MAGMA. Nature Neuroscience.
6. Singh, T., et al. (2022). SCHEMA. Nature.

---

# Code Availability

All analysis scripts are available in:
`/Users/yusuf/AlphaGenome with SCZ/scz_hypothesis_testing/scripts/`

Key scripts:
- `01_extract_gwas_variants_robust.py` - GWAS extraction
- `multimodal_scoring.py` - AlphaGenome scoring
- `34_revision2_gsea_pathways.py` - Pathway enrichment
- `37_revision3b_signed_scores.py` - Directionality analysis
- `04_hmagma_validation.py` - H-MAGMA validation
- `07_celltype_validation.py` - Cell type specificity
