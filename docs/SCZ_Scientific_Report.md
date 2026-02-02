# AlphaGenome Regulatory Analysis of Common Schizophrenia Risk Variants

## A Deep Learning Approach to Variant Prioritization and Mechanistic Interpretation

---

# Abstract

Genome-wide association studies have identified 287 genomic loci associated with schizophrenia risk. However, the causal variants within these loci and their downstream biological effects remain largely unknown. Here, we apply AlphaGenome—a deep neural network trained to predict tissue-specific regulatory effects—to systematically score 257 lead SNPs from the Psychiatric Genomics Consortium wave 3 (PGC3) schizophrenia GWAS.

We identify significant enrichment of intracellular calcium signaling genes (ATP2A2, ITPR3) among predicted regulatory targets (GSEA P=0.009, FDR=0.027). Digital validation against H-MAGMA 3D chromatin interaction maps confirms 22.4% recall of predicted regulatory loops. Comparison with officially prioritized PGC3 genes shows significant overlap (29/120 genes, OR=5.03, P=9.8×10⁻¹¹). Cell-type specificity analysis demonstrates significant depletion of targets in somatic tissues (OR=0.22-0.47), confirming central nervous system specificity.

This study demonstrates that deep learning models can bridge the gap between GWAS associations and biological mechanism, providing quantitative predictions for therapeutic target prioritization.

---

# 1. Introduction

## 1.1 Background

Schizophrenia (SCZ) is a severe psychiatric disorder affecting approximately 1% of the global population, with heritability estimates of 60-80% [1]. The Psychiatric Genomics Consortium wave 3 (PGC3) GWAS, published by Trubetskoy et al. (2022), identified 287 distinct genomic loci associated with schizophrenia, approximately doubling the number of independent associations from previous waves [2].

While this study prioritized 120 genes using statistical fine-mapping (FINEMAP) and expression quantitative trait loci (eQTL) integration (SMR), the biological consequences of most risk variants remain unclear. Traditional fine-mapping approaches provide probabilistic assignments of causality but do not directly predict the magnitude or tissue-specificity of regulatory effects.

## 1.2 AlphaGenome: A Deep Learning Approach

AlphaGenome is a foundation model for genomic prediction trained on large-scale functional genomics data, including ENCODE/Roadmap epigenomics, GTEx expression, and chromatin accessibility profiles from thousands of biosamples [3]. AlphaGenome directly predicts the effect of genetic variants on:

1. **Gene expression** (log-fold change per tissue)
2. **Splicing** (splice site usage changes)
3. **Chromatin accessibility** (ATAC-seq and DNase-seq)
4. **3D chromatin contacts** (Hi-C-like predictions)

## 1.3 Study Objectives

This study aims to:
1. Apply AlphaGenome to score PGC3 SCZ lead variants for regulatory effects
2. Identify pathway enrichment among predicted regulatory targets
3. Validate predictions against independent biological datasets (H-MAGMA, GTEx)
4. Compare predicted targets with officially prioritized PGC3 genes

---

# 2. Methods

## 2.1 GWAS Variant Extraction

### 2.1.1 Source Data
We obtained PGC3 schizophrenia GWAS summary statistics from the PGC data portal. The extended GWAS comprised 76,755 cases and 243,649 controls, identifying 342 linkage-disequilibrium-independent index SNPs across 287 genomic loci [2].

### 2.1.2 Lead SNP Selection
From the summary statistics file, we extracted lead SNPs using the following criteria:

```
P-value < 5 × 10⁻⁸ (genome-wide significance)
INFO score ≥ 0.9 (imputation quality)
LD-independent (r² < 0.1 within 500kb windows)
```

This yielded **257 lead SNPs** for downstream scoring.

> **Implementation**: `01_extract_gwas_variants.py`

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

### 2.2.3 Brain-Tissue Filtering
To prioritize CNS-relevant effects, we filtered expression predictions to 13 GTEx brain tissues.

> **Implementation**: `02_alphagenome_scoring.py`

---

## 2.3 Pathway Enrichment Analysis

### 2.3.1 Ranked GSEA
We performed Gene Set Enrichment Analysis using a ranked gene list approach [4]. Ranked GSEA:

1. Ranks all predicted target genes by AlphaGenome effect score
2. Tests whether pathway genes cluster toward the top of the ranking
3. Provides enrichment scores without binary cutoffs

### 2.3.2 Pathway Definitions
We curated neurobiologically-relevant pathway gene sets including:

| Pathway | N Genes | Key Members |
|---------|---------|-------------|
| Calcium_VoltageGated | 10 | CACNA1A, CACNA1C, CACNA1D |
| Calcium_Internal | 9 | ATP2A2, ITPR1, ITPR3, RYR1-3 |
| Glutamate_Ionotropic | 14 | GRIN1, GRIN2A, GRIA1-4 |

### 2.3.3 Statistical Testing
We used the Mann-Whitney U test to compare pathway gene ranks vs. non-pathway gene ranks. Effect size was quantified using rank-biserial correlation. Multiple testing correction was applied using Benjamini-Hochberg FDR.

> **Implementation**: `03_pathway_enrichment.py`

---

## 2.4 Validation Suite

### 2.4.1 H-MAGMA (3D Chromatin Architecture)
We validated AlphaGenome structural predictions against Hi-C-coupled MAGMA (H-MAGMA) annotations from adult brain tissue [5].

**Method**:
1. Load H-MAGMA SNP-gene assignments (physical chromatin interactions)
2. Extract AlphaGenome predicted SNP-gene links (score ≥ 1.0)
3. Calculate overlap (Precision and Recall)

### 2.4.2 PGC3 Gene Comparison
We compared AlphaGenome predicted targets with the 120 genes prioritized by Trubetskoy et al.

**Method**:
1. Fisher's exact test for overlap significance
2. Odds ratio calculation

### 2.4.3 Cell Type Specificity (GTEx snRNA-seq)
We tested whether AlphaGenome targets are depleted in non-brain cell types using the GTEx 8-tissue snRNA-seq atlas.

> **Implementation**: `05_hmagma_validation.py`, `06_celltype_validation.py`, `08_pgc3_comparison.py`

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

## 3.3 Validation Results

### 3.3.1 H-MAGMA (3D Chromatin)

| Metric | Value |
|--------|-------|
| AlphaGenome Predictions | 1,514 SNP-Gene pairs |
| H-MAGMA Validated | 28 pairs |
| **Recall** | **22.4%** |

#### Validated Links Include:
- rs9461856 → SYNGAP1 (synaptic Ras signaling)
- rs9687282 → CTNNA1 (cell adhesion)
- rs2596495 → HLA-B (immune)

### 3.3.2 PGC3 Gene Overlap

| Metric | Value |
|--------|-------|
| PGC3 Prioritized Genes | 120 |
| AlphaGenome Targets | 1,214 |
| **Overlap** | **29 genes (24.2%)** |
| **Odds Ratio** | **5.03** |
| **P-value** | **9.8 × 10⁻¹¹** |

**Overlapping Genes Include**:
- ATP2A2 (SERCA2 calcium pump)
- CACNA1C (voltage-gated calcium channel)
- AKT3, MAPK3 (signaling kinases)

**Key Finding**: Highly significant overlap validates AlphaGenome's ability to recover officially prioritized genes while identifying additional candidates.

### 3.3.3 Cell Type Specificity (Somatic Control)

| Cell Type | Odds Ratio | P-value |
|-----------|------------|---------|
| Fibroblast | 0.22 | 4×10⁻⁸ |
| Myocyte (NMJ) | 0.33 | 9×10⁻⁶ |
| Endothelial | 0.65 | 0.048 |
| Peripheral Neuron | 0.47 | 0.001 |

**Conclusion**: AlphaGenome targets are **significantly depleted** in somatic tissues (42/44 cell types show OR < 1), confirming CNS-specificity of predictions.

### 3.3.4 ATAC Peak Enrichment

| Cell Type | OR vs iPSC | P-value |
|-----------|-----------|---------|
| Glutamatergic | 0.27 | 0.001 |
| Dopaminergic | 0.44 | 0.017 |
| NPC | 0.40 | 0.010 |

**Key Finding**: SCZ lead SNPs are depleted in differentiated neuronal ATAC peaks relative to iPSC, suggesting variants act via distal regulatory mechanisms rather than direct promoter effects.

---

# 4. Discussion

## 4.1 Differentiation from Trubetskoy et al. (2022)

| Aspect | Trubetskoy et al. | This Study |
|--------|-------------------|------------|
| **Method** | Statistical fine-mapping + eQTL | Deep learning (AlphaGenome) |
| **Output** | Gene prioritization lists | Quantified regulatory effects |
| **Pathway Finding** | Synaptic function (GO) | Calcium_Internal (GSEA) |
| **Validation** | Rare variant convergence | H-MAGMA + Cell Type + PGC3 overlap |

### 4.1.1 Novel Contributions
1. **Quantified Effects**: AlphaGenome provides continuous scores rather than binary gene lists
2. **Intracellular Calcium**: We identify a specific sub-pathway (ATP2A2/ITPR3) not highlighted by GO enrichment
3. **Digital Validation**: We validate predictions against multiple independent datasets
4. **PGC3 Concordance**: Significant overlap (OR=5.03) with official gene list while identifying additional targets

## 4.2 Limitations

1. **AlphaGenome Training Bias**: The model was trained on adult GTEx tissues
2. **Lead SNP Only**: We scored only lead SNPs, not credible set variants
3. **MHC Region**: Several validated links were in the MHC region, which has complex LD

## 4.3 Future Directions

1. **CRISPR Validation**: Test top predictions (ATP2A2, ITPR3) in iPSC-derived neurons
2. **Drug Repurposing**: Screen for SERCA2 modulators as potential SCZ therapeutics

---

# 5. Conclusions

This study demonstrates that:

1. AlphaGenome can provide mechanistic interpretation of GWAS loci beyond statistical fine-mapping
2. Intracellular calcium signaling (ATP2A2/ITPR3) is an enriched pathway in SCZ (FDR=0.027)
3. Predictions are validated by H-MAGMA (22% recall) and PGC3 gene overlap (OR=5.03)
4. Targets show CNS specificity (depleted in 42/44 somatic cell types)

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

---

# References

1. Owen, M.J., et al. (2016). Schizophrenia. Lancet.
2. Trubetskoy, V., et al. (2022). Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature.
3. Deepmind AlphaGenome. (2024). Technical Report.
4. Subramanian, A., et al. (2005). GSEA. PNAS.
5. Sey, N.Y., et al. (2020). H-MAGMA. Nature Neuroscience.

---

# Code Availability

All analysis scripts are available at: https://github.com/dryusufcicek/SCZ_AlphaGenome_Study

Key scripts:
- `01_extract_gwas_variants.py` - GWAS extraction
- `02_alphagenome_scoring.py` - AlphaGenome scoring
- `03_pathway_enrichment.py` - Pathway enrichment (GSEA)
- `05_hmagma_validation.py` - H-MAGMA validation
- `06_celltype_validation.py` - Cell type specificity
- `08_pgc3_comparison.py` - PGC3 gene comparison

---

# Contact

- **Author**: Yusuf Cicek
- **Email**: yusuf.cicek@iuc.edu.tr
- **Institution**: Istanbul University - Cerrahpasa, Cerrahpasa Faculty of Medicine, Department of Psychiatry
