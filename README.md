# Multi-Layer Regulatory Architecture of Schizophrenia Risk Genes

**Deep Learning and 3D Chromatin Integration for Functional Genomics**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Analysis Complete](https://img.shields.io/badge/Status-Analysis%20Complete-success)]()

---

## Overview

This repository contains analysis code and results for functional genomic prioritization of schizophrenia risk genes using:
- **Fine-mapped GWAS variants** (PGC3, 20,760 variants, 255 loci)
- **Deep learning regulatory predictions** (AlphaGenome, 6 chromatin modalities)
- **3D chromatin architecture** (PsychENCODE fetal brain Hi-C)
- **Population-level validation** (GTEx v10 brain eQTL, 13 tissues)

### Key Findings
- **1,617 genes** with significant regulatory enrichment (FDR<0.10)
- **69% brain eQTL concordance** in top genes (2.33-fold enrichment, p<0.0001)
- **219 Hi-C-exclusive genes** accessible only through distal chromatin contacts
- **84% eQTL validation** for Hi-C genes (vs 44.6% linear-only)
- **Brain region specificity**: Strongest enrichment in frontal cortex, cerebellum, nucleus accumbens

---

## 📁 Repository Structure

```
SCZ_AlphaGenome_Study/
├── README.md                          # This file
├── results/                           # Final analysis results (CSV files)
├── tables/                            # Supplementary tables
├── manuscript/                        # Publication materials (SI, legends)
├── scripts/                           # Analysis pipeline
├── data/                              # Data source references
└── figures/                           # Figure specifications
```

---

## Methods Summary

### 1. Fine-Mapped Variants
- **Source**: PGC3 schizophrenia GWAS (Trubetskoy et al., Nature 2022)
- **Variants**: 20,760 variants across 255 loci
- **Fine-mapping**: FINEMAP with 95% credible sets
- **Posterior probability normalization**: Locus-wise scaling

### 2. Regulatory Predictions
- **Tool**: AlphaGenome (Enformer-based deep learning)
- **Modalities**: DNase, ATAC, CAGE, RNA, PROCAP, Splice_Junction
- **Window**: ±256 kb (512 kb total)
- **Normalization**: Z-score per modality

### 3. 3D Chromatin Integration
- **Data**: PsychENCODE fetal brain Hi-C (Won et al., Science 2016)
- **Loops**: 149,097 significant chromatin loops (FDR<0.05)
- **Gene mapping**: Variant-anchor overlap → TSS assignment
- **Distance range**: 10 kb - 1,019 kb (median: 556 kb)

### 4. Statistical Testing
- **Gene universe**: 19,966 GENCODE v38 protein-coding genes
- **Empirical p-values**: Proportion of universe genes ≥ observed score
- **FDR correction**: Benjamini-Hochberg across 1,617 tested genes
- **Threshold**: FDR q<0.10

### 5. Validation
- **eQTL data**: GTEx v10, 13 brain tissues (Aguet et al., Nat Genet 2020)
- **Enrichment tests**: Fisher's exact (top genes), Mann-Whitney U (tissue-specific)
- **Technical validation**: H-MAGMA concordance (93.2%)

---

## 📊 Key Results Files

### Main Results
| File | Description | Rows | Key Columns |
|------|-------------|------|-------------|
| `results/final_gene_results.csv` | All 1,617 significant genes | 1,617 | Gene, Score, FDR, eQTL, Modalities |
| `tables/supplementary_table_s4_top50_genes.csv` | Top 50 genes | 50 | Rank, Gene, Score, Variants, Source |
| `tables/supplementary_table_s5_all_genes.csv` | Complete gene list | 1,617 | All metrics |

### Validation Results
| File | Description |
|------|-------------|
| `results/brain_eqtl_enrichment.csv` | Tissue-specific eQTL enrichment (13 tissues) |
| `results/hic_genes.csv` | 502 Hi-C connected genes |
| `results/hic_exclusive_genes.csv` | 219 Hi-C-exclusive genes |

### Pathway Analysis
| File | Description |
|------|-------------|
| `results/pathway_enrichment_results.csv` | GSEA results (GO Biological Process) |
| `results/top_pathways_nominal.csv` | Pathways with p<0.05 (47 pathways) |

---

## Analysis Pipeline

### Prerequisites
- Python 3.11+
- Required packages: `pandas`, `numpy`, `scipy`, `statsmodels`
- AlphaGenome API access (for regulatory predictions)

### Pipeline Steps

```bash
# 1. Data Preparation
cd scripts/01_data_preparation/
python normalize_posterior_probs.py

# 2. Gene-Level Aggregation
cd ../02_gene_aggregation/
python aggregate_scores_to_genes.py

# 3. Hi-C Integration
cd ../03_hic_integration/
python map_variants_to_hic_genes.py
python run_hmagma_validation.py

# 4. Statistical Testing
cd ../04_statistical_testing/
python compute_empirical_zscores.py
python compute_fdr_qvalues.py

# 5. eQTL Validation
cd ../05_validation/
python gtex_eqtl_validation.py
python tissue_specific_enrichment.py

# 6. Pathway Analysis
cd ../06_pathway_analysis/
python pathway_enrichment_gsea.py
```

See `scripts/README.md` for detailed usage.

---

## Data Sources

All data sources are publicly available:

| Data Type | Source | Citation |
|-----------|--------|----------|
| **Schizophrenia GWAS** | PGC3 fine-mapping | Trubetskoy et al., Nature 2022 |
| **Regulatory predictions** | AlphaGenome API | Avsec et al., Nature 2026 |
| **Hi-C chromatin loops** | PsychENCODE | PsychENCODE Resources |
| **Brain eQTL** | GTEx v10 | GTEx Portal |
| **Gene annotations** | GENCODE v38 | - |
| **GO pathways** | Gene Ontology | GO Consortium 2023 |

**Data access instructions**: See `data/README.md`

---

## Key Statistics

### Gene Discovery
- **Total genes tested**: 1,617
- **Significantly enriched (FDR<0.10)**: 1,617 (100%)
- **With brain eQTL**: 1,116 (69%)
- **Hi-C genes**: 502 (31%)
- **Hi-C-exclusive**: 219 (5.6%)
- **Pan-brain genes (≥10 tissues)**: 367 (23%)

### Validation
- **Top gene eQTL enrichment**: 2.33-fold (p<0.0001)
- **Hi-C gene eQTL rate**: 84% (vs 44.6% linear-only)
- **Tissues with enrichment**: 12/13 (FDR<0.05)
- **H-MAGMA concordance**: 93.2%

### Brain Region Specificity
- **Frontal cortex (BA9)**: p=0.000366 ✓
- **Cerebellum**: p=0.000465 ✓
- **Nucleus accumbens**: p=0.000701 ✓

---

## Publication Materials

**Location**: `manuscript/`

Available files (NO manuscript document included):
- `SUPPLEMENTARY_INFORMATION_FINAL.md` - Supplementary methods, results, notes
- `FIGURE_LEGENDS_FINAL.md` - Complete figure legends (4 main + 7 extended data)
- `FINAL_ABSTRACT.md` - Manuscript abstract
- `METHODS_FINAL_REVISED.md` - Methods section only

**Note**: The manuscript document itself is not included in this repository.

---

## Citation

If you use this code or data, please cite:

```
[Citation to be added upon publication]
```

---

## Contact

**Investigator**: Yusuf Cicek, MD
**Email**: yusuf.cicek@iuc.edu.tr
**Institution**: Istanbul University-Cerrahpasa

---

## License

This project is licensed under the MIT License - see LICENSE file for details.

---

## Acknowledgments

- Psychiatric Genomics Consortium (PGC3 GWAS data)
- GTEx Consortium (brain eQTL data)
- PsychENCODE Consortium (Hi-C data)
- AlphaGenome developers

---

**Last updated**: February 2025
**Status**: Analysis complete, manuscript in preparation

