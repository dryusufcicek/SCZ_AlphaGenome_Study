# AlphaGenome Regulatory Analysis of Schizophrenia GWAS Variants

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A computational pipeline for predicting regulatory effects of schizophrenia risk variants using the AlphaGenome deep learning model.

## Overview

This repository contains the analysis code for:

> **Deep Learning-Based Regulatory Effect Prediction of Common Schizophrenia Risk Variants**
> 
> We apply AlphaGenome to score 257 lead SNPs from the PGC3 schizophrenia GWAS, identifying significant enrichment of intracellular calcium signaling (ATP2A2, ITPR3) and validating predictions against chromatin interaction maps.

## Key Findings

| Finding | Value | Significance |
|---------|-------|--------------|
| Pathway Enrichment | Calcium_Internal | FDR = 0.027 |
| PGC3 Gene Overlap | 29/120 (24.2%) | P = 9.8×10⁻¹¹ |
| H-MAGMA Validation | 22% recall | 28 SNP-gene pairs |
| Somatic Depletion | 42/44 cell types | OR = 0.22-0.47 |

## Repository Structure

```
SCZ_AlphaGenome_Study/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── config.py                          # Configuration and API keys
├── scripts/
│   ├── 01_extract_gwas_variants.py   # GWAS data extraction
│   ├── 02_alphagenome_scoring.py     # AlphaGenome API scoring
│   ├── 03_pathway_enrichment.py      # Ranked GSEA analysis
│   ├── 05_hmagma_validation.py       # H-MAGMA validation
│   ├── 06_celltype_validation.py     # Cell-type specificity
│   ├── 07_hic_loop_overlap.py        # Hi-C loop analysis
│   └── 08_pgc3_comparison.py         # Official gene list comparison
├── data/
│   ├── raw/                          # Original GWAS files (not tracked)
│   └── processed/                    # Processed variant lists
├── results/
│   ├── figures/                      # Publication figures
│   └── tables/                       # Supplementary tables
└── docs/
    └── SCZ_Scientific_Report.md      # Full methods and results
```

## Installation

```bash
# Clone repository
git clone https://github.com/dryusufcicek/SCZ_AlphaGenome_Study.git
cd SCZ_AlphaGenome_Study

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Set AlphaGenome API key
export ALPHAGENOME_API_KEY="your_api_key_here"
```

## Usage

### Step 1: Extract GWAS Variants
```bash
python scripts/01_extract_gwas_variants.py
```

### Step 2: Score Variants with AlphaGenome
```bash
python scripts/02_alphagenome_scoring.py
```

### Step 3: Run Pathway Enrichment
```bash
python scripts/03_pathway_enrichment.py
```

### Step 4: Validate Results
```bash
python scripts/05_hmagma_validation.py
python scripts/06_celltype_validation.py
```

## Data Sources

| Dataset | Source | Reference |
|---------|--------|-----------|
| PGC3 SCZ GWAS | [PGC Data Portal](https://pgc.unc.edu) | Trubetskoy et al., 2022 |
| H-MAGMA Annotations | [H-MAGMA GitHub](https://github.com/thewonlab/H-MAGMA) | Sey et al., 2020 |
| GTEx snRNA-seq | [GTEx Portal](https://gtexportal.org) | GTEx Consortium |

## Citation

If you use this code, please cite:

```bibtex
@article{scz_alphagenome_2026,
  title={Deep Learning-Based Regulatory Effect Prediction of Common Schizophrenia Risk Variants},
  author={Cicek, Yusuf},
  year={2026},
  note={Manuscript in preparation}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Yusuf Cicek
- **Email**: yusuf.cicek@iuc.edu.tr
- **Institution**: Istanbul University - Cerrahpasa, Cerrahpasa Faculty of Medicine, Department of Psychiatry

## Acknowledgments

- Psychiatric Genomics Consortium for GWAS data
- Google DeepMind for AlphaGenome API access
- Won Lab for H-MAGMA annotations
