# AlphaGenome Reveals the Regulatory Architecture of Schizophrenia Risk

**Authors:** Yusuf Cicek, et al.
**Repository:** AlphaGenome SCZ Analysis

## 1. Project Overview
Schizophrenia (SCZ) polygenic risk remains poorly understood at the mechanistic level. This study applies **AlphaGenome**, a sequence-based deep learning framework, to decode the regulatory impact of 10,832 candidate causal variants prioritized from the PGC3 GWAS. By implementing a rigorous **Posterior-Weighted Credible Set** analysis and **Unbiased Whole-Genome Enrichment** strategy, we demonstrate that SCZ risk converges fundamentally on the disruption of **neuronal transcriptional regulation** ($P < 10^{-30}$), which is accompanied by specific downstream deficits in **synaptic vesicle transport** and **intracellular calcium homeostasis**.

## 2. Analysis Workflow
The analysis pipeline is modularized into conceptual steps, mapped to the `scripts/` directory:

1.  **Locus Definition (`scripts/01_locus_definition`):**
    *   Input: PGC3 GWAS Summary Statistics.
    *   Action: Define independent genomic loci ($P < 5 \times 10^{-8}$, clumping).
    *   Output: Locus BED files.

2.  **Fine-Mapping (`scripts/02_finemapping`):**
    *   Input: Defined Loci + 1000G Reference Panel.
    *   Action: Calculate Approximate Bayes Factors (ABF) and Posterior Probabilities (PP).
    *   Output: 95% Credible Sets (`scz_credible_sets_proxy.csv`).

3.  **AlphaGenome Scoring (`scripts/03_alphagenome_scoring`):**
    *   Input: Variants from Credible Sets.
    *   Action: Predict regulatory impact (Log-Fold Change) across 53 biosamples using AlphaGenome.
    *   Output: Variant-level LFCs.

4.  **Variant-to-Gene Aggregation (`scripts/04_variant_to_gene`):**
    *   Input: Variant Scores + Posterior Probabilities.
    *   Action: Aggregate scores to genes using the formula $GeneScore = \sum (PP \times LFC)$.
    *   Output: Weighted Gene Scores.

5.  **Empirical Calibration (`scripts/05_null_model`):**
    *   Input: Weighted Gene Scores + Gene Universe.
    *   Action: Calculate Z-Scores against a null distribution of 47,808 genes.
    *   Output: Calibrated Z-Scores (`gene_z_scores.csv`).

6.  **Unbiased Enrichment (`scripts/06_enrichment`):**
    *   Input: Ranked Gene List.
    *   Action: GSEA (Mann-Whitney U) against 5,000 GO terms.
    *   Output: Pathway Enrichment Tables.

7.  **Cell-Type Analysis (`scripts/07_celltype_analysis`):**
    *   Input: Genomics Footprints (ATAC-seq).
    *   Action: Footprint-aware binomial testing.
    *   Output: Normalized Cell-Type Specificity (`celltype_enrichment_normalized.csv`).

## 3. Reproducibility
### Environment Setup
To replicate the analysis environment:
```bash
conda env create -f environment.yml
conda activate scz-alphagenome
```

### Dependencies
*   Python 3.10
*   `numpy`, `pandas`, `scipy` (Core stats)
*   `statsmodels` (FDR correction)
*   `seaborn` (Visualization)

## 4. How to Reproduce Key Figures
To generate the core figures presented in the manuscript:

**Figure 2: Genome-Wide Regulatory Convergence**
```bash
python scripts/04_variant_to_gene/compute_gene_scores.py
python scripts/06_enrichment/run_gsea.py
# Outputs: Fig2A_Regulatory_Manhattan.png, Fig2B_GSEA_Dotplot.png
```

**Figure 4: Cell-Type Specificity**
```bash
python scripts/07_celltype_analysis/footprint_normalization.py
# Outputs: Fig4C_CellType_Normalized.png
```

## 5. Data Availability
*   **Included:** All processed data (credible sets, scores, enrichment results) are in `data/processed/`.
*   **External:** 
    *   PGC3 GWAS Summary Statistics: Included in `data/raw/` (or download from PGC portal).
    *   GO Database: Automatically downloaded by `run_gsea.py`.
