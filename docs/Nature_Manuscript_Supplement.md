# Supplementary Information

**Title:** Sequence-Based Regulatory Scoring of Schizophrenia Credible-Set Variants with AlphaGenome  
**Repository:** `AlphaGenome_SCZ`

## Supplementary Note 1: Data Sources and Provenance
This analysis uses local schizophrenia fine-mapping data derived from the PGC3 framework (`data/processed/credible_sets/scz_credible_sets_official.csv`). The file contains 20,591 variants in 249 loci and provides per-variant posterior probabilities used for downstream weighting. Sequence-based regulatory predictions were obtained from local AlphaGenome output tables and normalized into a canonical schema before re-analysis to prevent row-loss caused by mixed historical output formats.

Cell-type accessibility resources were consumed from local external files corresponding to adult brain cluster peaks (Corces dataset naming convention: `Cluster*.idr.optimal.narrowPeak.gz`) and fetal consensus peaks (`GSE162170_atac_consensus_peaks.bed.gz`). Hi-C validation used the local adult brain BEDPE file (`adultbrain_hic.bedpe`) resolved from the shared validation directory. SCHEMA convergence used the SCHEMA flag column from the local extended PGC3 table in `data/validation/scz2022-Extended-Data-Table1.xlsx`.

GTEx brain eQTL validation was configured as best-effort download and parse. In the current execution, scripted URLs returned HTTP 404, so the pipeline wrote an explicit stub report (`results/tables/gtex_validation_stub.csv`) and did not report fabricated concordance statistics.

## Supplementary Note 2: Schema Normalization and Variant Coverage Recovery
Before correction, the active variant-score file lost 300 variants relative to the official credible-set input because mixed-width rows were dropped during a prior schema-repair pass. The revised schema normalization script (`scripts/fix_schema.py`) now parses legacy 5-column rows, 12-column rows, and malformed 13-column rows into a canonical 12-column table with standard modality headers. After this repair, `data/processed/alphagenome_variant_scores.csv` recovered full SNP coverage (20,591/20,591 variants).

The repaired file preserves all variants for global weighted analyses. Modality values are still unavailable for a subset of variants due historical upstream outputs; this limitation is propagated transparently into downstream modality analyses.

## Supplementary Note 3: Gene Aggregation and Empirical Calibration
Gene-level burden is defined as posterior-weighted accumulation of variant effects at assigned target genes. Aggregation is implemented in `scripts/04_variant_to_gene/compute_gene_scores.py` and outputs `data/processed/gene_scores/gene_weighted_scores.csv`. In this run, 483 genes had non-zero weighted burden.

Empirical calibration against the full local gene universe is implemented in `scripts/05_null_model/compute_empirical_zscores.py`. The calibrated matrix (`data/processed/gene_scores/gene_z_scores.csv`) includes 19,452 genes, preserving both global `z_score` and modality-specific z-columns. The weighted-score distribution figure is available at `results/figures/score_distribution_z.png`.

## Supplementary Note 4: Pathway Enrichment Implementation
Pathway testing is implemented in `scripts/06_enrichment/run_gsea.py` using one-sided Mann-Whitney U tests against non-pathway genes with Benjamini-Hochberg correction. The global weighted-score results are stored in `results/tables/gsea_enrichment_score.csv`; modality-stratified results are stored in separate files (`gsea_enrichment_DNase.csv`, `gsea_enrichment_H3K27ac.csv`, and others), and a combined top-50 summary per modality is stored in `results/tables/gsea_enrichment_ALL_MODALITIES_summary.csv`.

## Supplementary Note 5: Cell-Type Footprint Normalization
Cell-type overlap is implemented in `scripts/07_celltype_analysis/footprint_normalization.py` by merging peak intervals per class, estimating genomic footprint fractions, and testing observed variant overlaps against binomial expectation. Coordinate matching now explicitly converts 1-based SNP positions to the 0-based BED convention during overlap checks.

The resulting table (`data/processed/celltype_enrichment/cell_type_footprint_enrichment.csv`) reports seven classes. In this run, fetal progenitors reached nominal significance (P = 0.010, enrichment = 1.04), whereas mature excitatory and inhibitory classes were not significant.

## Supplementary Note 6: Orthogonal Validation Layers
Hi-C overlap (`scripts/09_validation/hic_chromatin_loops.py`) uses merged loop anchors to avoid overestimation from overlapping anchors. The hit table (`results/tables/hic_loop_overlap.csv`) contains 9,321 variants.

SCHEMA convergence (`scripts/09c_convergence_schema.py`) calculates hypergeometric overlap between top weighted genes and SCHEMA-flagged genes resolved from local extended data. The overlap table (`results/tables/schema_convergence_overlap.csv`) contains GRIN2A and SP4 in this run, with P = 0.00364.

GTEx validation (`scripts/09_validation/gtex_eqtl_validation.py`) currently returns a transparent stub status because downloadable brain files were not resolved from scripted URLs.

## Supplementary Tables Produced in This Run
`results/tables/Supplementary_Table_1_Variant_Scores.csv` contains the merged credible-set and AlphaGenome variant-level matrix (20,591 rows). `results/tables/Supplementary_Table_2_Gene_Scores.csv` contains calibrated genome-wide gene scores (19,452 rows). `results/tables/Supplementary_Table_3_Pathway_GSEA.csv` contains full global pathway enrichment output (2,717 rows). `results/tables/Supplementary_Table_4_CellType_Enrichment.csv` contains footprint-normalized cell-type enrichment values (7 rows). `results/tables/Supplementary_Table_5_MultiModal_GSEA_Top50.csv` contains the top-50 pathways for each of eight score spaces (400 rows). `results/tables/Supplementary_Table_6_SCHEMA_Overlap.csv` contains convergent genes from SCHEMA overlap (2 rows). `results/tables/Supplementary_Table_7_HiC_Loop_Overlap.csv` contains variant-level Hi-C anchor overlaps (9,321 rows).

## References
Trubetskoy V, Pardiñas AF, Qi T, et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature*. 2022;604:502-508. doi:10.1038/s41586-022-04434-5.

Singh T, Poterba T, Curtis D, et al. Rare coding variants in ten genes confer substantial risk for schizophrenia. *Nature*. 2022;604:509-516. doi:10.1038/s41586-022-04556-w.

Avsec Z, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature*. 2026;649. doi:10.1038/s41586-025-10014-0.

Corces MR, Shcherbina A, Kundu S, et al. Single-cell epigenomic analyses implicate candidate causal variants at inherited risk loci for Alzheimer's and Parkinson's diseases. *Nat Genet*. 2020;52(11):1158-1168. doi:10.1038/s41588-020-00721-x.

Trevino AE, Müller F, Andersen J, et al. Chromatin and gene-regulatory dynamics of the developing human cerebral cortex at single-cell resolution. *Cell*. 2021;184(19):5053-5069.e23. doi:10.1016/j.cell.2021.07.039.

GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776.
