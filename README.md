# AlphaGenome_SCZ

Sequence-based regulatory prioritization workflow for schizophrenia credible-set variants.

## What This Repository Contains
This repository runs a reproducible post-scoring analysis pipeline for schizophrenia GWAS credible-set variants using local AlphaGenome outputs. The current pipeline performs posterior-weighted variant-to-gene aggregation, empirical calibration, pathway enrichment, cell-type footprint testing, and orthogonal validation layers (Hi-C overlap, SCHEMA overlap, GTEx status check).

## Current Run Snapshot
The latest corrected run in this repository used 20,591 credible-set variants across 249 loci and produced 483 non-zero weighted genes after aggregation. Calibrated output contains 19,452 genes in the empirical universe. Main tables and figures are written to `results/tables` and `results/figures`.

## Pipeline Entry Point
Run from the repository root:

```bash
bash run_pipeline.sh
```

The script executes schema normalization (if backup exists), Steps 4-9, and figure/table generation.

## Key Outputs
- `data/processed/alphagenome_variant_scores.csv`
- `data/processed/gene_scores/gene_weighted_scores.csv`
- `data/processed/gene_scores/gene_z_scores.csv`
- `results/tables/gsea_enrichment_score.csv`
- `results/tables/gsea_enrichment_ALL_MODALITIES_summary.csv`
- `data/processed/celltype_enrichment/cell_type_footprint_enrichment.csv`
- `results/tables/hic_loop_overlap.csv`
- `results/tables/schema_convergence_overlap.csv`
- `docs/Nature_Manuscript_Draft.md`
- `docs/Nature_Manuscript_Supplement.md`

## Important Limitations
Some modality fields in historical AlphaGenome outputs are missing for a subset of variants. GTEx validation currently writes an explicit stub report when public file retrieval fails, rather than reporting unsupported concordance.

## License
MIT (see `LICENSE`).
