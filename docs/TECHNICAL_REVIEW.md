# Technical Review Report (Code + Data Integrity)

## Scope
A full pass was performed over pipeline scripts, key processed datasets, and manuscript-facing outputs in this repository state.

## Critical Findings Resolved
The enrichment script previously contained two `main()` implementations in one file, which allowed conflicting execution paths and stale output behavior. The script is now replaced with one implementation and explicit modality-wise outputs.

The variant-score table previously dropped 300 variants because mixed-width historical rows were discarded. Schema repair was rewritten to recover legacy and malformed rows into a canonical 12-column format, restoring coverage from 20,291 to 20,591 variants.

Variant-to-gene aggregation previously expected modality columns with different casing than the active data file. Column resolution is now case-robust, and modality values are propagated when available.

Empirical calibration previously returned only non-zero genes, which broke full-universe assumptions for downstream enrichment. The Z-score step now writes the complete gene universe with global and modality-specific z-columns.

Hi-C validation previously estimated expected overlap using potentially double-counted anchors. Anchor intervals are now merged before genomic coverage estimation.

Cell-type overlap now explicitly converts SNP coordinates (1-based) to BED coordinate space (0-based) before interval intersection.

SCHEMA convergence previously depended on hardcoded fallback genes. The script now resolves genes from local tabular sources and writes explicit failure if a gene list cannot be parsed.

GTEx validation previously had crash-prone branches and ambiguous behavior under missing files. It now always returns an explicit status report (`PARTIAL`/`NOT_COMPLETED`) without fabricated concordance values.

## Residual Risks (Not Auto-Resolvable Without New Upstream Data)
Historical AlphaGenome outputs contain missing modality values for a substantial subset of variants. The pipeline now handles this transparently, but modality-specific inference remains constrained by upstream data completeness.

GTEx cross-validation cannot be completed automatically in this run because scripted retrieval endpoints were unavailable. Coordinate-level harmonization against independently acquired eQTL files remains a pending task.

SCHEMA overlap depends on the locally available SCHEMA flag representation in extended tables. If a richer per-gene SCHEMA source is provided, convergence sensitivity should be re-estimated.

## Verification Performed
The corrected pipeline was executed end-to-end with updated scripts and produced internally consistent outputs in `data/processed`, `results/tables`, and `results/figures`.

All patched Python scripts passed AST syntax parsing.

Coverage checks confirmed `alphagenome_variant_scores.csv` now matches official credible-set SNP coverage (20,591/20,591).
