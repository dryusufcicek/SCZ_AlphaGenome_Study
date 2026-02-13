#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

echo "======================================================================"
echo "AlphaGenome SCZ: Analysis Pipeline"
echo "======================================================================"

# Optional schema repair if mixed historical rows exist
if [ -f "data/processed/alphagenome_variant_scores.csv.bak_schema_fix" ]; then
  echo "[PREP] Normalizing variant score schema..."
  python scripts/fix_schema.py
fi

# Step 4: Gene Aggregation (Multi-Modal)
echo "[STEP 4] Aggregating variant scores..."
python scripts/04_variant_to_gene/compute_gene_scores.py

# Step 5: Empirical Z-scoring
echo "[STEP 5] Calculating empirical Z-scores..."
python scripts/05_null_model/compute_empirical_zscores.py

# Step 6: Pathway enrichment
echo "[STEP 6] Running GSEA..."
python scripts/06_enrichment/run_gsea.py

# Step 7: Cell-type specificity
echo "[STEP 7] Running cell-type footprint enrichment..."
python scripts/07_celltype_analysis/footprint_normalization.py

# Step 8: Sensitivity analysis
echo "[STEP 8] Running sensitivity analyses..."
python scripts/08_sensitivity_analysis/run_sensitivity.py

# Step 9a: Hi-C loop overlap
echo "[STEP 9a] Running Hi-C validation..."
python scripts/09_validation/hic_chromatin_loops.py

# Step 9b: GTEx eQTL validation (best effort + explicit stub if unavailable)
echo "[STEP 9b] Running GTEx validation..."
python scripts/09_validation/gtex_eqtl_validation.py

# Step 9c: SCHEMA convergence
echo "[STEP 9c] Running SCHEMA convergence analysis..."
python scripts/09c_convergence_schema.py

# Figures/tables
echo "[ASSETS] Generating manuscript tables and figures..."
python scripts/generate_figures.py

echo "======================================================================"
echo "PIPELINE COMPLETED"
echo "Outputs: data/processed, results/tables, results/figures"
echo "======================================================================"
