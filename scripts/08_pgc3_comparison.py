"""
Script 11: AlphaGenome vs PGC3 Prioritized Gene Comparison

Compares AlphaGenome predicted target genes with the official
PGC3 prioritized gene list from Trubetskoy et al. 2022.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
import sys

# Setup
BASE_DIR = Path(".").resolve()
if "AlphaGenome with SCZ" not in str(BASE_DIR):
    BASE_DIR = Path("/Users/yusuf/AlphaGenome with SCZ")

VALIDATION_DIR = BASE_DIR / "data" / "validation"
PROCESSED_DIR = BASE_DIR / "scz_hypothesis_testing/data/processed"
RESULTS_DIR = BASE_DIR / "scz_hypothesis_testing/results"
PGC3_FILE = VALIDATION_DIR / "scz2022-Extended-Data-Table1.xlsx"

def load_pgc3_genes():
    """Load official PGC3 prioritized genes."""
    if not PGC3_FILE.exists():
        return None
    
    xlsx = pd.ExcelFile(PGC3_FILE)
    
    # Load Extended Data Table 1 (120 prioritized genes)
    df = pd.read_excel(xlsx, sheet_name='Extended.Data.Table.1', header=None, skiprows=1)
    # Columns: Index.SNP, Ensembl.ID, Symbol.ID, gene_biotype, ...
    df.columns = ['SNP', 'Ensembl', 'Symbol', 'Biotype'] + [f'col{i}' for i in range(4, len(df.columns))]
    
    genes = set(df['Symbol'].dropna().str.upper())
    snps = set(df['SNP'].dropna())
    
    print(f"PGC3 Prioritized: {len(genes)} genes, {len(snps)} loci")
    return genes, snps

def load_alphagenome_genes():
    """Load AlphaGenome predicted target genes."""
    scores_file = PROCESSED_DIR / "robust_modality_scores.csv"
    if not scores_file.exists():
        return None
    
    df = pd.read_csv(scores_file)
    
    genes = set()
    for genes_str in df['genes'].dropna():
        for entry in str(genes_str).split(';'):
            if ':' in entry:
                gene, score = entry.split(':')
                if float(score) >= 1.0:  # Only high-scoring genes
                    genes.add(gene.split('.')[0].upper())
    
    print(f"AlphaGenome Targets (score >= 1.0): {len(genes)} genes")
    return genes

def main():
    print("=" * 70)
    print("ALPHAGENOME vs PGC3 PRIORITIZED GENE COMPARISON")
    print("=" * 70)
    
    # Load gene lists
    pgc3_genes, pgc3_snps = load_pgc3_genes()
    ag_genes = load_alphagenome_genes()
    
    if pgc3_genes is None or ag_genes is None:
        print("Error loading gene lists")
        return
    
    # Calculate overlap
    overlap = pgc3_genes.intersection(ag_genes)
    
    print("\n" + "=" * 70)
    print("OVERLAP ANALYSIS")
    print("=" * 70)
    
    print(f"\nOverlapping genes: {len(overlap)} / {len(pgc3_genes)} ({len(overlap)/len(pgc3_genes)*100:.1f}%)")
    print(f"\nGenes in BOTH lists:")
    for gene in sorted(overlap):
        print(f"  - {gene}")
    
    # Genes in PGC3 but not AlphaGenome
    pgc3_only = pgc3_genes - ag_genes
    print(f"\nPGC3 only ({len(pgc3_only)} genes): {', '.join(sorted(list(pgc3_only)[:10]))}...")
    
    # Genes in AlphaGenome but not PGC3
    ag_only = ag_genes - pgc3_genes
    print(f"\nAlphaGenome only ({len(ag_only)} genes): {', '.join(sorted(list(ag_only)[:10]))}...")
    
    # Fisher's exact test
    # Estimate genome background (20,000 protein-coding genes)
    genome_size = 20000
    
    # Contingency table
    #             In AG    Not in AG
    # In PGC3       a         b
    # Not PGC3      c         d
    
    a = len(overlap)
    b = len(pgc3_only)
    c = len(ag_only)
    d = genome_size - a - b - c
    
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
    
    print("\n" + "=" * 70)
    print("STATISTICAL ENRICHMENT")
    print("=" * 70)
    print(f"\nFisher's Exact Test:")
    print(f"  Odds Ratio: {odds_ratio:.2f}")
    print(f"  P-value: {p_value:.2e}")
    
    if p_value < 0.05:
        print("  ✅ Significant overlap (P < 0.05)")
    else:
        print("  → Not significantly enriched")
    
    # Save results
    results = {
        'PGC3_Genes': len(pgc3_genes),
        'AlphaGenome_Genes': len(ag_genes),
        'Overlap': len(overlap),
        'Overlap_Percent': len(overlap)/len(pgc3_genes)*100,
        'OR': odds_ratio,
        'P': p_value
    }
    
    output_file = RESULTS_DIR / "pgc3_comparison.csv"
    pd.DataFrame([results]).to_csv(output_file, index=False)
    
    # Save overlap gene list
    overlap_df = pd.DataFrame({'Gene': sorted(list(overlap))})
    overlap_df.to_csv(RESULTS_DIR / "pgc3_overlap_genes.csv", index=False)
    
    print(f"\nResults saved to {output_file}")
    
    return results

if __name__ == "__main__":
    main()
