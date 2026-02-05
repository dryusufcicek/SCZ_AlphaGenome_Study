#!/usr/bin/env python3
"""
03b_aggregate_scores_weighted.py
================================

Step 3b: Aggregating Variant Scores into Gene Scores.

Addressing Reviewer Critique #2: "Aggregation Bias: Winner's Curse in Gene Scoring...
The analysis summarizes gene scores by taking the maximum value...
Solution: PP-Weighted Scoring."

Methodology:
1. Load Variant-Level Scores (from Step 3).
2. Load Proxy Posterior Probabilities (from Step 1).
3. Compute Weighted Gene Score:
   GeneScore = Sum( Variant_Effect_on_Gene * Variant_PP )
   
   - If a gene is hit by one noisy variant (high score, low PP) -> Score is suppressed.
   - If a gene is hit by the causal variant (high score, high PP) -> Score is retained.
   - If a gene is hit by multiple likely variants -> Score accumulates.

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR

def parse_gene_scores(row):
    """
    Parse the 'all_genes' column string into a dictionary of {gene: score}.
    Format: "GENE1:0.5;GENE2:0.9"
    """
    gene_map = {}
    val = row.get('all_genes')
    
    if pd.isna(val) or not val:
        # Fallback: if 'target_gene' is set but all_genes missing
        tgt = row.get('target_gene')
        score = row.get('target_gene_score', 0)
        if pd.notna(tgt) and score > 0:
            gene_map[tgt] = score
        return gene_map
        
    for entry in str(val).split(';'):
        if ':' in entry:
            parts = entry.split(':')
            g = parts[0]
            try:
                s = float(parts[1])
                gene_map[g] = s
            except:
                pass
    return gene_map

def aggregate_scores(scores_file):
    print(f"Loading scores from {scores_file}...")
    
    # Read CSV in chunks if huge, but 10k variants fits in memory
    try:
        df = pd.read_csv(scores_file)
    except:
        print("❌ Could not read scores file.")
        return None
        
    print(f"Loaded {len(df)} variant scores.")
    
    # Check for Proxy_PP column
    if 'Proxy_PP' not in df.columns:
        print("⚠️ 'Proxy_PP' column missing. Merging from extraction file...")
        # Load Step 1 file
        meta_file = PROCESSED_DIR / "scz_credible_sets_proxy.csv"
        if meta_file.exists():
            meta = pd.read_csv(meta_file, usecols=['SNP', 'Proxy_PP', 'Locus_ID'])
            # Merge on SNP (and Locus if needed)
            df = df.merge(meta, on='SNP', how='left')
            print(f"Merged PP data. variants with PP: {df['Proxy_PP'].notna().sum()}")
        else:
            print("❌ Critical: No PP data available. Cannot perform weighted scoring.")
            return None
            
    # Fill missing PP with 0 (safe fallback) or equal weight
    df['Proxy_PP'] = df['Proxy_PP'].fillna(0.0)
    
    # Aggregate
    gene_weighted_sums = {}
    gene_max_scores = {} # Keep max for comparison
    
    print("Aggregating Gene Scores (PP-Weighted)...")
    for _, row in tqdm(df.iterrows(), total=len(df)):
        pp = row['Proxy_PP']
        locus = row.get('Locus_ID', -1)
        
        # Get all genes affected by this variant
        Variant_Genes = parse_gene_scores(row)
        
        for gene, effect in Variant_Genes.items():
            # Weighted Sum: Score * Probability
            contribution = effect * pp
            
            if gene not in gene_weighted_sums:
                gene_weighted_sums[gene] = 0.0
                gene_max_scores[gene] = 0.0
                
            gene_weighted_sums[gene] += contribution
            gene_max_scores[gene] = max(gene_max_scores[gene], effect)
            
    # Format Results
    results = []
    for gene in gene_weighted_sums:
        results.append({
            'gene': gene,
            'weighted_score': gene_weighted_sums[gene],
            'max_score': gene_max_scores[gene],
            # 'n_variants': ... could add this
        })
        
    res_df = pd.DataFrame(results).sort_values('weighted_score', ascending=False)
    
    print("\nTop 10 Genes (Weighted):")
    print(res_df.head(10).to_string(index=False))
    
    return res_df

def main():
    print("=" * 70)
    print("STEP 4: ROBUST SCORE AGGREGATION")
    print("=" * 70)
    
    # Input: Result of Step 3 (which is running)
    scores_file = PROCESSED_DIR / "alphgenome_variant_scores_proxy.csv"
    
    if not scores_file.exists():
        print(f"⚠️ Input file {scores_file} not ready yet.")
        print("Please wait for Step 2 (Scoring) to complete.")
        return

    res_df = aggregate_scores(scores_file)
    
    if res_df is not None:
        out_file = PROCESSED_DIR / "gene_scores_weighted.csv"
        res_df.to_csv(out_file, index=False)
        print(f"\n✅ Saved weighted scores to {out_file}")
        print("Winner's Curse corrected.")

if __name__ == "__main__":
    main()
