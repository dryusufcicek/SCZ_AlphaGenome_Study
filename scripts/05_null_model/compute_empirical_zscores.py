#!/usr/bin/env python3
"""
03c_empirical_z_scores.py
=========================

Step 2 of Critique Roadmap: Empirical Calibration of AlphaGenome Scores.

Problem: "AlphaGenome scores ... are model-internal units ... not calibrated."
Solution: "Introduce empirical nulls and Z-scoring."

Methodology:
1. Load Weighted Gene Scores (from Step 3b).
2. Load the "Background Universe" (from Step 4).
3. Construct the Null Distribution:
   - The scores of the ~40,000 non-target genes (Score = 0 or near 0 noise).
   - The scores of the ~2,000 target genes.
4. Calculate Z-Scores:
   - Mean/SD derived from the Whole Genome distribution.
   - Z = (Score - Mean_genome) / SD_genome.
5. Calculate Empirical P-values:
   - Rank-based P-value matching the Z-score.

Output:
   - `gene_z_scores.csv` with columns: [gene, raw_score, z_score, p_empirical]

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, FIGURES_DIR

def calculate_empirical_stats(weighted_scores_file, universe_file):
    print(f"Loading weighted scores from {weighted_scores_file}...")
    if not weighted_scores_file.exists():
        print("⚠️ Weighted scores file not found. Please run Step 3b first.")
        return None

    scores_df = pd.read_csv(weighted_scores_file)
    print(f"Loaded {len(scores_df)} gene scores.")

    # Load Universe (using JSI file as proxy for 'all genes')
    uni_df = pd.read_csv(universe_file)
    # Find gene column
    col = next((c for c in uni_df.columns if 'gene' in c.lower() or 'symbol' in c.lower()), None)
    all_genes = set(uni_df[col].dropna().unique())
    print(f"Loaded Universe: {len(all_genes)} genes.")

    # Merge: Create Full Genome DataFrame
    # Genes present in scores_df keep their score
    # Genes NOT present get Score = 0.0 (Background noise assumption)
    
    score_map = dict(zip(scores_df['gene'], scores_df['weighted_score']))
    
    full_data = []
    for gene in all_genes:
        s = score_map.get(gene, 0.0)
        full_data.append({'gene': gene, 'raw_score': s})
        
    # Add any genes in scores but not in universe (rare but possible)
    for gene, s in score_map.items():
        if gene not in all_genes:
            full_data.append({'gene': gene, 'raw_score': s})
            
    df = pd.DataFrame(full_data)
    print(f"Full Genome Distribution: {len(df)} genes.")
    
    # Calculate Statistics
    # We use the whole genome as the null (conservative approach for "Global Specificity")
    # Alternatively, we could sample random SNPs, but "Whole Genome" is the standard background for GSEA.
    
    mu = df['raw_score'].mean()
    sigma = df['raw_score'].std()
    
    print(f"\nGenome-wide Statistics:")
    print(f"  Mean: {mu:.4f}")
    print(f"  Std : {sigma:.4f}")
    
    # Calculate Z-Scores
    df['z_score'] = (df['raw_score'] - mu) / sigma
    
    # Calculate Empirical P-values (Rank / N)
    # Sort descending
    df = df.sort_values('raw_score', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    df['p_empirical'] = df['rank'] / len(df)
    
    # Filter to only genes with non-zero scores for the output (to save space, but keeping Z-context)
    # Actually, let's keep top 500 or all non-zeros
    output_df = df[df['raw_score'] > 0].copy()
    
    print(f"\nTop 5 Genes by Z-Score:")
    print(output_df[['gene', 'raw_score', 'z_score', 'p_empirical']].head(5))
    
    return output_df, df['raw_score']

def plot_distribution(all_scores, output_path):
    """Plot the distribution of scores with Z-score cutoffs."""
    plt.figure(figsize=(10, 6))
    
    # Log scale often helps if scores are power-law distributed
    # But Z-scores imply normal-ish. Let's inspect raw first.
    # Usually regulatory scores are heavy-tailed.
    
    sns.histplot(all_scores, bins=100, log_scale=(False, True))
    plt.title("Distribution of Genome-Wide Regulatory Burden (AlphaGenome)")
    plt.xlabel("Weighted Score")
    plt.ylabel("Count (Log Scale)")
    
    plt.savefig(output_path)
    print(f"Saved distribution plot to {output_path}")

def main():
    print("=" * 70)
    print("STEP 2 (NEW): EMPIRICAL Z-SCORING")
    print("=" * 70)
    
    scores_file = PROCESSED_DIR / "gene_scores_weighted.csv"
    universe_file = PROCESSED_DIR / "gene_jsi.csv"
    
    if not scores_file.exists():
        print(f"Scores file {scores_file} missing. Run Step 3b first.")
        # Create dummy for verification if needed
        # return 
        
    result_df, all_scores = calculate_empirical_stats(scores_file, universe_file)
    
    if result_df is not None:
        out_file = PROCESSED_DIR / "gene_z_scores.csv"
        result_df.to_csv(out_file, index=False)
        print(f"\n✅ Saved Z-scores to {out_file}")
        
        # Plot
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        plot_distribution(all_scores, FIGURES_DIR / "score_distribution_z.png")

if __name__ == "__main__":
    main()
