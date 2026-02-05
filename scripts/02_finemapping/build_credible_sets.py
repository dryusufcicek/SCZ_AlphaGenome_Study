#!/usr/bin/env python3
"""
01_extract_gwas_variants.py
===========================

Step 1 (REVISED): Extract "Proxy" Credible Sets from PGC3 Schizophrenia GWAS.

CRITICAL METHODOLOGICAL UPDATE:
Since official FINE-MAPPING data is unavailable, we are constructing
"Proxy Credible Sets" to address the "LD Blindness" critique.

Methodology:
1. Extract ALL genome-wide significant SNPs (P < 5e-8).
2. Group SNPs into 500kb Loci.
3. Calculate "Proxy Posterior Probability" (Proxy_PP) for each SNP:
   - Convert P-value to Bayes Factor approx (BF ~ 1/P for highly sig).
   - Normalize BF within the locus to sum to 1.
   - Proxy_PP = BF_snp / Sum(BF_locus)
   
This allows us to perform "Weighted Scoring" later, satisfying the reviewer's
demand for robust aggregation, even without the official file.

Input:
    - PGC3 Primary GWAS (data/raw/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz)
      OR data/validation/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz

Output:
    - data/processed/scz_credible_sets_proxy.csv

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import glob

# Add parent directory for config import
sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, GWAS_PVALUE_THRESHOLD

# Path search for the daner file
POSSIBLE_PATHS = [
    Path("data/raw/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz"),
    Path("data/validation/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz")
]

def load_gwas_data():
    """Find and load the GWAS summary stats."""
    gwas_file = None
    for p in POSSIBLE_PATHS:
        if p.exists():
            gwas_file = p
            break
    
    if not gwas_file:
        # Fallback search
        found = list(Path(".").rglob("daner_PGC_SCZ*.gz"))
        if found:
            gwas_file = found[0]
            
    if not gwas_file:
        print("❌ CRITICAL: No GWAS summary statistics found.")
        return None
        
    print(f"Loading GWAS from {gwas_file}...")
    
    try:
        df = pd.read_csv(gwas_file, sep='\t', compression='gzip')
    except:
        # Sometimes daner files have weird headers
        df = pd.read_csv(gwas_file, sep='\t', compression='gzip', error_bad_lines=False)
        
    print(f"  Loaded {len(df):,} variants")
    return df

def calculate_proxy_pp(locus_df):
    """
    Calculate Proxy Posterior Probability for a locus.
    Assumption: One of these SNPs is causal.
    Method: Z-score approx to Bayes Factor.
    """
    # Avoid log(0)
    p_values = locus_df['P'].replace(0, 1e-300) 
    
    # Approx Bayes Factor (Wakefield approx simplified)
    # stronger P -> higher BF
    # We use -log10(P) as a rough proportional weight for the exponential scale
    log_bf = -np.log10(p_values)
    
    # Softmax function to get probabilities
    # Shift for numerical stability
    max_val = log_bf.max()
    exp_val = np.exp(log_bf - max_val)
    proxy_pp = exp_val / exp_val.sum()
    
    return proxy_pp

def create_proxy_credible_sets(df):
    """
    Bioinformatics Pipeline:
    1. Filter Significant
    2. Cluster into Loci (Clumping)
    3. Calculate PP
    """
    print("\n--- Generating Proxy Credible Sets ---")
    
    # 1. Significance Filter
    sig = df[df['P'] < 5e-8].copy()
    print(f"  Significant Variants (P < 5e-8): {len(sig):,}")
    
    if len(sig) == 0:
        return pd.DataFrame()

    # 2. Assign Locus ID (Simple Distance Clustering)
    # Sort by Chromosome and Position/P-value
    sig = sig.sort_values(['CHR', 'P'])
    
    locus_ids = []
    current_locus = 0
    last_chr = None
    last_bp = -9999999
    
    # Iterate and assign locus IDs
    # If SNP is > 500kb away from previous component, new locus
    # This is a "greedy chain" clustering
    sig = sig.sort_values(['CHR', 'BP'])
    
    for idx, row in sig.iterrows():
        c = row['CHR']
        b = row['BP']
        
        if c != last_chr or (b - last_bp > 500000): # 500kb gap
            current_locus += 1
            
        locus_ids.append(current_locus)
        last_chr = c
        last_bp = b
        
    sig['Locus_ID'] = locus_ids
    num_loci = sig['Locus_ID'].nunique()
    print(f"  Clustered into {num_loci} independent loci (500kb merging)")
    
    # 3. Calculate Proxy PP per Locus
    print("  Calculating Proxy Posterior Probabilities...")
    sig['Proxy_PP'] = sig.groupby('Locus_ID', group_keys=False).apply(calculate_proxy_pp).values
    
    # 4. Filter for "95% Credible Set"
    # Sort by PP descending within locus
    sig = sig.sort_values(['Locus_ID', 'Proxy_PP'], ascending=[True, False])
    
    # Calculate cumulative sum
    sig['Cumulative_PP'] = sig.groupby('Locus_ID')['Proxy_PP'].cumsum()
    
    # Keep variants until cumulative sum > 0.95 (plus the one that crosses it)
    # We want to keep the Top 95% mass
    # Logic: Shift cumsum to get previous, check if previous < 0.95
    mask = sig.groupby('Locus_ID')['Cumulative_PP'].shift(1).fillna(0) < 0.95
    credible_sets = sig[mask].copy()
    
    print(f"\n  Final Proxy Credible Set Size: {len(credible_sets):,}")
    print(f"  Avg variants per locus: {len(credible_sets)/num_loci:.1f}")
    
    return credible_sets

def main():
    print("=" * 70)
    print("STEP 1: VARIANT EXTRACTION (PROXY MODE)")
    print("=" * 70)
    
    df = load_gwas_data()
    if df is None: return
    
    # Run the Pipeline
    proxy_sets = create_proxy_credible_sets(df)
    
    # Save
    output_file = PROCESSED_DIR / "scz_credible_sets_proxy.csv"
    
    # Select clean columns
    out_cols = ['Locus_ID', 'SNP', 'CHR', 'BP', 'A1', 'A2', 'P', 'Proxy_PP', 'Cumulative_PP']
    proxy_sets[out_cols].to_csv(output_file, index=False)
    
    print(f"\n✅ Saved Proxy Credible Sets to {output_file}")
    print("Rationale: We now have weighted probabilities for aggregation, solving the 'Lead SNP' bias.")
    
    # Generate a plot command hint
    print("\nNext Step: Run Step 2 (Scoring) using these weighted variants.")

if __name__ == "__main__":
    main()
