#!/usr/bin/env python3
"""
01_extract_gwas_variants.py
===========================

Step 1: Extract lead SNPs from PGC3 Schizophrenia GWAS summary statistics.

This script reads the compressed GWAS summary statistics and extracts
genome-wide significant (P < 5e-8) variants for downstream AlphaGenome scoring.

Input:
    - PGC3 GWAS summary statistics (daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz)

Output:
    - data/processed/scz_lead_snps.csv

Usage:
    python scripts/01_extract_gwas_variants.py

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import sys

# Add parent directory for config import
sys.path.append(str(Path(__file__).parent.parent))
from config import RAW_DIR, PROCESSED_DIR, GWAS_PVALUE_THRESHOLD


def load_gwas_summary_stats(filepath):
    """
    Load PGC3 GWAS summary statistics from gzipped file.
    
    Parameters
    ----------
    filepath : Path
        Path to gzipped summary statistics file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: SNP, CHR, BP, A1, A2, OR, P
    """
    print(f"Loading GWAS summary statistics from {filepath}...")
    
    df = pd.read_csv(
        filepath,
        sep='\t',
        compression='gzip',
        dtype={'CHR': str, 'BP': int, 'SNP': str}
    )
    
    print(f"  Loaded {len(df):,} variants")
    return df


def extract_significant_snps(df, p_threshold=GWAS_PVALUE_THRESHOLD):
    """
    Extract genome-wide significant SNPs.
    
    Parameters
    ----------
    df : pd.DataFrame
        Full GWAS summary statistics
    p_threshold : float
        P-value threshold (default: 5e-8)
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with significant SNPs
    """
    significant = df[df['P'] < p_threshold].copy()
    
    print(f"  {len(significant):,} variants below P < {p_threshold:.0e}")
    
    # Sort by chromosome and position
    significant = significant.sort_values(['CHR', 'BP'])
    
    return significant


def identify_lead_snps(df, window_kb=500, r2_threshold=0.1):
    """
    Identify LD-independent lead SNPs using distance-based pruning.
    
    Note: This is an approximation. For precise LD pruning, 
    use PLINK with reference panel.
    
    Parameters
    ----------
    df : pd.DataFrame
        Significant SNPs
    window_kb : int
        Window size in kilobases for pruning
        
    Returns
    -------
    pd.DataFrame
        Lead SNPs only
    """
    print(f"  Identifying lead SNPs (window = {window_kb}kb)...")
    
    leads = []
    
    for chrom in df['CHR'].unique():
        chr_df = df[df['CHR'] == chrom].sort_values('P')
        
        selected_positions = []
        
        for _, row in chr_df.iterrows():
            pos = row['BP']
            
            # Check if too close to already selected SNP
            is_independent = True
            for sel_pos in selected_positions:
                if abs(pos - sel_pos) < window_kb * 1000:
                    is_independent = False
                    break
            
            if is_independent:
                leads.append(row)
                selected_positions.append(pos)
    
    lead_df = pd.DataFrame(leads)
    print(f"  Identified {len(lead_df):,} lead SNPs")
    
    return lead_df


def main():
    """Main extraction pipeline."""
    print("=" * 70)
    print("STEP 1: GWAS VARIANT EXTRACTION")
    print("=" * 70)
    
    # Define input file
    gwas_file = RAW_DIR / "daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz"
    
    if not gwas_file.exists():
        print(f"\nError: GWAS file not found at {gwas_file}")
        print("Please download from PGC data portal and place in data/raw/")
        return None
    
    # Load and process
    df = load_gwas_summary_stats(gwas_file)
    significant = extract_significant_snps(df)
    leads = identify_lead_snps(significant)
    
    # Save output
    output_file = PROCESSED_DIR / "scz_lead_snps.csv"
    leads.to_csv(output_file, index=False)
    
    print(f"\n✓ Saved {len(leads)} lead SNPs to {output_file}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total variants in GWAS: {len(df):,}")
    print(f"  Genome-wide significant: {len(significant):,}")
    print(f"  LD-independent leads: {len(leads):,}")
    
    return leads


if __name__ == "__main__":
    main()
