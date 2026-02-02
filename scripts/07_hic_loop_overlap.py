"""
Script 10: Direct Hi-C Loop Overlap Validation

Tests if SCZ GWAS variants physically overlap with Hi-C chromatin loops
from adult brain tissue.

Input: adultbrain_hic.bedpe (BEDPE format: chr1 start1 end1 chr2 start2 end2)
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
HIC_FILE = VALIDATION_DIR / "adultbrain_hic.bedpe"

def load_snps():
    """Load SCZ lead SNPs with coordinates."""
    snps_file = PROCESSED_DIR / "scz_lead_snps_robust.csv"
    if not snps_file.exists():
        return None
    return pd.read_csv(snps_file)

def load_hic_loops():
    """Load Hi-C loops from BEDPE file."""
    if not HIC_FILE.exists():
        print(f"Error: Hi-C file not found: {HIC_FILE}")
        return None
    
    # File is space-separated, not tab-separated
    df = pd.read_csv(HIC_FILE, sep=r'\s+')
    # Rename columns to standard BEDPE names if needed
    if 'chrom1' not in df.columns and len(df.columns) >= 6:
        df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'] + list(df.columns[6:])
    print(f"Loaded {len(df)} Hi-C loops")
    print(f"Columns: {list(df.columns)}")
    return df

def snp_in_loop(chrom, pos, loops_df):
    """Check if SNP falls within a loop anchor."""
    chrom_str = f'chr{chrom}' if not str(chrom).startswith('chr') else str(chrom)
    
    # Check anchor 1
    anchor1 = loops_df[
        (loops_df['chrom1'] == chrom_str) & 
        (loops_df['start1'] <= pos) & 
        (loops_df['end1'] >= pos)
    ]
    
    # Check anchor 2
    anchor2 = loops_df[
        (loops_df['chrom2'] == chrom_str) & 
        (loops_df['start2'] <= pos) & 
        (loops_df['end2'] >= pos)
    ]
    
    if len(anchor1) > 0 or len(anchor2) > 0:
        # Return the interacting coordinates
        interacting = []
        for _, row in anchor1.iterrows():
            interacting.append((row['chrom2'], row['start2'], row['end2']))
        for _, row in anchor2.iterrows():
            interacting.append((row['chrom1'], row['start1'], row['end1']))
        return True, interacting
    
    return False, []

def main():
    print("=" * 70)
    print("DIRECT Hi-C LOOP OVERLAP VALIDATION")
    print("=" * 70)
    
    # Load data
    snps_df = load_snps()
    if snps_df is None:
        return
    
    loops_df = load_hic_loops()
    if loops_df is None:
        return
    
    print(f"Testing {len(snps_df)} SCZ SNPs against {len(loops_df)} loops\n")
    
    # Test each SNP
    results = []
    loop_hits = 0
    
    for _, row in snps_df.iterrows():
        chrom = row['CHR']
        pos = int(row['BP'])
        snp = row['SNP']
        
        in_loop, interacting = snp_in_loop(chrom, pos, loops_df)
        
        if in_loop:
            loop_hits += 1
            results.append({
                'SNP': snp,
                'CHR': chrom,
                'BP': pos,
                'InLoop': True,
                'N_Interacting': len(interacting),
                'Interacting_Coords': ';'.join([f"{c}:{s}-{e}" for c, s, e in interacting[:3]])  # First 3
            })
        else:
            results.append({
                'SNP': snp,
                'CHR': chrom,
                'BP': pos,
                'InLoop': False,
                'N_Interacting': 0,
                'Interacting_Coords': ''
            })
    
    # Summary
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    
    rate = loop_hits / len(snps_df) * 100
    print(f"\nSNPs in Hi-C loop anchors: {loop_hits} / {len(snps_df)} ({rate:.1f}%)")
    
    # Expected by chance (genome coverage by loops)
    # Estimate loop anchor coverage
    total_anchor_bp = 0
    for _, row in loops_df.iterrows():
        total_anchor_bp += (row['end1'] - row['start1'])
        total_anchor_bp += (row['end2'] - row['start2'])
    
    genome_size = 3.1e9
    expected_rate = (total_anchor_bp / genome_size) * 100
    enrichment = rate / expected_rate if expected_rate > 0 else float('inf')
    
    print(f"Expected by chance: {expected_rate:.2f}%")
    print(f"Enrichment: {enrichment:.2f}x")
    
    # Save results
    results_df = pd.DataFrame(results)
    
    # Save hits only
    hits_df = results_df[results_df['InLoop'] == True]
    output_file = RESULTS_DIR / "hic_loop_overlap.csv"
    hits_df.to_csv(output_file, index=False)
    print(f"\n{len(hits_df)} SNPs in loops saved to {output_file}")
    
    return results_df

if __name__ == "__main__":
    main()
