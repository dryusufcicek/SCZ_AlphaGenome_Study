"""
Script 08: Cell-Type-Specific Peak Enrichment

Tests if SCZ GWAS variants are enriched in open chromatin peaks
from specific iPSC-derived neuronal subtypes.

Cell Types:
- iN_Glut: Glutamatergic neurons
- iN_GABA: GABAergic interneurons  
- iN_Dopa: Dopaminergic neurons
- NPC: Neural Progenitor Cells
- iPSC: Induced Pluripotent Stem Cells (negative control)
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

# Cell type peak files
PEAK_FILES = {
    'Glutamatergic': 'iN_Glut.txt',
    'GABAergic': 'iN_GABA.txt',
    'Dopaminergic': 'iN_Dopa.txt',
    'NPC': 'NPC.txt',
    'iPSC': 'iPSC.txt'
}

def load_variants():
    """Load SCZ lead SNPs with coordinates."""
    snps_file = PROCESSED_DIR / "scz_lead_snps_robust.csv"
    if not snps_file.exists():
        print(f"Error: SNPs file not found: {snps_file}")
        return None
    
    df = pd.read_csv(snps_file)
    # Ensure CHR is formatted correctly
    df['CHR_clean'] = df['CHR'].astype(str).str.replace('chr', '')
    return df

def load_peaks(peak_file):
    """Load ATAC-seq peaks from BED-like file."""
    filepath = VALIDATION_DIR / peak_file
    if not filepath.exists():
        print(f"Warning: Peak file not found: {filepath}")
        return None
    
    df = pd.read_csv(filepath, sep='\t')
    # Standardize CHR
    df['CHR_clean'] = df['CHR'].astype(str).str.replace('chr', '')
    return df

def count_snps_in_peaks(snps_df, peaks_df):
    """Count how many SNPs fall within peaks."""
    if peaks_df is None or snps_df is None:
        return 0
    
    hits = 0
    
    # Group peaks by chromosome for efficiency
    peaks_by_chr = {}
    for chrom in peaks_df['CHR_clean'].unique():
        chr_peaks = peaks_df[peaks_df['CHR_clean'] == chrom]
        peaks_by_chr[chrom] = list(zip(chr_peaks['START'], chr_peaks['END']))
    
    # Check each SNP
    for _, row in snps_df.iterrows():
        chrom = str(row['CHR_clean'])
        pos = int(row['BP'])
        
        if chrom in peaks_by_chr:
            for start, end in peaks_by_chr[chrom]:
                if start <= pos <= end:
                    hits += 1
                    break
    
    return hits

def main():
    print("=" * 70)
    print("CELL-TYPE-SPECIFIC PEAK ENRICHMENT")
    print("=" * 70)
    
    # Load SCZ variants
    snps_df = load_variants()
    if snps_df is None:
        return
    
    n_snps = len(snps_df)
    print(f"Loaded {n_snps} SCZ lead SNPs\n")
    
    # Test each cell type
    results = []
    
    print(f"{'Cell Type':<20} {'Peaks':<12} {'SNPs in Peaks':<15} {'Rate':<10}")
    print("-" * 60)
    
    for cell_type, peak_file in PEAK_FILES.items():
        peaks_df = load_peaks(peak_file)
        
        if peaks_df is None:
            continue
        
        n_peaks = len(peaks_df)
        hits = count_snps_in_peaks(snps_df, peaks_df)
        rate = hits / n_snps if n_snps > 0 else 0
        
        print(f"{cell_type:<20} {n_peaks:<12} {hits:<15} {rate:.1%}")
        
        results.append({
            'CellType': cell_type,
            'N_Peaks': n_peaks,
            'SNPs_in_Peaks': hits,
            'Rate': rate
        })
    
    # Statistical comparison (GABAergic vs iPSC as control)
    print("\n" + "=" * 70)
    print("ENRICHMENT ANALYSIS (vs iPSC control)")
    print("=" * 70)
    
    results_df = pd.DataFrame(results)
    
    if 'iPSC' in results_df['CellType'].values:
        ipsc_rate = results_df[results_df['CellType'] == 'iPSC']['Rate'].values[0]
        ipsc_hits = results_df[results_df['CellType'] == 'iPSC']['SNPs_in_Peaks'].values[0]
        
        for _, row in results_df.iterrows():
            if row['CellType'] == 'iPSC':
                continue
            
            # Fisher's exact test: is this cell type enriched vs iPSC?
            # Contingency table:
            #           In_Peak  Not_In_Peak
            # CellType    a         b
            # iPSC        c         d
            
            a = row['SNPs_in_Peaks']
            b = n_snps - a
            c = ipsc_hits
            d = n_snps - c
            
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
            
            sig = "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
            direction = "ENRICHED" if odds_ratio > 1 else "depleted"
            
            print(f"{row['CellType']:<20} OR={odds_ratio:.2f}  P={p_value:.3f} {sig} ({direction})")
            
            # Update results
            results_df.loc[results_df['CellType'] == row['CellType'], 'OR_vs_iPSC'] = odds_ratio
            results_df.loc[results_df['CellType'] == row['CellType'], 'P_vs_iPSC'] = p_value
    
    # Save results
    output_file = RESULTS_DIR / "celltype_peak_enrichment.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    
    return results_df

if __name__ == "__main__":
    main()
