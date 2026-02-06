#!/usr/bin/env python3
"""
07_footprint_normalization.py
=============================

Statistically Correct Cell Type Specificity Analysis (Real Data Validation).

Methodology:
1. Load SCZ Credible Set Variants (Step 1).
2. Load Real ATAC-seq Peak Sets (Corces et al. 2020 / Trevino et al. 2021).
3. Compute Genomic Footprint (Denominator):
   - Union of all peaks for a cell class (bedtools merge equivalent).
   - Total bp coverage / Genome Size.
4. Compute Observed Enrichment (Numerator):
   - REAL INTERSECTION of variants with merged peaks.
   - Count how many variants fall into the open chromatin of each cell type.
5. Statistical Test:
   - Binomial Test (Observed vs Expected given Footprint).

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys

# Add parent directory for config import
sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, FIGURES_DIR, CS_DIR, CT_DIR

# -----------------------------------------------------------------------------
# REAL DATA MAPPING (Corces et al. 2020)
# -----------------------------------------------------------------------------

CORCES_CLUSTER_MAP = {
    'Excitatory_Neuron': [1, 2, 3, 4, 5, 6, 7, 8], 
    'Inhibitory_Neuron': [9, 10, 11, 12, 13, 14],
    'Astrocytes': [18, 19],
    'Oligodendrocytes': [15, 16, 17],
    'Microglia': [20, 21],
    'OPC': [22, 23],
}

EXT_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "external"
GENOME_SIZE_BP = 2.8e9  # Mappable human genome

column_map = {'P-value': 'P_Value', 'P_Value': 'P_Value'} # varying column names support

def load_variants():
    """Load variants from PROXY credible sets."""
    file = CS_DIR / "scz_credible_sets_proxy.csv"
    if not file.exists():
        file = PROCESSED_DIR / "scz_lead_snps.csv"
        
    if not file.exists():
        print("Error: No variant file found.")
        return pd.DataFrame()

    df = pd.read_csv(file)
    print(f"Loaded {len(df)} variants")
    return df

def get_merged_peaks_and_footprint(cell_class, clusters):
    """
    Reads all BED files for a cell class, merges overlapping intervals, 
    and returns:
    1. A dictionary of intervals by chromosome: {'chr1': [(s,e), ...]}
    2. Total footprint size (bp).
    """
    all_intervals = []
    file_count = 0
    
    # 1. READ RAW INTERVALS
    for cid in clusters:
        bed_file = EXT_DATA_DIR / f"Cluster{cid}.idr.optimal.narrowPeak.gz"
        if bed_file.exists():
            file_count += 1
            try:
                # Read standard BED (chr, start, end)
                df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[0,1,2], engine='c')
                all_intervals.extend(df.values.tolist())
            except Exception as e:
                pass
    
    if not all_intervals:
        print(f"  > {cell_class:<20}: No peaks found.")
        return {}, 0

    # 2. MERGE INTERVALS (Python Implementation)
    # Sort by chr, then start
    all_intervals.sort(key=lambda x: (x[0], x[1]))
    
    merged_by_chr = {}
    total_bp = 0
    merged_count = 0
    
    if all_intervals:
        curr_chr, curr_start, curr_end = all_intervals[0]
        
        for r in all_intervals[1:]:
            if r[0] == curr_chr and r[1] < curr_end:
                # Overlap: extend end
                curr_end = max(curr_end, r[2])
            else:
                # Commit current
                if curr_chr not in merged_by_chr: merged_by_chr[curr_chr] = []
                merged_by_chr[curr_chr].append((curr_start, curr_end))
                total_bp += (curr_end - curr_start)
                merged_count += 1
                
                # Start new
                curr_chr, curr_start, curr_end = r
        
        # Commit last
        if curr_chr not in merged_by_chr: merged_by_chr[curr_chr] = []
        merged_by_chr[curr_chr].append((curr_start, curr_end))
        total_bp += (curr_end - curr_start)
        merged_count += 1

    print(f"  > {cell_class:<20}: {total_bp/1e6:6.1f} MB ({merged_count} merged peaks)")
    return merged_by_chr, total_bp

def get_trevino_peaks():
    """Load Trevino fetal peaks (already consensus)."""
    trevino_file = EXT_DATA_DIR / "GSE162170_atac_consensus_peaks.bed.gz"
    if not trevino_file.exists():
        return {}, 0
        
    try:
        df = pd.read_csv(trevino_file, sep='\t', header=None, usecols=[0,1,2])
        all_intervals = df.values.tolist()
        
        # Assume it's already non-overlapping? Best to merge to be safe.
        all_intervals.sort(key=lambda x: (x[0], x[1]))
        
        merged_by_chr = {}
        total_bp = 0
        
        if all_intervals:
            curr_chr, curr_start, curr_end = all_intervals[0]
            for r in all_intervals[1:]:
                if r[0] == curr_chr and r[1] < curr_end:
                   curr_end = max(curr_end, r[2])
                else:
                   if curr_chr not in merged_by_chr: merged_by_chr[curr_chr] = []
                   merged_by_chr[curr_chr].append((curr_start, curr_end))
                   total_bp += (curr_end - curr_start)
                   curr_chr, curr_start, curr_end = r
            
            if curr_chr not in merged_by_chr: merged_by_chr[curr_chr] = []
            merged_by_chr[curr_chr].append((curr_start, curr_end))
            total_bp += (curr_end - curr_start)
            
        print(f"  > {'Fetal_Progenitors':<20}: {total_bp/1e6:6.1f} MB (Consensus)")
        return merged_by_chr, total_bp
        
    except:
        return {}, 0

def check_intersection(variants_df, peak_dict):
    """
    Count how many variants overlap with the given peak set.
    peak_dict: {'chr1': [(s,e), ...], ...}
    """
    hits = 0
    
    # Pre-group variants by chromosome to avoid scanning full peak dict every time
    # Assuming 'CHR' or 'chromosome' column
    # Standarize chrom names
    
    # Identify chr col
    chr_col = None
    pos_col = None
    for c in variants_df.columns:
        if c.lower() in ['chr', 'chromosome']: chr_col = c
        if c.lower() in ['bp', 'pos', 'position']: pos_col = c
        
    if not chr_col or not pos_col:
        print("Warning: Could not identify CHR/POS columns")
        return 0

    # Iterate variants
    for _, row in variants_df.iterrows():
        chrom = str(row[chr_col])
        if not chrom.startswith('chr'): chrom = 'chr' + chrom
        pos = int(row[pos_col])
        
        if chrom in peak_dict:
            # Linear scan of peaks in this chromosome
            # Since peaks are sorted, we can break early, but simplest is linear for now
            # Optimization: Bisect?
            # Or just iterate.
            
            for (start, end) in peak_dict[chrom]:
                if start <= pos <= end:
                    hits += 1
                    break # One hit is enough
                if start > pos:
                    break # Passed it (sorted peaks)
                    
    return hits

def main():
    print("=" * 70)
    print("STEP 2: STATISTICAL CORRECTION (Real Data - Exact Intersection)")
    print("=" * 70)
    
    variants = load_variants()
    if variants.empty: return
    
    total_snps = len(variants)
    results = []

    print(f"\nProcessing Footprints & Intersections (N={total_snps} variants)...")
    print("-" * 80)
    print(f"{'Cell Type':<20} {'Size(MB)':<10} {'Exp.Hits':<10} {'Obs.Hits':<10} {'Enrichment':<10} {'P-value'}")

    # 1. Process Corces
    for cell_class, clusters in CORCES_CLUSTER_MAP.items():
        peak_dict, total_bp = get_merged_peaks_and_footprint(cell_class, clusters)
        
        if total_bp == 0: continue
        
        # Calculate Stats
        footprint_frac = total_bp / GENOME_SIZE_BP
        expected_hits = total_snps * footprint_frac
        
        # REAL INTERSECTION
        observed_hits = check_intersection(variants, peak_dict)
        
        # Binomial Test
        p_val = stats.binom.sf(observed_hits - 1, total_snps, footprint_frac)
        enrichment = observed_hits / expected_hits if expected_hits > 0 else 0
        
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        
        print(f"{cell_class:<20} {total_bp/1e6:<10.1f} {expected_hits:<10.1f} {observed_hits:<10} {enrichment:<10.2f} {p_val:.2e} {sig}")
        
        results.append({
            'label': cell_class,
            'p_val': p_val,
            'enrichment': enrichment,
            'observed': observed_hits,
            'expected': expected_hits
        })

    # 2. Process Trevino
    peak_dict, total_bp = get_trevino_peaks()
    if total_bp > 0:
        footprint_frac = total_bp / GENOME_SIZE_BP
        expected_hits = total_snps * footprint_frac
        observed_hits = check_intersection(variants, peak_dict)
        p_val = stats.binom.sf(observed_hits - 1, total_snps, footprint_frac)
        enrichment = observed_hits / expected_hits if expected_hits > 0 else 0
        
        sig = "***" if p_val < 0.001 else ""
        print(f"{'Fetal_Progenitors':<20} {total_bp/1e6:<10.1f} {expected_hits:<10.1f} {observed_hits:<10} {enrichment:<10.2f} {p_val:.2e} {sig}")
        
        results.append({
            'label': 'Fetal_Progenitors',
            'p_val': p_val,
            'enrichment': enrichment,
            'observed': observed_hits,
            'expected': expected_hits
        })

    # Save
    df = pd.DataFrame(results)
    out_file = CT_DIR / "cell_type_footprint_enrichment.csv"
    df.to_csv(out_file, index=False)
    print(f"\nSaved verified results to {out_file}")

if __name__ == "__main__":
    main()
