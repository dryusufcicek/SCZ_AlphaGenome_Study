"""
Script 04: H-MAGMA 3D Chromatin Validation

Validates AlphaGenome structural predictions against H-MAGMA Hi-C data.
Input: Adultbrain.transcript.annot (formatted as Gene SNP1 SNP2 ...)
Requires: brainspan/rows_metadata.csv for ID Mapping (ENSG -> Symbol)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact
import sys

# Add path for config
sys.path.append(str(Path(__file__).parent.parent))
BASE_DIR = Path(".").resolve()
if "AlphaGenome with SCZ" not in str(BASE_DIR):
    match = [p for p in Path("/Users/yusuf/AlphaGenome with SCZ").glob("*") if p.is_dir()]
    if match:
        BASE_DIR = Path("/Users/yusuf/AlphaGenome with SCZ")

DATA_DIR = BASE_DIR / "data"
VALIDATION_DIR = DATA_DIR / "validation"
ANNOT_FILE = VALIDATION_DIR / "Annotation_Files" / "Adultbrain.transcript.annot"
PROCESSED_DIR = BASE_DIR / "scz_hypothesis_testing/data/processed"
METADATA_FILE = BASE_DIR / "scz_hypothesis_testing/data/brainspan/rows_metadata.csv"

def load_gene_mapping():
    """Load ENSG -> Symbol mapping from BrainSpan metadata."""
    if not METADATA_FILE.exists():
        print(f"Error: Metadata file not found at {METADATA_FILE}")
        return None
    
    print("Loading ID mapping from BrainSpan metadata...")
    try:
        df = pd.read_csv(METADATA_FILE)
        # Expected cols: ensembl_gene_id, gene_symbol
        # Handle Clean Ensembl IDs (remove version)
        
        mapping = {}
        for _, row in df.iterrows():
            ensg = str(row['ensembl_gene_id']).split('.')[0]
            symbol = str(row['gene_symbol']).upper()
            mapping[ensg] = symbol
            
        print(f"Loaded {len(mapping)} gene mappings")
        return mapping
    except Exception as e:
        print(f"Error loading mapping: {e}")
        return None

def load_alphagenome_predictions():
    """Load AlphaGenome predicted target genes."""
    scores_file = PROCESSED_DIR / "robust_modality_scores.csv"
    if not scores_file.exists():
        print(f"Error: AlphaGenome scores not found at {scores_file}")
        return None, None
        
    df = pd.read_csv(scores_file)
    
    # Extract SNP-Gene pairs
    pairs = set()
    snps_scored = set()
    
    for _, row in df.iterrows():
        snp = row['SNP']
        snps_scored.add(snp)
        if pd.isna(row.get('genes', '')):
            continue
            
        genes_str = str(row['genes'])
        for entry in genes_str.split(';'):
            if ':' in entry:
                gene, score = entry.split(':')
                if float(score) >= 1.0: 
                    clean_gene = gene.split('.')[0].upper()
                    pairs.add((snp, clean_gene))
    
    return pairs, snps_scored

def load_hmagma_interactions(id_map):
    """Load H-MAGMA annotations (Adult Brain) and convert to Symbols."""
    if not ANNOT_FILE.exists():
        print(f"\n⚠️  H-MAGMA file not found: {ANNOT_FILE}")
        return None
    
    if id_map is None:
        print("Warning: No ID mapping provided. Using raw IDs.")
    
    print("Loading H-MAGMA interactions...")
    hmagma_pairs = set()
    
    mapped_count = 0
    total_genes = 0
    
    with open(ANNOT_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            
            gene_info = parts[0] # ENSG ID
            clean_ensg = gene_info.split('.')[0]
            total_genes += 1
            
            # Map to Symbol
            if id_map and clean_ensg in id_map:
                gene_symbol = id_map[clean_ensg]
                mapped_count += 1
            else:
                # Fallback: keep ENSG if no map (unlikely to match but valid)
                gene_symbol = clean_ensg
            
            snps = parts[2:]
            
            for snp in snps:
                hmagma_pairs.add((snp, gene_symbol))
    
    print(f"Mapped {mapped_count}/{total_genes} H-MAGMA genes to Symbols")
    return hmagma_pairs

def main():
    print("=" * 70)
    print("H-MAGMA 3D CHROMATIN VALIDATION (With ID Mapping)")
    print("=" * 70)
    
    # Load Mapping
    id_map = load_gene_mapping()
    
    # Load AlphaGenome predictions
    print("Loading AlphaGenome predictions...")
    ag_pairs, ag_snps_universe = load_alphagenome_predictions()
    if ag_pairs is None: return
        
    print(f"AlphaGenome: {len(ag_pairs)} predicted SNP-Gene links")
    
    # Load H-MAGMA
    hmagma_pairs = load_hmagma_interactions(id_map)
    if hmagma_pairs is None: return
        
    print(f"H-MAGMA:     {len(hmagma_pairs)} total validated links")
    
    # Validation Statistics
    print("\nCalculating overlap...")
    
    overlap = ag_pairs.intersection(hmagma_pairs)
    
    precision = len(overlap) / len(ag_pairs) if len(ag_pairs) > 0 else 0
    
    hmagma_subset = {p for p in hmagma_pairs if p[0] in ag_snps_universe}
    recall = len(overlap) / len(hmagma_subset) if len(hmagma_subset) > 0 else 0
    
    print("-" * 30)
    print(f"VALIDATION RESULTS")
    print("-" * 30)
    print(f"Overlap (Confirmed Links): {len(overlap)}")
    print(f"Precision (AG -> HMAGMA):  {precision:.2%}")
    print(f"Recall (HMAGMA -> AG):     {recall:.2%}")
    
    if len(overlap) > 0:
        out_file = DATA_DIR.parent / "scz_hypothesis_testing/results/hmagma_validated_links.csv"
        with open(out_file, "w") as f:
            f.write("SNP,Gene\n")
            for s, g in overlap:
                f.write(f"{s},{g}\n")
        print(f"\nValidated links saved to {out_file}")

if __name__ == "__main__":
    main()
