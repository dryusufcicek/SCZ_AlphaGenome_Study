"""
Gene ID Mapping Script

Maps Ensembl gene IDs to gene symbols using GENCODE annotation.
This fixes the 62.5% data loss from gene name mismatch.

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-14
"""

import sys
from pathlib import Path
import pandas as pd
import gzip
import re

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            parts = item.split(' ', 1)
            if len(parts) == 2:
                key, value = parts
                attrs[key] = value.strip('"')
    return attrs

def create_gene_id_mapping():
    """
    Create mapping from Ensembl gene IDs to gene symbols.

    Returns
    -------
    pd.DataFrame
        Mapping with columns: gene_id (Ensembl), gene_name (symbol)
    """
    print("="*80)
    print("CREATING GENE ID MAPPING")
    print("="*80)

    gtf_file = config.DATA_DIR / "external" / "gencode.v43.basic.annotation.gtf.gz"

    print(f"\nParsing GENCODE GTF: {gtf_file.name}")

    mappings = []

    with gzip.open(gtf_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            attributes = parse_gtf_attributes(fields[8])

            gene_id = attributes.get('gene_id', '')
            gene_name = attributes.get('gene_name', '')

            if gene_id and gene_name:
                # Remove version from Ensembl ID (e.g., ENSG00000123456.7 → ENSG00000123456)
                gene_id_no_version = gene_id.split('.')[0]

                mappings.append({
                    'gene_id': gene_id,
                    'gene_id_no_version': gene_id_no_version,
                    'gene_name': gene_name
                })

            if line_num % 100000 == 0:
                print(f"  Processed {line_num:,} lines, found {len(mappings):,} genes...", end='\r')

    print(f"\n  ✓ Parsed {line_num:,} lines total")

    df = pd.DataFrame(mappings)
    print(f"  ✓ Created {len(df):,} gene ID mappings")

    # Remove duplicates (keep first)
    df = df.drop_duplicates(subset=['gene_id_no_version'], keep='first')
    print(f"  ✓ {len(df):,} unique gene IDs")

    return df

def apply_mapping_to_gene_scores():
    """
    Apply gene ID mapping to AlphaGenome gene scores.

    This will:
    1. Parse genes in AlphaGenome scores (mix of Ensembl IDs and symbols)
    2. Map Ensembl IDs to symbols
    3. Create standardized gene names
    """
    print("\n" + "="*80)
    print("APPLYING MAPPING TO GENE SCORES")
    print("="*80)

    # Load gene scores
    gene_scores_file = config.PROCESSED_DIR / "gene_scores" / "gene_scores_aggregated.csv"
    print(f"\nLoading gene scores from: {gene_scores_file.name}")
    gene_scores = pd.read_csv(gene_scores_file)
    print(f"  ✓ Loaded {len(gene_scores):,} genes")

    # Load mapping
    mapping_file = config.PROCESSED_DIR / "annotations" / "gene_id_to_symbol_mapping.csv"
    print(f"\nLoading gene ID mapping from: {mapping_file.name}")
    mapping = pd.read_csv(mapping_file)
    print(f"  ✓ Loaded {len(mapping):,} mappings")

    # Create lookup dictionary
    id_to_symbol = {}
    for idx, row in mapping.iterrows():
        # Map both versioned and unversioned IDs
        id_to_symbol[row['gene_id']] = row['gene_name']
        id_to_symbol[row['gene_id_no_version']] = row['gene_name']

    print(f"\n  Mapping dictionary size: {len(id_to_symbol):,}")

    # Apply mapping
    print(f"\nMapping gene identifiers...")

    gene_scores['gene_original'] = gene_scores['gene']
    gene_scores['gene_symbol'] = gene_scores['gene'].apply(
        lambda x: id_to_symbol.get(x, x)  # If not found in mapping, keep original
    )

    # Count how many were mapped
    ensembl_pattern = re.compile(r'^ENSG\d+')

    original_ensembl = gene_scores['gene_original'].apply(lambda x: bool(ensembl_pattern.match(str(x)))).sum()
    after_ensembl = gene_scores['gene_symbol'].apply(lambda x: bool(ensembl_pattern.match(str(x)))).sum()

    mapped = original_ensembl - after_ensembl

    print(f"\n  Original Ensembl IDs: {original_ensembl:,}")
    print(f"  Successfully mapped: {mapped:,} ({100*mapped/original_ensembl:.1f}%)")
    print(f"  Remaining Ensembl IDs: {after_ensembl:,}")

    # Use gene_symbol as the new gene column
    gene_scores['gene'] = gene_scores['gene_symbol']
    gene_scores = gene_scores.drop(columns=['gene_symbol'])

    # Save updated gene scores
    output_file = config.PROCESSED_DIR / "gene_scores" / "gene_scores_mapped.csv"
    print(f"\nSaving mapped gene scores to: {output_file}")
    gene_scores.to_csv(output_file, index=False)
    print(f"  ✓ Saved {len(gene_scores):,} genes")

    return gene_scores

def verify_improvement():
    """
    Verify that gene ID mapping improved overlap with gene universe.
    """
    print("\n" + "="*80)
    print("VERIFYING IMPROVEMENT")
    print("="*80)

    # Load data
    gene_scores_old = pd.read_csv(config.PROCESSED_DIR / "gene_scores" / "gene_scores_aggregated.csv")
    gene_scores_new = pd.read_csv(config.PROCESSED_DIR / "gene_scores" / "gene_scores_mapped.csv")
    gene_universe = pd.read_csv(config.PROCESSED_DIR / "annotations" / "gene_universe_protein_coding.csv")

    # Calculate overlaps
    universe_genes = set(gene_universe['gene_name'])

    overlap_old = len(set(gene_scores_old['gene']) & universe_genes)
    overlap_new = len(set(gene_scores_new['gene']) & universe_genes)

    print(f"\nOVERLAP WITH GENE UNIVERSE:")
    print(f"  Before mapping: {overlap_old:,} / {len(gene_scores_old):,} ({100*overlap_old/len(gene_scores_old):.1f}%)")
    print(f"  After mapping:  {overlap_new:,} / {len(gene_scores_new):,} ({100*overlap_new/len(gene_scores_new):.1f}%)")
    print(f"  Improvement:    +{overlap_new - overlap_old:,} genes ({100*(overlap_new-overlap_old)/len(gene_scores_old):.1f}%)")

    # Check CACNA1C
    print(f"\nCAC NA1C STATUS:")
    cacna1c_in_old = 'CACNA1C' in gene_scores_old['gene'].values
    cacna1c_in_new = 'CACNA1C' in gene_scores_new['gene'].values

    print(f"  In original scores: {cacna1c_in_old}")
    print(f"  In mapped scores:   {cacna1c_in_new}")

    if cacna1c_in_new:
        cacna1c_row = gene_scores_new[gene_scores_new['gene'] == 'CACNA1C'].iloc[0]
        if 'gene_original' in cacna1c_row:
            print(f"  Original ID: {cacna1c_row['gene_original']}")
        print(f"  Composite score: {cacna1c_row['composite_score']:.6f}")

def main():
    """Create gene ID mapping and apply it."""
    print("╔" + "="*78 + "╗")
    print("║" + " "*26 + "GENE ID MAPPING" + " "*37 + "║")
    print("╚" + "="*78 + "╝")

    # Step 1: Create mapping from GENCODE
    mapping = create_gene_id_mapping()

    # Save mapping
    output_file = config.PROCESSED_DIR / "annotations" / "gene_id_to_symbol_mapping.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving mapping to: {output_file}")
    mapping.to_csv(output_file, index=False)
    print(f"  ✓ Saved {len(mapping):,} mappings")

    # Step 2: Apply mapping to gene scores
    gene_scores_mapped = apply_mapping_to_gene_scores()

    # Step 3: Verify improvement
    verify_improvement()

    print("\n" + "="*80)
    print("✓ GENE ID MAPPING COMPLETE")
    print("="*80)
    print(f"\nNext: Use gene_scores_mapped.csv for downstream analyses")
    print(f"This will increase gene overlap and improve statistical power.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
