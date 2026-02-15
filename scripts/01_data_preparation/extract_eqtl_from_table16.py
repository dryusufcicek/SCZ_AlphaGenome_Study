"""
Extract eQTL Data from PGC3 Supplementary Table 16

Extracts TWAS/eQTL colocalization data from:
- eQTLGen (blood eQTLs, 15,659 gene-SNP pairs)
- PsychENCODE (brain eQTLs, 10,947 gene-SNP pairs)
- Fetal Brain (fetal brain eQTLs, 772 gene-SNP pairs)

Creates comprehensive variant-gene eQTL lookup table for prioritization.

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-15
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

# Path to supplementary tables
SUPP_DIR = Path("/Users/yusuf/Downloads/2020-08-14908C-s11 2")

def extract_table16_eqtls():
    """
    Extract eQTL/TWAS data from Supplementary Table 16.

    Returns
    -------
    dict
        Dictionary with keys 'eQTLGen', 'PsychENCODE', 'Fetal Brain'
        Each containing a DataFrame with eQTL evidence
    """
    print("="*80)
    print("EXTRACTING eQTL DATA FROM SUPPLEMENTARY TABLE 16")
    print("="*80)

    table_path = SUPP_DIR / "Supplementary Table 16.xlsx"

    if not table_path.exists():
        print(f"❌ Table 16 not found at {table_path}")
        return None

    print(f"\nReading: {table_path.name}")
    print("This contains TWAS/eQTL colocalization results")

    # Read all sheets
    sheets = pd.read_excel(table_path, sheet_name=None)

    eqtl_data = {}

    for sheet_name in ['eQTLGen', 'PsychENCODE', 'Fetal Brain']:
        if sheet_name not in sheets:
            print(f"⚠️  Sheet '{sheet_name}' not found")
            continue

        df = sheets[sheet_name]
        print(f"\n{'='*80}")
        print(f"SHEET: {sheet_name}")
        print(f"{'='*80}")
        print(f"Shape: {len(df):,} rows × {len(df.columns)} columns")

        # Display key columns
        print(f"\nColumns:")
        for col in df.columns:
            non_null = df[col].notna().sum()
            print(f"  - {col}: {non_null:,}/{len(df):,} non-null ({100*non_null/len(df):.1f}%)")

        # Show sample data
        print(f"\nSample data:")
        key_cols = ['Gene', 'topSNP', 'topSNP_chr', 'topSNP_bp', 'p_eQTL', 'b_eQTL', 'p_GWAS']
        available_cols = [c for c in key_cols if c in df.columns]
        print(df[available_cols].head(5).to_string())

        # Store
        eqtl_data[sheet_name] = df

    return eqtl_data

def create_comprehensive_eqtl_lookup(eqtl_data):
    """
    Create comprehensive eQTL lookup table combining all sources.

    Parameters
    ----------
    eqtl_data : dict
        Dictionary with eQTL data from different sources

    Returns
    -------
    pd.DataFrame
        Comprehensive lookup: variant → gene → eQTL evidence
    """
    print("\n" + "="*80)
    print("CREATING COMPREHENSIVE eQTL LOOKUP TABLE")
    print("="*80)

    # Load credible sets
    cred_file = config.CS_DIR / "scz_credsets_normalized.csv"
    print(f"\nLoading SCZ credible sets: {cred_file.name}")
    cred_df = pd.read_csv(cred_file)
    print(f"  {len(cred_df):,} credible set variants")

    # Combine all eQTL sources
    all_eqtl_records = []

    for source_name, df in eqtl_data.items():
        print(f"\nProcessing {source_name}...")

        # Standardize column names
        df_clean = df.copy()

        # Extract key information
        for idx, row in df_clean.iterrows():
            record = {
                'source': source_name,
                'gene_symbol': row['Gene'],
                'topSNP': row['topSNP'],
                'chr': row['topSNP_chr'],
                'bp': row['topSNP_bp'],
                'eqtl_pval': row['p_eQTL'],
                'eqtl_beta': row['b_eQTL'],
                'eqtl_se': row['se_eQTL'],
                'gwas_pval': row['p_GWAS'],
                'gwas_beta': row['b_GWAS'],
            }
            all_eqtl_records.append(record)

    # Create combined DataFrame
    eqtl_lookup = pd.DataFrame(all_eqtl_records)

    print(f"\n✓ Combined eQTL records: {len(eqtl_lookup):,}")
    print(f"  Unique genes: {eqtl_lookup['gene_symbol'].nunique():,}")
    print(f"  Unique SNPs: {eqtl_lookup['topSNP'].nunique():,}")

    # Distribution by source
    print(f"\nRecords by source:")
    print(eqtl_lookup['source'].value_counts())

    # Match with credible sets
    # Create rsID-based lookup (topSNP is in rsID format)
    print(f"\nMatching with credible set variants...")

    # Merge on SNP rsID
    cred_eqtl = cred_df.merge(
        eqtl_lookup,
        left_on='SNP',
        right_on='topSNP',
        how='left'
    )

    # Count matches
    n_matched = cred_eqtl['source'].notna().sum()
    print(f"  ✓ {n_matched:,}/{len(cred_df):,} credible set variants have eQTL evidence ({100*n_matched/len(cred_df):.1f}%)")

    # For matched variants, count genes
    matched_df = cred_eqtl[cred_eqtl['source'].notna()].copy()
    if len(matched_df) > 0:
        print(f"  Unique genes implicated: {matched_df['gene_symbol'].nunique():,}")
        print(f"\nTop 10 genes by number of eQTL variants:")
        top_genes = matched_df['gene_symbol'].value_counts().head(10)
        for gene, count in top_genes.items():
            print(f"    {gene}: {count} variants")

    # Save comprehensive lookup
    output_file = config.PROCESSED_DIR / "eqtl_data" / "eqtl_lookup_table16.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    eqtl_lookup.to_csv(output_file, index=False)
    print(f"\n✓ Saved comprehensive eQTL lookup: {output_file.name}")

    # Save matched credible set variants with eQTL info
    matched_file = config.PROCESSED_DIR / "eqtl_data" / "credsets_with_eqtl_evidence.csv"
    cred_eqtl.to_csv(matched_file, index=False)
    print(f"✓ Saved credible sets with eQTL evidence: {matched_file.name}")

    return eqtl_lookup, cred_eqtl

def analyze_eqtl_coverage():
    """Analyze eQTL coverage for key SCZ genes."""
    print("\n" + "="*80)
    print("ANALYZING eQTL COVERAGE FOR KEY SCZ GENES")
    print("="*80)

    # Load matched data
    matched_file = config.PROCESSED_DIR / "eqtl_data" / "credsets_with_eqtl_evidence.csv"
    cred_eqtl = pd.read_csv(matched_file)

    # Known SCZ genes
    known_scz_genes = [
        'CACNA1C', 'TCF4', 'ZNF804A', 'NRGN', 'GRIN2A', 'GRM3',
        'DRD2', 'SRR', 'CLCN3', 'SNAP91', 'FURIN', 'TSNARE1',
        'AS3MT', 'CNNM2', 'NT5C2', 'ITIH3', 'CNTN4', 'PCCB'
    ]

    print(f"\nChecking eQTL evidence for {len(known_scz_genes)} known SCZ genes:")

    for gene in known_scz_genes:
        gene_data = cred_eqtl[cred_eqtl['gene_symbol'] == gene]

        if len(gene_data) > 0:
            n_vars = len(gene_data)
            sources = gene_data['source'].unique()
            min_pval = gene_data['eqtl_pval'].min()
            print(f"  ✓ {gene:12s}: {n_vars:2d} eQTL variants, sources: {sources}, min p={min_pval:.2e}")
        else:
            print(f"  ✗ {gene:12s}: No eQTL evidence in Table 16")

def create_variant_gene_eqtl_weights():
    """
    Summarize eQTL coverage statistics.

    Since we have eQTL evidence for only 182/20,760 variants (0.9%),
    document the coverage for use in downstream gene aggregation.
    """
    print("\n" + "="*80)
    print("eQTL COVERAGE SUMMARY")
    print("="*80)

    # Load eQTL evidence
    eqtl_file = config.PROCESSED_DIR / "eqtl_data" / "credsets_with_eqtl_evidence.csv"
    print(f"\nLoading eQTL evidence: {eqtl_file.name}")
    eqtl_df = pd.read_csv(eqtl_file)

    # Summarize coverage
    total_variants = len(eqtl_df)
    variants_with_eqtl = eqtl_df['source'].notna().sum()
    pct_with_eqtl = 100 * variants_with_eqtl / total_variants

    print(f"\nCredible set variants with eQTL evidence:")
    print(f"  Total variants: {total_variants:,}")
    print(f"  With eQTL evidence: {variants_with_eqtl:,} ({pct_with_eqtl:.1f}%)")
    print(f"  Without eQTL evidence: {total_variants - variants_with_eqtl:,} ({100 - pct_with_eqtl:.1f}%)")

    # Get unique gene-variant pairs with eQTL
    eqtl_evidence = eqtl_df[eqtl_df['source'].notna()].copy()
    eqtl_pairs = eqtl_evidence[['SNP', 'gene_symbol']].drop_duplicates()

    print(f"\nUnique variant-gene eQTL pairs: {len(eqtl_pairs):,}")
    print(f"Unique genes with eQTL support: {eqtl_pairs['gene_symbol'].nunique():,}")
    print(f"Unique variants with eQTL support: {eqtl_pairs['SNP'].nunique():,}")

    # Distribution by source
    print(f"\neQTL evidence by source:")
    source_counts = eqtl_evidence['source'].value_counts()
    for source, count in source_counts.items():
        print(f"  {source:15s}: {count:4d} variant-gene pairs")

    # Save simple eQTL variant list
    eqtl_variants = eqtl_df[eqtl_df['source'].notna()][['SNP', 'CHR', 'BP', 'gene_symbol', 'source', 'eqtl_pval']].drop_duplicates()
    output_file = config.PROCESSED_DIR / "eqtl_data" / "variants_with_eqtl_support.csv"
    eqtl_variants.to_csv(output_file, index=False)

    print(f"\n✓ Saved eQTL variant list: {output_file.name}")

    return eqtl_pairs

def main():
    """Main execution."""
    print("╔" + "="*78 + "╗")
    print("║" + " "*18 + "eQTL EXTRACTION FROM SUPPLEMENTARY TABLE 16" + " "*17 + "║")
    print("╚" + "="*78 + "╝")

    # Step 1: Extract eQTL data from Table 16
    eqtl_data = extract_table16_eqtls()

    if eqtl_data is None:
        print("\n❌ Could not extract eQTL data")
        return 1

    # Step 2: Create comprehensive lookup
    eqtl_lookup, cred_eqtl = create_comprehensive_eqtl_lookup(eqtl_data)

    # Step 3: Analyze coverage for known SCZ genes
    analyze_eqtl_coverage()

    # Step 4: Summarize eQTL coverage
    eqtl_pairs = create_variant_gene_eqtl_weights()

    print("\n" + "="*80)
    print("✓ eQTL EXTRACTION COMPLETE")
    print("="*80)

    print(f"""
SUMMARY:
- Extracted eQTL data from 3 sources (eQTLGen, PsychENCODE, Fetal Brain)
- Total eQTL records: {len(eqtl_lookup):,}
- Credible set variants with eQTL evidence: {cred_eqtl['source'].notna().sum():,}/{len(cred_eqtl):,} (0.9%)
- Unique variant-gene eQTL pairs: {len(eqtl_pairs):,}

FILES CREATED:
- {config.PROCESSED_DIR / "eqtl_data" / "eqtl_lookup_table16.csv"}
- {config.PROCESSED_DIR / "eqtl_data" / "credsets_with_eqtl_evidence.csv"}
- {config.PROCESSED_DIR / "eqtl_data" / "variants_with_eqtl_support.csv"}

NEXT STEP:
Update aggregate_with_eqtl_priority.py to use real eQTL data from Table 16.
    """)

    return 0

if __name__ == "__main__":
    sys.exit(main())
