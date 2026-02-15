"""
Extract eQTL Information from PGC3 SCZ GWAS 2022 Supplementary Tables

Extracts brain eQTL annotations from the PGC3 SCZ GWAS 2022 supplementary tables
to create a comprehensive variant-gene eQTL lookup table.

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

def explore_supplementary_tables():
    """Explore supplementary tables to identify eQTL columns."""
    print("="*80)
    print("EXPLORING SUPPLEMENTARY TABLES FOR eQTL DATA")
    print("="*80)

    # Tables likely to contain eQTL info:
    # - Table 11: Likely detailed variant annotations (7.7 MB)
    # - Table 28: Very large, likely comprehensive annotations (12 MB)
    # - Table 7: Medium size, could be locus-level annotations (548 KB)

    tables_to_check = [
        ("Supplementary Table 7.xlsx", "Locus-level annotations"),
        ("Supplementary Table 11.xlsx", "Detailed variant annotations"),
        ("Supplementary Table 28.xlsx", "Comprehensive annotations"),
    ]

    for table_name, description in tables_to_check:
        table_path = SUPP_DIR / table_name
        if not table_path.exists():
            print(f"\n⚠️  {table_name} not found")
            continue

        print(f"\n{'='*80}")
        print(f"TABLE: {table_name}")
        print(f"Description: {description}")
        print(f"{'='*80}")

        try:
            # Read first few rows to understand structure
            df = pd.read_excel(table_path, nrows=5)

            print(f"\nShape: {df.shape[0]} rows (sample), {df.shape[1]} columns")
            print(f"\nColumn names:")
            for i, col in enumerate(df.columns, 1):
                print(f"  {i:2d}. {col}")

            # Look for eQTL-related columns
            eqtl_cols = [col for col in df.columns if any(
                keyword in str(col).lower()
                for keyword in ['eqtl', 'expression', 'qtl', 'tissue', 'gene', 'brain']
            )]

            if eqtl_cols:
                print(f"\n✓ Found {len(eqtl_cols)} potentially relevant columns:")
                for col in eqtl_cols:
                    print(f"  - {col}")
                    # Show sample values
                    sample_vals = df[col].dropna().head(3).tolist()
                    if sample_vals:
                        print(f"    Sample: {sample_vals}")
            else:
                print("\n❌ No obvious eQTL columns found")

        except Exception as e:
            print(f"❌ Error reading {table_name}: {e}")

def extract_table_11_eqtls():
    """
    Extract eQTL data from Supplementary Table 11.
    This is likely the detailed variant-level annotations.
    """
    print("\n" + "="*80)
    print("EXTRACTING eQTL DATA FROM SUPPLEMENTARY TABLE 11")
    print("="*80)

    table_path = SUPP_DIR / "Supplementary Table 11.xlsx"

    if not table_path.exists():
        print(f"❌ Table 11 not found at {table_path}")
        return None

    print(f"\nReading: {table_path.name}")
    print("This may take a moment (7.7 MB file)...")

    try:
        # Read the full table
        df = pd.read_excel(table_path)

        print(f"✓ Loaded: {len(df):,} rows × {len(df.columns)} columns")

        # Display structure
        print(f"\nColumn names:")
        for i, col in enumerate(df.columns, 1):
            print(f"  {i:2d}. {col}")

        # Check for key columns we need
        required_cols = ['SNP', 'rsid', 'CHR', 'BP', 'variant', 'A1', 'A2']
        found_id_cols = [col for col in df.columns if any(req.lower() in col.lower() for req in required_cols)]

        print(f"\nVariant ID columns found: {found_id_cols}")

        # Look for eQTL annotations
        eqtl_cols = [col for col in df.columns if any(
            keyword in str(col).lower()
            for keyword in ['eqtl', 'expression', 'qtl', 'tissue', 'gene_name', 'gene_id', 'brain']
        )]

        print(f"\neQTL-related columns ({len(eqtl_cols)}):")
        for col in eqtl_cols:
            non_null = df[col].notna().sum()
            print(f"  - {col}: {non_null:,}/{len(df):,} non-null ({100*non_null/len(df):.1f}%)")

        return df

    except Exception as e:
        print(f"❌ Error reading Table 11: {e}")
        return None

def extract_table_28_eqtls():
    """
    Extract eQTL data from Supplementary Table 28.
    This is the largest table and likely most comprehensive.
    """
    print("\n" + "="*80)
    print("EXTRACTING eQTL DATA FROM SUPPLEMENTARY TABLE 28")
    print("="*80)

    table_path = SUPP_DIR / "Supplementary Table 28.xlsx"

    if not table_path.exists():
        print(f"❌ Table 28 not found at {table_path}")
        return None

    print(f"\nReading: {table_path.name}")
    print("This may take a moment (12 MB file)...")

    try:
        # Read the full table
        df = pd.read_excel(table_path)

        print(f"✓ Loaded: {len(df):,} rows × {len(df.columns)} columns")

        # Display structure
        print(f"\nColumn names:")
        for i, col in enumerate(df.columns, 1):
            print(f"  {i:2d}. {col}")

        # Look for eQTL annotations
        eqtl_cols = [col for col in df.columns if any(
            keyword in str(col).lower()
            for keyword in ['eqtl', 'expression', 'qtl', 'tissue', 'gene', 'brain']
        )]

        print(f"\neQTL-related columns ({len(eqtl_cols)}):")
        for col in eqtl_cols:
            non_null = df[col].notna().sum()
            print(f"  - {col}: {non_null:,}/{len(df):,} non-null ({100*non_null/len(df):.1f}%)")

            # Show sample values
            if non_null > 0:
                sample_vals = df[col].dropna().head(3).tolist()
                print(f"    Sample: {sample_vals}")

        return df

    except Exception as e:
        print(f"❌ Error reading Table 28: {e}")
        return None

def create_eqtl_lookup_from_supplementary(supp_df, table_name):
    """
    Create eQTL lookup table from supplementary table data.

    Parameters
    ----------
    supp_df : pd.DataFrame
        Supplementary table with eQTL annotations
    table_name : str
        Name of the source table
    """
    print("\n" + "="*80)
    print(f"CREATING eQTL LOOKUP FROM {table_name}")
    print("="*80)

    # Load our credible sets for matching
    cred_file = config.CS_DIR / "scz_credsets_normalized.csv"
    print(f"\nLoading SCZ credible sets: {cred_file.name}")
    cred_df = pd.read_csv(cred_file)
    print(f"  {len(cred_df):,} credible set variants")

    # Try to identify variant ID columns in supplementary table
    # Common column names: SNP, rsID, variant_id, MarkerName, etc.
    variant_id_cols = [col for col in supp_df.columns if any(
        keyword in str(col).lower()
        for keyword in ['snp', 'rsid', 'variant', 'marker']
    )]

    print(f"\nVariant ID columns in {table_name}:")
    for col in variant_id_cols:
        print(f"  - {col}")
        print(f"    Sample: {supp_df[col].dropna().head(3).tolist()}")

    # Try to identify gene columns
    gene_cols = [col for col in supp_df.columns if any(
        keyword in str(col).lower()
        for keyword in ['gene', 'symbol', 'ensembl']
    )]

    print(f"\nGene columns in {table_name}:")
    for col in gene_cols:
        non_null = supp_df[col].notna().sum()
        print(f"  - {col}: {non_null:,} non-null")
        if non_null > 0:
            print(f"    Sample: {supp_df[col].dropna().head(3).tolist()}")

    # Try to identify eQTL evidence columns
    eqtl_evidence_cols = [col for col in supp_df.columns if any(
        keyword in str(col).lower()
        for keyword in ['eqtl', 'expression', 'brain', 'tissue', 'pval', 'p-value', 'fdr']
    )]

    print(f"\neQTL evidence columns in {table_name}:")
    for col in eqtl_evidence_cols:
        non_null = supp_df[col].notna().sum()
        if non_null > 0:
            print(f"  - {col}: {non_null:,} non-null")
            print(f"    Sample: {supp_df[col].dropna().head(3).tolist()}")

    # Note: Actual merging will depend on the specific column structure
    # This exploration will help us identify the right approach

    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("""
Based on the exploration above, we need to:
1. Identify the correct variant ID column for matching
2. Identify gene annotation columns (gene symbol, Ensembl ID)
3. Identify eQTL evidence columns (tissue, p-value, effect size)
4. Merge with our credible sets on variant ID
5. Create lookup table: variant → gene → eQTL_evidence
    """)

    return supp_df

def main():
    """Main execution."""
    print("╔" + "="*78 + "╗")
    print("║" + " "*20 + "eQTL EXTRACTION FROM SUPPLEMENTARY TABLES" + " "*17 + "║")
    print("╚" + "="*78 + "╝")

    # Step 1: Explore tables to understand structure
    print("\nSTEP 1: EXPLORING TABLE STRUCTURES")
    explore_supplementary_tables()

    # Step 2: Extract from most promising tables
    print("\n\nSTEP 2: EXTRACTING DATA FROM KEY TABLES")

    # Try Table 11 first (likely variant-level annotations)
    df_table11 = extract_table_11_eqtls()

    # Try Table 28 (comprehensive annotations)
    df_table28 = extract_table_28_eqtls()

    # Step 3: Create lookup table from the best source
    if df_table11 is not None:
        create_eqtl_lookup_from_supplementary(df_table11, "Table 11")
    elif df_table28 is not None:
        create_eqtl_lookup_from_supplementary(df_table28, "Table 28")
    else:
        print("\n❌ Could not extract data from supplementary tables")
        return 1

    print("\n" + "="*80)
    print("✓ EXPLORATION COMPLETE")
    print("="*80)
    print("\nNext: Based on the column structure identified above,")
    print("we can create a targeted extraction script to build the eQTL lookup table.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
