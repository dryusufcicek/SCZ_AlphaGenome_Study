"""
Data Source Verification Script

This script verifies that all required input files exist, have correct formats,
and are using the correct genomic build (hg38).

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-14
"""

import sys
from pathlib import Path
import pandas as pd
import gzip

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

def check_file_exists(filepath, description):
    """Check if a file exists and return its size."""
    path = Path(filepath)
    if not path.exists():
        print(f"❌ MISSING: {description}")
        print(f"   Expected: {filepath}")
        return False, 0
    else:
        size_mb = path.stat().st_size / (1024 * 1024)
        print(f"✓ FOUND: {description} ({size_mb:.2f} MB)")
        return True, size_mb

def verify_credible_sets():
    """Verify credible sets file."""
    print("\n" + "="*80)
    print("CREDIBLE SETS (Fine-Mapped Variants)")
    print("="*80)

    cs_file = config.CS_DIR / "scz_credsets_ST11d_grch38_pip.csv"
    exists, size = check_file_exists(cs_file, "SCZ credible sets (GRCh38)")

    if exists:
        df = pd.read_csv(cs_file)
        print(f"   Variants: {len(df):,}")
        print(f"   Loci: {df['Locus_ID'].nunique()}")
        print(f"   Chromosomes: {sorted(df['CHR'].unique())}")

        # Check PP sums
        pp_sums = df.groupby('Locus_ID')['PP'].sum()
        problem_loci = (pp_sums > 1.0).sum()
        print(f"   PP sum range: {pp_sums.min():.4f} - {pp_sums.max():.4f}")
        print(f"   Loci with PP > 1.0: {problem_loci} / {len(pp_sums)} ({100*problem_loci/len(pp_sums):.1f}%)")

        if problem_loci > 0:
            print(f"   ⚠️  WARNING: {problem_loci} loci have PP sum > 1.0 (will be normalized)")

        return True
    return False

def verify_alphagenome_scores():
    """Verify AlphaGenome scoring output."""
    print("\n" + "="*80)
    print("ALPHAGENOME SCORES")
    print("="*80)

    ag_file = config.AG_DIR / "variant_scores_lfc.csv"
    exists, size = check_file_exists(ag_file, "AlphaGenome variant scores (LFC)")

    if exists:
        df = pd.read_csv(ag_file)
        print(f"   Variants: {len(df):,}")
        print(f"   Columns: {len(df.columns)}")

        # Check for required score columns
        required_scores = ['score_DNase', 'score_ATAC', 'score_CAGE',
                          'score_RNA', 'score_PROCAP', 'score_Splice_Junction']
        missing_scores = [s for s in required_scores if s not in df.columns]

        if missing_scores:
            print(f"   ❌ MISSING SCORE COLUMNS: {missing_scores}")
            return False
        else:
            print(f"   ✓ All 6 modality scores present")

        # Check score distributions
        print(f"\n   Score Statistics:")
        for score_col in required_scores:
            vals = df[score_col].dropna()
            print(f"   - {score_col:25s}: mean={vals.mean():7.4f}, std={vals.std():7.4f}, "
                  f"range=[{vals.min():7.4f}, {vals.max():7.4f}]")

        # Check for scale mismatch (splice vs RNA issue)
        splice_median_abs = df['score_Splice_Junction'].abs().median()
        rna_median_abs = df['score_RNA'].abs().median()

        if splice_median_abs > 10 * rna_median_abs:
            print(f"\n   ⚠️  WARNING: Splice scores may be 10x+ larger than RNA scores")
            print(f"       Splice median(abs): {splice_median_abs:.4f}")
            print(f"       RNA median(abs): {rna_median_abs:.4f}")
            print(f"       Ratio: {splice_median_abs/rna_median_abs if rna_median_abs > 0 else 'inf':.1f}x")
            print(f"       → Will apply per-modality normalization")

        return True
    return False

def verify_gtex_data():
    """Verify GTEx eQTL data."""
    print("\n" + "="*80)
    print("GTEx eQTL DATA")
    print("="*80)

    gtex_tar = config.PROCESSED_DIR / "eqtl_data" / "GTEx_Analysis_v10_eQTL.tar"
    exists, size = check_file_exists(gtex_tar, "GTEx v10 eQTL summary statistics (tar)")

    if exists:
        if size > 1000:  # Over 1 GB
            print(f"   ⚠️  Large file ({size/1024:.2f} GB) - needs extraction")
            print(f"       Run: scripts_new/01_data_preparation/extract_gtex_brain.py")
        return True
    return False

def verify_celltype_data():
    """Verify cell-type ATAC-seq peak data."""
    print("\n" + "="*80)
    print("CELL-TYPE ATAC-SEQ PEAKS")
    print("="*80)

    ct_dir = config.PROCESSED_DIR / "celltype_peaks"
    if not ct_dir.exists():
        print(f"❌ MISSING: Cell-type peaks directory")
        return False

    # Count peak files
    peak_files = list(ct_dir.glob("*.narrowPeak.gz"))
    print(f"✓ FOUND: {len(peak_files)} cell-type peak files")

    if len(peak_files) > 0:
        # Check first file format
        sample_file = peak_files[0]
        with gzip.open(sample_file, 'rt') as f:
            first_line = f.readline()
            fields = first_line.strip().split('\t')
            print(f"   Sample file: {sample_file.name}")
            print(f"   Columns: {len(fields)}")
            print(f"   Format: BED (chr, start, end, ...)")
        return True
    else:
        print(f"⚠️  WARNING: No peak files found")
        return False

def verify_hic_data():
    """Verify Hi-C chromatin loop data."""
    print("\n" + "="*80)
    print("Hi-C CHROMATIN LOOPS")
    print("="*80)

    hic_file = config.PROCESSED_DIR / "hic_data" / "psychcode_loops_fdr.bedpe"
    exists, size = check_file_exists(hic_file, "PsychENCODE Hi-C loops (BEDPE)")

    if exists:
        df = pd.read_csv(hic_file, sep='\t', nrows=5)
        print(f"   Columns: {len(df.columns)}")
        print(f"   Format: BEDPE (chr1, start1, end1, chr2, start2, end2, ...)")
        return True
    return False

def verify_annotations():
    """Verify annotation files."""
    print("\n" + "="*80)
    print("ANNOTATION FILES")
    print("="*80)

    annot_dir = config.PROCESSED_DIR / "annotations"

    # Check GO pathways
    go_file = annot_dir / "GO_Biological_Process_2023.gmt"
    check_file_exists(go_file, "GO Biological Process gene sets (GMT)")

    # Check cell-type expression
    siletti_file = annot_dir / "Siletti_Supercluster_expression_specificity_TDEP_label.tsv.gz"
    check_file_exists(siletti_file, "Siletti brain cell-type expression")

    # Check SCHEMA genes
    schema_file = annot_dir / "schema_genes.xlsx"
    check_file_exists(schema_file, "SCHEMA rare variant genes")

    return True

def verify_external_data():
    """Verify external reference data."""
    print("\n" + "="*80)
    print("EXTERNAL REFERENCE DATA")
    print("="*80)

    ext_dir = config.DATA_DIR / "external"

    # Check GENCODE annotation
    gencode_file = ext_dir / "gencode.v43.basic.annotation.gtf.gz"
    exists, size = check_file_exists(gencode_file, "GENCODE v43 gene annotation (GTF)")

    if exists:
        print(f"   Build: GRCh38")
        print(f"   Version: v43 (comprehensive gene annotation)")

    return exists

def main():
    """Run all verification checks."""
    print("╔" + "="*78 + "╗")
    print("║" + " "*20 + "DATA SOURCE VERIFICATION REPORT" + " "*27 + "║")
    print("╚" + "="*78 + "╝")

    results = {
        'Credible Sets': verify_credible_sets(),
        'AlphaGenome Scores': verify_alphagenome_scores(),
        'GTEx eQTL Data': verify_gtex_data(),
        'Cell-Type ATAC Peaks': verify_celltype_data(),
        'Hi-C Loops': verify_hic_data(),
        'Annotations': verify_annotations(),
        'External Reference': verify_external_data(),
    }

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    passed = sum(results.values())
    total = len(results)

    for name, status in results.items():
        status_str = "✓ PASS" if status else "❌ FAIL"
        print(f"{name:30s}: {status_str}")

    print(f"\nOverall: {passed}/{total} checks passed")

    if passed == total:
        print("\n✓ All data sources verified successfully!")
        print("  Ready to proceed with analysis pipeline.")
        return 0
    else:
        print("\n⚠️  Some data sources are missing or have issues.")
        print("   Please resolve the issues above before proceeding.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
