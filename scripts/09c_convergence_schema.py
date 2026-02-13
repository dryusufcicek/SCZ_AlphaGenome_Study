#!/usr/bin/env python3
"""Step 9c: Convergence analysis with SCHEMA rare variant genes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
EXT_DIR = BASE_DIR / "data" / "external"
GENE_SCORES_FILE = DATA_DIR / "gene_scores" / "gene_weighted_scores.csv"
SCHEMA_FILE = EXT_DIR / "schema_genes.xlsx"
SCHEMA_FALLBACK_FILE = BASE_DIR.parent / "data" / "validation" / "scz2022-Extended-Data-Table1.xlsx"
OUT_DIR = DATA_DIR / "validation" / "schema_convergence"
TABLE_OUT = BASE_DIR / "results" / "tables" / "schema_convergence_overlap.csv"
UNIVERSE_FILE = DATA_DIR / "gene_universe.csv"


def find_column(df: pd.DataFrame, patterns: list[str]) -> str | None:
    cols = {c.lower(): c for c in df.columns}
    for pattern in patterns:
        for low, original in cols.items():
            if pattern in low:
                return original
    return None


def genes_from_schema_workbook(filepath: Path) -> set[str] | None:
    if not filepath.exists():
        return None

    print(f"Loading SCHEMA workbook: {filepath}")
    xls = pd.ExcelFile(filepath)
    for sheet in xls.sheet_names:
        df = pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl")
        gene_col = find_column(df, ["gene", "symbol"])
        if gene_col is None:
            continue

        fdr_col = find_column(df, ["fdr", "q_"])
        p_col = find_column(df, ["p", "pvalue", "p_value", "p-meta", "p_meta"])

        work = df.copy()
        if fdr_col is not None:
            work[fdr_col] = pd.to_numeric(work[fdr_col], errors="coerce")
            sig = work[work[fdr_col] < 0.05]
            if not sig.empty:
                genes = set(sig[gene_col].dropna().astype(str))
                print(f"Using sheet '{sheet}' with FDR<0.05 ({len(genes)} genes)")
                return genes

        if p_col is not None:
            work[p_col] = pd.to_numeric(work[p_col], errors="coerce")
            sig = work[work[p_col] < 0.05]
            if not sig.empty:
                genes = set(sig[gene_col].dropna().astype(str))
                print(f"Using sheet '{sheet}' with P<0.05 ({len(genes)} genes)")
                return genes

        genes = set(work[gene_col].dropna().astype(str))
        if genes:
            print(f"Using sheet '{sheet}' full gene list ({len(genes)} genes)")
            return genes

    return None


def genes_from_extended_table(filepath: Path) -> set[str] | None:
    if not filepath.exists():
        return None

    print(f"Loading SCHEMA fallback from extended table: {filepath}")
    xls = pd.ExcelFile(filepath)
    # Prefer detailed sheet if available.
    candidates = [s for s in xls.sheet_names if "ST12" in s or "criteria" in s.lower()] + xls.sheet_names

    for sheet in candidates:
        df = pd.read_excel(filepath, sheet_name=sheet)
        schema_col = find_column(df, ["schema"])
        gene_col = find_column(df, ["symbol", "gene"])
        if schema_col is None or gene_col is None:
            continue

        schema_flag = pd.to_numeric(df[schema_col], errors="coerce").fillna(0)
        genes = set(df.loc[schema_flag > 0, gene_col].dropna().astype(str))
        if genes:
            print(f"Using sheet '{sheet}' SCHEMA-flagged genes ({len(genes)} genes)")
            return genes

    return None


def load_schema_genes() -> set[str] | None:
    # Prefer extended data table with explicit per-gene SCHEMA flags when available.
    genes = genes_from_extended_table(SCHEMA_FALLBACK_FILE)
    if genes:
        return genes

    genes = genes_from_schema_workbook(SCHEMA_FILE)
    if genes:
        # Guard against summary tables where "gene" columns are numeric counts.
        parsed = pd.Series(list(genes), dtype=str)
        numeric_frac = pd.to_numeric(parsed, errors="coerce").notna().mean()
        if numeric_frac < 0.2:
            return genes
        print("Resolved SCHEMA workbook content appears numeric-only; skipping.")

    print("No SCHEMA gene list could be resolved from available files.")
    return None


def main() -> None:
    print("=" * 72)
    print("STEP 9c: CONVERGENCE ANALYSIS (AlphaGenome x SCHEMA)")
    print("=" * 72)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_OUT.parent.mkdir(parents=True, exist_ok=True)

    if not GENE_SCORES_FILE.exists():
        print(f"Gene score file missing: {GENE_SCORES_FILE}")
        return

    our_df = pd.read_csv(GENE_SCORES_FILE)
    if "gene" not in our_df.columns or "weighted_score" not in our_df.columns:
        print("gene_weighted_scores.csv must include gene and weighted_score columns.")
        return

    our_df = our_df.sort_values("weighted_score", ascending=False)
    top_n = min(500, len(our_df))
    our_top_genes = set(our_df.head(top_n)["gene"].astype(str))
    print(f"Top regulatory gene set size: {len(our_top_genes)}")

    schema_genes = load_schema_genes()
    if not schema_genes:
        print("SCHEMA genes unavailable; aborting convergence analysis.")
        return

    overlap = sorted(our_top_genes.intersection(schema_genes))
    k = len(overlap)

    if UNIVERSE_FILE.exists():
        uni_df = pd.read_csv(UNIVERSE_FILE)
        uni_col = next((c for c in uni_df.columns if "gene" in c.lower() or "symbol" in c.lower()), None)
        M = int(uni_df[uni_col].dropna().nunique()) if uni_col else 20000
    else:
        M = 20000

    n = len(schema_genes)
    N = len(our_top_genes)
    p_val = float(hypergeom.sf(k - 1, M, n, N)) if k > 0 else 1.0

    print(f"SCHEMA genes (n): {n}")
    print(f"Top AlphaGenome genes (N): {N}")
    print(f"Overlap (k): {k}")
    print(f"Hypergeometric p-value: {p_val:.4e}")

    stats_file = OUT_DIR / "convergence_stats.txt"
    with open(stats_file, "w") as f:
        f.write(f"Universe_M={M}\n")
        f.write(f"SCHEMA_n={n}\n")
        f.write(f"TopAlphaGenome_N={N}\n")
        f.write(f"Overlap_k={k}\n")
        f.write(f"Hypergeometric_p={p_val}\n")
        f.write("Overlap_genes=" + ",".join(overlap) + "\n")

    pd.DataFrame({"gene": overlap}).to_csv(TABLE_OUT, index=False)

    print(f"Saved stats to {stats_file}")
    print(f"Saved overlap table to {TABLE_OUT}")


if __name__ == "__main__":
    main()
