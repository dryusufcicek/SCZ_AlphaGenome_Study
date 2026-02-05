"""
SCZ Hypothesis Testing - Configuration and Utilities
"""
import os
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
GWAS_DIR = DATA_DIR / "gwas"
BRAINSPAN_DIR = DATA_DIR / "brainspan"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# Ensure directories exist
for dir_path in [PROCESSED_DIR, RESULTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# API Keys
ALPHAGENOME_API_KEY = os.getenv("ALPHAGENOME_API_KEY")

# BrainSpan age bins (per SCZ Consensus paper - Pergola et al. 2023)
AGE_BINS = {
    'perinatal': (-0.7, 6),    # Fetal to 6 years
    'juvenile': (6, 25),       # 6-25 years (includes adolescence) ← KEY PERIOD
    'adult': (25, 50),         # 25-50 years
    'older_adult': (50, 100)   # 50+ years
}


# PGC3 genome-wide significance threshold
GWAS_PVALUE_THRESHOLD = 5e-8

def pcw_to_years(pcw):
    """Convert post-conception weeks to years (negative for prenatal)."""
    # 40 weeks = 0 years (birth)
    return (pcw - 40) / 52

def parse_brainspan_age(age_str):
    """Parse BrainSpan age string to years.
    
    Examples: '8 pcw', '4 mos', '2 yrs'
    """
    parts = age_str.lower().strip().split()
    if len(parts) != 2:
        return None
    
    value = float(parts[0])
    unit = parts[1]
    
    if 'pcw' in unit:
        return pcw_to_years(value)
    elif 'mos' in unit:
        return value / 12
    elif 'yrs' in unit or 'yr' in unit:
        return value
    else:
        return None

def get_age_bin(age_years):
    """Assign age in years to a developmental bin."""
    for bin_name, (low, high) in AGE_BINS.items():
        if low <= age_years < high:
            return bin_name
    return 'adult'  # Default for very old ages
