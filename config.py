"""
Configuration file for SCZ AlphaGenome Analysis Pipeline

This file contains paths and API configuration.
Users should set their AlphaGenome API key as an environment variable
or directly in this file (not recommended for public repositories).
"""

import os
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

# Base directory (adjust if running from different location)
BASE_DIR = Path(__file__).parent.resolve()

# Data directories
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Results directories
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

# Create directories if they don't exist
for dir_path in [RAW_DIR, PROCESSED_DIR, FIGURES_DIR, TABLES_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# =============================================================================
# API CONFIGURATION
# =============================================================================

# AlphaGenome API Key
# Set via environment variable: export ALPHAGENOME_API_KEY="your_key"
ALPHAGENOME_API_KEY = os.environ.get("ALPHAGENOME_API_KEY", None)

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# GWAS significance threshold
GWAS_PVALUE_THRESHOLD = 5e-8

# AlphaGenome scoring window (bases on each side of variant)
SCORING_WINDOW_HALF = 262144  # 256kb each side = 512kb total

# Brain tissues for filtering (GTEx nomenclature)
BRAIN_GTEX_TISSUES = [
    'Brain_Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia',
    'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum',
    'Brain_Cortex',
    'Brain_Frontal_Cortex_BA9',  # DLPFC
    'Brain_Hippocampus',
    'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c-1',
    'Brain_Substantia_nigra'
]

# Gene score threshold for "high-impact" classification
HIGH_IMPACT_SCORE_THRESHOLD = 1.0

# =============================================================================
# PATHWAY DEFINITIONS
# =============================================================================

PATHWAY_SETS = {
    'Calcium_VoltageGated': ['CACNA1A', 'CACNA1B', 'CACNA1C', 'CACNA1D', 'CACNA1E', 
                            'CACNA1F', 'CACNA1G', 'CACNA1H', 'CACNA1I', 'CACNA1S'],
    'Calcium_Internal': ['ATP2A2', 'ATP2A1', 'ATP2A3', 'RYR1', 'RYR2', 'RYR3', 
                        'ITPR1', 'ITPR2', 'ITPR3'],
    'Glutamate_Ionotropic': ['GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D', 
                            'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4', 
                            'GRIK1', 'GRIK2', 'GRIK3', 'GRIK4', 'GRIK5'],
    'Glutamate_Metabotropic': ['GRM1', 'GRM2', 'GRM3', 'GRM4', 'GRM5', 'GRM6', 'GRM7', 'GRM8'],
    'GABA_Receptors': ['GABRA1', 'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 
                      'GABRB1', 'GABRB2', 'GABRB3', 'GABRG1', 'GABRG2', 
                      'GABBR1', 'GABBR2'],
    'Dopamine_System': ['DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5', 'SLC6A3', 'TH', 'DDC'],
    'Synaptic_Scaffold': ['DLG1', 'DLG2', 'DLG3', 'DLG4', 'SHANK1', 'SHANK2', 'SHANK3',
                         'HOMER1', 'HOMER2', 'HOMER3'],
}
