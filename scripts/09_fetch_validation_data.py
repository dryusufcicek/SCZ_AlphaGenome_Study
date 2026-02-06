
import os
import subprocess
from pathlib import Path

# Config
DATA_DIR = Path("AlphaGenome_SCZ/data/external")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Files to download
URLS = {
    "Corces_Peaks_tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147672/suppl/GSE147672_scATAC_idr_peaks.tar.gz",
    "Trevino_Peaks_bed": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162170/suppl/GSE162170_atac_consensus_peaks.bed.gz"
}

def download_file(url, dest_path):
    print(f"Downloading {url}...")
    subprocess.run(["curl", "-L", "-o", str(dest_path), url], check=True)
    print(f"Saved to {dest_path}")

def main():
    print("Fetching Validation Data...")
    
    # 1. Corces (Adult)
    corces_dest = DATA_DIR / "GSE147672_scATAC_idr_peaks.tar.gz"
    if not corces_dest.exists():
        download_file(URLS["Corces_Peaks_tar"], corces_dest)
    
    # Extract Corces
    print("Extracting Corces...")
    subprocess.run(["tar", "-xzf", str(corces_dest), "-C", str(DATA_DIR)], check=True)
    
    # 2. Trevino (Fetal)
    trevino_dest = DATA_DIR / "GSE162170_atac_consensus_peaks.bed.gz"
    if not trevino_dest.exists():
        download_file(URLS["Trevino_Peaks_bed"], trevino_dest)
        
    print("\nDownload Complete. Listing files:")
    subprocess.run(["ls", "-R", str(DATA_DIR)])

if __name__ == "__main__":
    main()
