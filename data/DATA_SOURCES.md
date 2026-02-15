# Data Sources

## Primary Data Sources

All data used in this analysis are from publicly available sources.

### 1. Schizophrenia GWAS Fine-Mapping
**Source**: Psychiatric Genomics Consortium (PGC3)
**Citation**: Trubetskoy, V. et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature 604, 502–508 (2022).
**Data**: Fine-mapped credible sets (95% credible intervals)
**Access**: Supplementary Table 11d from Nature publication
**URL**: https://www.nature.com/articles/s41586-022-04434-5
**License**: Publicly available for research use

**Data specs**:
- 20,760 variants across 255 independent loci
- Posterior probabilities from FINEMAP
- Coordinates: GRCh38/hg38

---

### 2. Brain Hi-C Chromatin Loops
**Source**: PsychENCODE Consortium
**Citation**: Won, H. et al. Chromosome conformation elucidates regulatory relationships in developing human brain. Nature 538, 523–527 (2016).
**Data**: Fetal brain Hi-C chromatin loops
**Access**: PsychENCODE data portal
**URL**: http://psychencode.org/
**License**: Publicly available for research use

**Data specs**:
- 149,097 significant chromatin loops (FDR<0.05)
- Tissue: Fetal brain (germinal zone)
- Gestational age: 17-18 weeks
- Format: BEDPE (paired genomic coordinates)
- Resolution: ~40 kb

---

### 3. Brain Expression QTL (eQTL)
**Source**: Genotype-Tissue Expression (GTEx) Consortium v10
**Citation**: Aguet, F. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318-1330 (2020).
**Data**: Brain tissue eQTL (13 tissues)
**Access**: GTEx Portal
**URL**: https://gtexportal.org/
**License**: Open access (dbGaP for individual-level data)

**Brain tissues used**:
1. Amygdala
2. Anterior cingulate cortex (BA24)
3. Caudate (basal ganglia)
4. Cerebellar hemisphere
5. Cerebellum
6. Cortex
7. Frontal cortex (BA9)
8. Hippocampus
9. Hypothalamus
10. Nucleus accumbens (basal ganglia)
11. Putamen (basal ganglia)
12. Spinal cord (cervical c-1)
13. Substantia nigra

**Data specs**:
- Significant eGenes (q<0.05 by permutation)
- GTEx v10 release
- Sample sizes: 114-209 donors per tissue

---

### 4. Gene Annotations
**Source**: GENCODE
**Version**: v38 (GRCh38.p13)
**Citation**: Frankish, A. et al. GENCODE reference annotation for the human genome. Nucleic Acids Research 47, D766–D773 (2019).
**Access**: GENCODE portal
**URL**: https://www.gencodegenes.org/
**License**: Open access

**Data specs**:
- 19,966 protein-coding genes analyzed
- Coordinates: GRCh38/hg38
- Includes: Gene IDs, symbols, TSS positions

---

### 5. Gene Ontology Pathways
**Source**: Gene Ontology Consortium
**Release**: 2023
**Citation**: Gene Ontology Consortium. The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Research 49, D325–D334 (2021).
**Access**: Gene Ontology portal
**URL**: http://geneontology.org/
**License**: CC BY 4.0

**Data specs**:
- Biological Process ontology
- 3,623 pathways tested (filtered to 10-500 genes)
- GO annotations from GENCODE v38

---

### 6. AlphaGenome Regulatory Predictions
**Tool**: AlphaGenome
**Method**: Enformer-based deep learning
**Access**: API (https://alphagenome.ai)
**Citation**: [Tool citation - add if published]

**Prediction modalities**:
1. DNase-seq (chromatin accessibility)
2. ATAC-seq (chromatin accessibility)
3. CAGE (transcription start sites)
4. RNA-seq (steady-state RNA levels)
5. PRO-CAP (active transcription)
6. Splice junctions (alternative splicing)

**Specifications**:
- Sequence window: ±256 kb (512 kb total)
- Genome build: GRCh38/hg38
- Output: Log-fold change (reference vs alternative allele)

---

## Data Processing

### 1. Variant Filtering
- Started with PGC3 fine-mapped credible sets (95% CI)
- Retained all 20,760 variants across 255 loci
- No additional filtering applied
- Posterior probabilities normalized per locus

### 2. Gene Universe
- GENCODE v38 protein-coding genes
- Chromosomes 1-22, X only (excluded Y, MT, contigs)
- Total: 19,966 genes
- Tested: 1,617 genes with ≥1 variant connection

### 3. eQTL Data
- GTEx v10 significant eGenes (q<0.05)
- 13 brain tissues
- Used for validation only (not for gene filtering)

---

## Data Availability

### Included in This Repository
- ✅ Final analysis results (CSV files)
- ✅ Supplementary tables
- ✅ Analysis scripts

### Not Included (Available from Original Sources)
- ❌ Raw PGC3 GWAS summary statistics (see citation above)
- ❌ Raw Hi-C loop data (see PsychENCODE portal)
- ❌ Raw GTEx eQTL data (see GTEx portal)
- ❌ AlphaGenome predictions (generate via API)

### Intermediate Files
Intermediate analysis files are available upon request.

---

## Reproducibility

To reproduce the analysis:

1. Download PGC3 fine-mapped variants (Supplementary Table 11d)
2. Download PsychENCODE Hi-C loops
3. Download GTEx v10 brain eGenes
4. Download GENCODE v38 annotations
5. Run AlphaGenome predictions via API
6. Execute analysis pipeline (see `scripts/README.md`)

---

## Citations

If using these data sources, please cite the original publications:

```
Trubetskoy et al., Nature 2022 (PGC3 GWAS)
Won et al., Nature 2016 (PsychENCODE Hi-C)
Aguet et al., Science 2020 (GTEx v10)
Frankish et al., NAR 2019 (GENCODE)
GO Consortium, NAR 2021 (Gene Ontology)
```

---

**Last updated**: February 2025
**Contact**: Yusuf Cicek (yusuf.cicek@iuc.edu.tr)

