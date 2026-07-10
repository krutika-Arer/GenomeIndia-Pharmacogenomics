# Indian Pharmacogenomics Analysis Pipeline using Genome India Project Data

---
# Overview

This repository contains a comprehensive bioinformatics pipeline developed for large-scale pharmacogenomic variant analysis using Genome India Project whole-genome sequencing data.

The pipeline integrates multiple genomic resources to identify clinically relevant pharmacogenomic variants enriched in the Indian population. Variant annotation is performed using Ensembl Variant Effect Predictor (VEP), followed by integration with ClinVar, PharmGKB, gnomAD, and 1000 Genomes databases.

The primary objective is to identify population-specific variants that may influence drug response, efficacy, toxicity, and precision medicine in Indian populations.

---

# Objectives

- Process Genome India Project variant datasets
- Annotate variants using Ensembl VEP
- Identify clinically significant variants
- Integrate ClinVar clinical annotations
- Integrate PharmGKB pharmacogenomic annotations
- Compare Indian allele frequencies with global populations
- Identify rare and Indian-enriched pharmacogenomic variants
- Generate analysis-ready datasets for downstream precision medicine research

---

# Workflow

```text
Genome India VCF Files
          │
          ▼
Variant Filtering
          │
          ▼
Coordinate Extraction
          │
          ▼
Ensembl VEP Annotation
          │
          ▼
ClinVar Annotation
          │
          ▼
PharmGKB Annotation
          │
          ▼
Population Frequency Analysis
          │
          ▼
Variant Classification
          │
          ▼
Final Pharmacogenomic Dataset
```

---

# Project Structure

```
Indian-PGX-GenomeIndia-Pipeline/

│
├── scripts/
│   ├── 01_merge_variants.py
│   ├── 02_clean_variants.py
│   ├── 03_extract_coordinates.py
│   ├── 04_run_vep.sh
│   ├── 05_parse_vep.py
│   ├── 06_clinvar_merge.py
│   ├── 07_pharmgkb_merge.py
│   ├── 08_frequency_analysis.py
│   └── plots.py
│
├── data/
│   ├── sample_data/
│   └── README.md
│
├── output/
│
├── figures/
│
├── docs/
│
├── environment.yml
├── requirements.txt
├── LICENSE
└── README.md
```

---

# Datasets Used

## Genome India Project

Whole Genome Sequencing variants from the Genome India Project.

Contains allele frequencies from Indian populations.

---

## Ensembl Variant Effect Predictor (VEP)

Version:
- Ensembl VEP 116
- GRCh38

Used for

- Variant consequences
- Gene annotation
- HGVS nomenclature
- SIFT
- PolyPhen
- Canonical transcripts

---

## ClinVar

Clinical significance annotation including

- Pathogenic
- Likely Pathogenic
- Benign
- Drug Response
- Risk Factor
- Protective variants

---

## PharmGKB

Used for

- Drug-Gene relationships
- Clinical annotations
- Pharmacogenomic evidence
- Variant-drug associations

---

## gnomAD

Population allele frequencies

Including

- AFR
- AMR
- ASJ
- EAS
- FIN
- NFE
- SAS

---

## 1000 Genomes

Used for comparison of Indian variants against global populations.

---

# Software Requirements

- Python 3.10+
- Bash
- Ensembl VEP
- samtools
- bcftools
- tabix
- pandas
- numpy
- matplotlib

---

# Installation

Clone repository

```bash
git clone https://github.com/yourusername/Indian-PGX-GenomeIndia-Pipeline.git

cd Indian-PGX-GenomeIndia-Pipeline
```

Create environment

```bash
conda env create -f environment.yml

conda activate pgx
```

---

# Running the Pipeline

Merge variants

```bash
python scripts/01_merge_variants.py
```

Extract coordinates

```bash
python scripts/03_extract_coordinates.py
```

Run VEP

```bash
bash scripts/04_run_vep.sh
```

Merge ClinVar

```bash
python scripts/06_clinvar_merge.py
```

Merge PharmGKB

```bash
python scripts/07_pharmgkb_merge.py
```

Population frequency analysis

```bash
python scripts/08_frequency_analysis.py
```

---

# Pipeline Output

The pipeline generates

- Annotated variants
- ClinVar integrated dataset
- PharmGKB integrated dataset
- Rare variant dataset
- Indian enriched variants
- Population allele frequency tables
- Drug-response variants
- Publication-ready tables

---

# Variant Annotation

Each variant is annotated with

- Chromosome
- Position
- Reference allele
- Alternate allele
- Gene
- Transcript
- Consequence
- Variant class
- Canonical transcript
- HGVS
- Protein change
- SIFT
- PolyPhen
- ClinVar significance
- PharmGKB evidence
- Global allele frequency
- Indian allele frequency

---

# Example Pipeline

```
VCF

↓

Quality Filtering

↓

Coordinate Extraction

↓

Variant Annotation (VEP)

↓

Clinical Annotation (ClinVar)

↓

Drug Annotation (PharmGKB)

↓

Population Comparison

↓

Final PGx Dataset
```

---

# Results

The pipeline enables identification of

- Pathogenic pharmacogenomic variants
- Drug-response variants
- Indian-specific alleles
- Rare global variants
- Clinically actionable variants
- Precision medicine biomarkers

---

# Applications

- Precision Medicine
- Clinical Genomics
- Population Genetics
- Pharmacogenomics
- Drug Discovery
- Personalized Medicine
- Translational Research

---

# Reproducibility

The analysis is fully reproducible using

- Python scripts
- Bash workflows
- Ensembl VEP
- Public annotation databases

---

# Citation

If you use this repository, please cite

Genome India Project

Ensembl Variant Effect Predictor (VEP)

ClinVar

PharmGKB

gnomAD

1000 Genomes Project

---

# Author

**Krutika Arer**

M.Sc. Bioinformatics

Shri Shankara Cancer Hospital and Research Centre

Department of Molecular Oncology

Bengaluru, Karnataka, India

---

# Acknowledgements

Genome India Project

Ensembl

ClinVar

PharmGKB

gnomAD Consortium

1000 Genomes Consortium

Open-source bioinformatics community

---

# License

This project is licensed under the MIT License.

---

## Future Work

- Machine learning prediction of pharmacogenomic variants
- Clinical decision support integration
- Structural variant annotation
- Copy number variation analysis
- AI-assisted variant prioritization
- Interactive web dashboard
- Population-specific pharmacogenomic database
