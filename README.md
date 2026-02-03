# MDD_RNA_snRNA

This repository contains analysis pipelines for bulk RNA-seq and single-nucleus RNA-seq (snRNA-seq) data.

---

## System requirements

### OS
- Linux (Ubuntu 18.04 LTS)
- Windows (R-based downstream analyses only)

### R versions
- R 4.3.3 → Integration_analysis_Public_data.R
- R 4.2.1 → Subclustering_and_Trajectory.R
- R 4.0.2 → all other analyses

### R packages
- preprocessCore (1.52.1)
- limma (3.46.0)
- edgeR (3.32.1, 4.0.16)
- WGCNA (1.71)
- Seurat (4.3.0.1, 5.0.3)
- DoubletFinder (2.0.3, 2.0.4)
- monocle3 (1.3.1, 1.3.7)
- CellChat (1.6.1)
- circlize (0.4.15)
- car (3.1-2)

---

## Installation
Human reference (GRCh38) and annotation (Ensembl release 100) for RNA-seq
```bash
curl -O ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
```
which should install in about 2-5 minutes.

Human reference (GRCh38) for snRNA-seq
```bash
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```
which should install in about 3-5 minutes.
