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

### Software
- Trim Galore (v0.6.6)
- FastQC (v0.11.9)
- MultiQC (v1.0.dev0)
- TopHat (v2.1.1)
- Cufflinks (v2.2.1)
- Cell Ranger (v5.0.1)
- Cytoscape (v3.8.0)
- GraphPad Prism 7 (v7.0.4)

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
Human reference (GRCh38) for RNA-seq
```bash
curl -O ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
Human reference (GRCh38) for snRNA-seq
```bash
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```
