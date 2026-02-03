# MDD_RNA_snRNA
This repository contains analysis pipelines for bulk RNA-seq and single-nucleus RNA-seq (snRNA-seq) data.

## System requirements
- R = 4.3.3 for Integration_analysis_Public_data.R
- R = 4.2.1 for Subclustering_and_Trajectory.R
- R = 4.0.2 for all other analyses
- Ubuntu 18.04 LTS
- Software used for data analysis: Trim Galore (v0.6.6), FastQC (v0.11.9), MultiQC (v1.0.dev0), TopHat (v2.1.1), Cufflinks (v2.2.1), R (v4.0.2, v4.2.1, v4.3.3), GraphPad Prism 7 (v7.0.4), Cytoscape (v3.8.0), Cell Ranger (v5.0.1).
- R packages: preprocessCore (v1.52.1), car (v3.1-2), limma (v3.46.0), WGCNA (v1.71), circlize (v0.4.15), Seurat (v4.3.0.1, v5.0.3), DoubletFinder (v2.0.3, v2.0.4), edgeR (v3.32.1, v4.0.16), monocle3 (v1.3.1, v1.3.7), CellChat (v1.6.1).

## Installation
Human reference (GRCh38) for RNA-seq
```bash
curl -O ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
Human reference (GRCh38) for snRNA-seq
```bash
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```
