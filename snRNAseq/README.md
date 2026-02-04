# Create a small random demo dataset (snRNA-seq–like)

Because the original sequencing data are not publicly available, we provide instructions to generate a small synthetic dataset for testing the pipeline.
This dataset mimics sparse single-nucleus RNA-seq (snRNA-seq) count data and allows the full workflow to run quickly on a standard desktop computer.

Run the following R code to generate a demo Seurat object:
```bash
library(Seurat)
set.seed(1)
counts <- matrix(rpois(1000 * 400, lambda = 5), nrow = 1000, ncol = 400)
rownames(counts) <- paste0("Gene", 1:1000)
colnames(counts) <- paste0("Cell", 1:400)
obj <- CreateSeuratObject(counts = counts)
saveRDS(obj, "demo_seurat.rds")
```

This will create a small Seurat object (demo_seurat.rds) that can be used to run the complete pipeline in under few minutes.

- QC, Clustering and Annotation: 'snRNAseq_dlPFC_dataset.R' or 'snRNAseq_CBC_dataset.R'
- Subclustering and Trajectory analysis: 'Subclustering_and_Trajectory.R'
- Pseudobulk DEG analysis: 'Pseudobulk_DEG_EdgeR.R'
- Enrichment_analysis: 'Enrichment_analysis_cluster.R'
- Cell-cell interaction: 'CellChat.R'
- Integration analysis with Public data: 'Integration_analysis_Public_data.R'
