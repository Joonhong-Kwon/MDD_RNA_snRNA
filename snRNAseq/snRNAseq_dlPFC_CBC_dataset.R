library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(WGCNA)
library(ggplot2)
library(patchwork)
library(ape)
library(dittoSeq)
library(fields)
library(ggdendro)
library(metap)
suppressMessages(require(DoubletFinder))

### Load data
k1 <- Read10X(data.dir = "kang1/outs/filtered_feature_bc_matrix/")
k2 <- Read10X(data.dir = "kang2/outs/filtered_feature_bc_matrix/")
k3 <- Read10X(data.dir = "kang3/outs/filtered_feature_bc_matrix/")
k4 <- Read10X(data.dir = "kang4/outs/filtered_feature_bc_matrix/")
k17 <- Read10X(data.dir = "kang17/outs/k17_filtered_feature_bc_matrix/")
k18 <- Read10X(data.dir = "kang18/outs/k18_filtered_feature_bc_matrix/")
k19 <- Read10X(data.dir = "kang19/outs/filtered_feature_bc_matrix/")
k20 <- Read10X(data.dir = "kang20/outs/k20_filtered_feature_bc_matrix/")
k23 <- Read10X(data.dir = "kang23/outs/filtered_feature_bc_matrix/")
k24 <- Read10X(data.dir = "kang24/outs/filtered_feature_bc_matrix/")
k29 <- Read10X(data.dir = "kang29/outs/filtered_feature_bc_matrix/")
k30 <- Read10X(data.dir = "kang30/outs/filtered_feature_bc_matrix/")
k33 <- Read10X(data.dir = "kang33/outs/filtered_feature_bc_matrix/")
k34 <- Read10X(data.dir = "kang34/outs/filtered_feature_bc_matrix/")
k35 <- Read10X(data.dir = "kang35/outs/filtered_feature_bc_matrix/")
k36 <- Read10X(data.dir = "kang36/outs/filtered_feature_bc_matrix/")
k37 <- Read10X(data.dir = "kang37/outs/filtered_feature_bc_matrix/")
k38 <- Read10X(data.dir = "kang38/outs/filtered_feature_bc_matrix/")
k39 <- Read10X(data.dir = "kang39/outs/filtered_feature_bc_matrix/")
k40 <- Read10X(data.dir = "kang40/outs/filtered_feature_bc_matrix/")
k41 <- Read10X(data.dir = "kang41/outs/filtered_feature_bc_matrix/")
k42 <- Read10X(data.dir = "kang42/outs/filtered_feature_bc_matrix/")
k43 <- Read10X(data.dir = "kang43/outs/filtered_feature_bc_matrix/")
k44 <- Read10X(data.dir = "kang44/outs/filtered_feature_bc_matrix/")

### Create Seurat object
dlPFC_C1 <- CreateSeuratObject(counts=k19, min.cells=3, project="dlPFC_C1")
dlPFC_C2 <- CreateSeuratObject(counts=k3, min.cells=3, project="dlPFC_C2")
dlPFC_C3 <- CreateSeuratObject(counts=k43, min.cells=3, project="dlPFC_C3")
dlPFC_C4 <- CreateSeuratObject(counts=k23, min.cells=3, project="dlPFC_C4")
dlPFC_C5 <- CreateSeuratObject(counts=k35, min.cells=3, project="dlPFC_C5")
dlPFC_C6 <- CreateSeuratObject(counts=k39, min.cells=3, project="dlPFC_C6")

dlPFC_M1 <- CreateSeuratObject(counts=k17, min.cells=3, project="dlPFC_M1")
dlPFC_M2 <- CreateSeuratObject(counts=k1, min.cells=3, project="dlPFC_M2")
dlPFC_M3 <- CreateSeuratObject(counts=k41, min.cells=3, project="dlPFC_M3")
dlPFC_M4 <- CreateSeuratObject(counts=k29, min.cells=3, project="dlPFC_M4")
dlPFC_M5 <- CreateSeuratObject(counts=k33, min.cells=3, project="dlPFC_M5")
dlPFC_M6 <- CreateSeuratObject(counts=k37, min.cells=3, project="dlPFC_M6")

CBC_C1 <- CreateSeuratObject(counts=k20, min.cells=3, project="CBC_C1")
CBC_C2 <- CreateSeuratObject(counts=k4, min.cells=3, project="CBC_C2")
CBC_C3 <- CreateSeuratObject(counts=k44, min.cells=3, project="CBC_C3")
CBC_C4 <- CreateSeuratObject(counts=k24, min.cells=3, project="CBC_C4")
CBC_C5 <- CreateSeuratObject(counts=k36, min.cells=3, project="CBC_C5")
CBC_C6 <- CreateSeuratObject(counts=k40, min.cells=3, project="CBC_C6")

CBC_M1 <- CreateSeuratObject(counts=k18, min.cells=3, project="CBC_M1")
CBC_M2 <- CreateSeuratObject(counts=k2, min.cells=3, project="CBC_M2")
CBC_M3 <- CreateSeuratObject(counts=k42, min.cells=3, project="CBC_M3")
CBC_M4 <- CreateSeuratObject(counts=k30, min.cells=3, project="CBC_M4")
CBC_M5 <- CreateSeuratObject(counts=k34, min.cells=3, project="CBC_M5")
CBC_M6 <- CreateSeuratObject(counts=k38, min.cells=3, project="CBC_M6")

### remove matrix
rm (k1, k2, k3, k4, k17, k18, k19, k20, k23, k24, k29, k30, k33, k34, k35, k36, k37, k38, k39, k40,
    k41, k42, k43, k44)

### QC metrics and each sample cell filtering
### Mitochondrial genes
mito_genes <- rownames(dlPFC_C1)[grep("^MT-", rownames(dlPFC_C1))]
dlPFC_C1[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C1, pattern="^MT-")
dlPFC_C2[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C2, pattern="^MT-")
dlPFC_C3[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C3, pattern="^MT-")
dlPFC_C4[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C4, pattern="^MT-")
dlPFC_C5[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C5, pattern="^MT-")
dlPFC_C6[["percent.mt"]] <- PercentageFeatureSet(dlPFC_C6, pattern="^MT-")

dlPFC_M1[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M1, pattern="^MT-")
dlPFC_M2[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M2, pattern="^MT-")
dlPFC_M3[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M3, pattern="^MT-")
dlPFC_M4[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M4, pattern="^MT-")
dlPFC_M5[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M5, pattern="^MT-")
dlPFC_M6[["percent.mt"]] <- PercentageFeatureSet(dlPFC_M6, pattern="^MT-")

CBC_C1[["percent.mt"]] <- PercentageFeatureSet(CBC_C1, pattern="^MT-")
CBC_C2[["percent.mt"]] <- PercentageFeatureSet(CBC_C2, pattern="^MT-")
CBC_C3[["percent.mt"]] <- PercentageFeatureSet(CBC_C3, pattern="^MT-")
CBC_C4[["percent.mt"]] <- PercentageFeatureSet(CBC_C4, pattern="^MT-")
CBC_C5[["percent.mt"]] <- PercentageFeatureSet(CBC_C5, pattern="^MT-")
CBC_C6[["percent.mt"]] <- PercentageFeatureSet(CBC_C6, pattern="^MT-")

CBC_M1[["percent.mt"]] <- PercentageFeatureSet(CBC_M1, pattern="^MT-")
CBC_M2[["percent.mt"]] <- PercentageFeatureSet(CBC_M2, pattern="^MT-")
CBC_M3[["percent.mt"]] <- PercentageFeatureSet(CBC_M3, pattern="^MT-")
CBC_M4[["percent.mt"]] <- PercentageFeatureSet(CBC_M4, pattern="^MT-")
CBC_M5[["percent.mt"]] <- PercentageFeatureSet(CBC_M5, pattern="^MT-")
CBC_M6[["percent.mt"]] <- PercentageFeatureSet(CBC_M6, pattern="^MT-")

### Ribosomal genes
ribo_genes <- rownames(dlPFC_C1)[grep("^RP[SL]", rownames(dlPFC_C1))]
dlPFC_C1[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C1, pattern="^RP[SL]")
dlPFC_C2[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C2, pattern="^RP[SL]")
dlPFC_C3[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C3, pattern="^RP[SL]")
dlPFC_C4[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C4, pattern="^RP[SL]")
dlPFC_C5[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C5, pattern="^RP[SL]")
dlPFC_C6[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_C6, pattern="^RP[SL]")

dlPFC_M1[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M1, pattern="^RP[SL]")
dlPFC_M2[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M2, pattern="^RP[SL]")
dlPFC_M3[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M3, pattern="^RP[SL]")
dlPFC_M4[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M4, pattern="^RP[SL]")
dlPFC_M5[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M5, pattern="^RP[SL]")
dlPFC_M6[["percent.ribo"]] <- PercentageFeatureSet(dlPFC_M6, pattern="^RP[SL]")

CBC_C1[["percent.ribo"]] <- PercentageFeatureSet(CBC_C1, pattern="^RP[SL]")
CBC_C2[["percent.ribo"]] <- PercentageFeatureSet(CBC_C2, pattern="^RP[SL]")
CBC_C3[["percent.ribo"]] <- PercentageFeatureSet(CBC_C3, pattern="^RP[SL]")
CBC_C4[["percent.ribo"]] <- PercentageFeatureSet(CBC_C4, pattern="^RP[SL]")
CBC_C5[["percent.ribo"]] <- PercentageFeatureSet(CBC_C5, pattern="^RP[SL]")
CBC_C6[["percent.ribo"]] <- PercentageFeatureSet(CBC_C6, pattern="^RP[SL]")

CBC_M1[["percent.ribo"]] <- PercentageFeatureSet(CBC_M1, pattern="^RP[SL]")
CBC_M2[["percent.ribo"]] <- PercentageFeatureSet(CBC_M2, pattern="^RP[SL]")
CBC_M3[["percent.ribo"]] <- PercentageFeatureSet(CBC_M3, pattern="^RP[SL]")
CBC_M4[["percent.ribo"]] <- PercentageFeatureSet(CBC_M4, pattern="^RP[SL]")
CBC_M5[["percent.ribo"]] <- PercentageFeatureSet(CBC_M5, pattern="^RP[SL]")
CBC_M6[["percent.ribo"]] <- PercentageFeatureSet(CBC_M6, pattern="^RP[SL]")

### Hemoglobin genes
hemo_genes <- rownames(dlPFC_C1)[grep("^HB[^(P)]", rownames(dlPFC_C1))]
dlPFC_C1[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C1, pattern="^HB[^(P)]")
dlPFC_C2[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C2, pattern="^HB[^(P)]")
dlPFC_C3[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C3, pattern="^HB[^(P)]")
dlPFC_C4[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C4, pattern="^HB[^(P)]")
dlPFC_C5[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C5, pattern="^HB[^(P)]")
dlPFC_C6[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_C6, pattern="^HB[^(P)]")

dlPFC_M1[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M1, pattern="^HB[^(P)]")
dlPFC_M2[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M2, pattern="^HB[^(P)]")
dlPFC_M3[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M3, pattern="^HB[^(P)]")
dlPFC_M4[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M4, pattern="^HB[^(P)]")
dlPFC_M5[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M5, pattern="^HB[^(P)]")
dlPFC_M6[["percent.hemo"]] <- PercentageFeatureSet(dlPFC_M6, pattern="^HB[^(P)]")

CBC_C1[["percent.hemo"]] <- PercentageFeatureSet(CBC_C1, pattern="^HB[^(P)]")
CBC_C2[["percent.hemo"]] <- PercentageFeatureSet(CBC_C2, pattern="^HB[^(P)]")
CBC_C3[["percent.hemo"]] <- PercentageFeatureSet(CBC_C3, pattern="^HB[^(P)]")
CBC_C4[["percent.hemo"]] <- PercentageFeatureSet(CBC_C4, pattern="^HB[^(P)]")
CBC_C5[["percent.hemo"]] <- PercentageFeatureSet(CBC_C5, pattern="^HB[^(P)]")
CBC_C6[["percent.hemo"]] <- PercentageFeatureSet(CBC_C6, pattern="^HB[^(P)]")

CBC_M1[["percent.hemo"]] <- PercentageFeatureSet(CBC_M1, pattern="^HB[^(P)]")
CBC_M2[["percent.hemo"]] <- PercentageFeatureSet(CBC_M2, pattern="^HB[^(P)]")
CBC_M3[["percent.hemo"]] <- PercentageFeatureSet(CBC_M3, pattern="^HB[^(P)]")
CBC_M4[["percent.hemo"]] <- PercentageFeatureSet(CBC_M4, pattern="^HB[^(P)]")
CBC_M5[["percent.hemo"]] <- PercentageFeatureSet(CBC_M5, pattern="^HB[^(P)]")
CBC_M6[["percent.hemo"]] <- PercentageFeatureSet(CBC_M6, pattern="^HB[^(P)]")

### Add cell id
dlPFC_C1 <- RenameCells(object=dlPFC_C1, add.cell.id="dlPFC_C1")
dlPFC_C2 <- RenameCells(object=dlPFC_C2, add.cell.id="dlPFC_C2")
dlPFC_C3 <- RenameCells(object=dlPFC_C3, add.cell.id="dlPFC_C3")
dlPFC_C4 <- RenameCells(object=dlPFC_C4, add.cell.id="dlPFC_C4")
dlPFC_C5 <- RenameCells(object=dlPFC_C5, add.cell.id="dlPFC_C5")
dlPFC_C6 <- RenameCells(object=dlPFC_C6, add.cell.id="dlPFC_C6")

dlPFC_M1 <- RenameCells(object=dlPFC_M1, add.cell.id="dlPFC_M1")
dlPFC_M2 <- RenameCells(object=dlPFC_M2, add.cell.id="dlPFC_M2")
dlPFC_M3 <- RenameCells(object=dlPFC_M3, add.cell.id="dlPFC_M3")
dlPFC_M4 <- RenameCells(object=dlPFC_M4, add.cell.id="dlPFC_M4")
dlPFC_M5 <- RenameCells(object=dlPFC_M5, add.cell.id="dlPFC_M5")
dlPFC_M6 <- RenameCells(object=dlPFC_M6, add.cell.id="dlPFC_M6")

CBC_C1 <- RenameCells(object=CBC_C1, add.cell.id="CBC_C1")
CBC_C2 <- RenameCells(object=CBC_C2, add.cell.id="CBC_C2")
CBC_C3 <- RenameCells(object=CBC_C3, add.cell.id="CBC_C3")
CBC_C4 <- RenameCells(object=CBC_C4, add.cell.id="CBC_C4")
CBC_C5 <- RenameCells(object=CBC_C5, add.cell.id="CBC_C5")
CBC_C6 <- RenameCells(object=CBC_C6, add.cell.id="CBC_C6")

CBC_M1 <- RenameCells(object=CBC_M1, add.cell.id="CBC_M1")
CBC_M2 <- RenameCells(object=CBC_M2, add.cell.id="CBC_M2")
CBC_M3 <- RenameCells(object=CBC_M3, add.cell.id="CBC_M3")
CBC_M4 <- RenameCells(object=CBC_M4, add.cell.id="CBC_M4")
CBC_M5 <- RenameCells(object=CBC_M5, add.cell.id="CBC_M5")
CBC_M6 <- RenameCells(object=CBC_M6, add.cell.id="CBC_M6")

### merge datasets for plot
alldata_before <- merge(dlPFC_C1, c(dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                             dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6,
                             CBC_C1, CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                             CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6))
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
VlnPlot(alldata_before, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()
FeatureScatter(alldata_before, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
FeatureScatter(alldata_before, "nCount_RNA", "percent.mt", group.by = "orig.ident", pt.size = 0.5)

### Histogram
metadata1 <- alldata_before@meta.data
metadata1$sample <- NA
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C1"))] <- "dlPFC_C1"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C2"))] <- "dlPFC_C2"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C3"))] <- "dlPFC_C3"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C4"))] <- "dlPFC_C4"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C5"))] <- "dlPFC_C5"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_C6"))] <- "dlPFC_C6"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M1"))] <- "dlPFC_M1"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M2"))] <- "dlPFC_M2"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M3"))] <- "dlPFC_M3"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M4"))] <- "dlPFC_M4"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M5"))] <- "dlPFC_M5"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^dlPFC_M6"))] <- "dlPFC_M6"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C1"))] <- "CBC_C1"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C2"))] <- "CBC_C2"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C3"))] <- "CBC_C3"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C4"))] <- "CBC_C4"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C5"))] <- "CBC_C5"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_C6"))] <- "CBC_C6"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M1"))] <- "CBC_M1"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M2"))] <- "CBC_M2"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M3"))] <- "CBC_M3"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M4"))] <- "CBC_M4"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M5"))] <- "CBC_M5"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^CBC_M6"))] <- "CBC_M6"
alldata_before@meta.data <- metadata1

metadata1$group <- NA
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
metadata1$group[which(str_detect(metadata1$orig.ident, "^CBC_C"))] <- "CBC_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^CBC_M"))] <- "CBC_M"
alldata_before@meta.data <- metadata1

### Visualize the number UMIs/transcripts per cell
metadata1 %>% 
  ggplot(aes(color=group, x=nCount_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 9000)

# Visualize the distribution of genes detected per cell via histogram
metadata1 %>% 
  ggplot(aes(color=group, x=nFeature_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(400))

VlnPlot(alldata_before, group.by = "group", features = feats, pt.size = 0, ncol = 5) +
  NoLegend()

# save
saveRDS(object=alldata_before, file="Alldata_before_QC.RDS")


### Check top 5% nCount by quantile function
#quantile(dlPFC_C2$nCount_RNA,probs = c(0.75, 0.95, 0.995))
dlPFC_C1 <- subset(x=dlPFC_C1, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_C2 <- subset(x=dlPFC_C2, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_C3 <- subset(x=dlPFC_C3, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_C4 <- subset(x=dlPFC_C4, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_C5 <- subset(x=dlPFC_C5, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_C6 <- subset(x=dlPFC_C6, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)

dlPFC_M1 <- subset(x=dlPFC_M1, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_M2 <- subset(x=dlPFC_M2, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_M3 <- subset(x=dlPFC_M3, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_M4 <- subset(x=dlPFC_M4, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_M5 <- subset(x=dlPFC_M5, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
dlPFC_M6 <- subset(x=dlPFC_M6, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)

CBC_C1 <- subset(x=CBC_C1, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_C2 <- subset(x=CBC_C2, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_C3 <- subset(x=CBC_C3, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_C4 <- subset(x=CBC_C4, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_C5 <- subset(x=CBC_C5, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_C6 <- subset(x=CBC_C6, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)

CBC_M1 <- subset(x=CBC_M1, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_M2 <- subset(x=CBC_M2, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_M3 <- subset(x=CBC_M3, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_M4 <- subset(x=CBC_M4, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_M5 <- subset(x=CBC_M5, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)
CBC_M6 <- subset(x=CBC_M6, subset=nFeature_RNA > 400 & nCount_RNA < 9000 & percent.mt < 5)

alldata_after1 <- merge(dlPFC_C1, c(dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                             dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6,
                             CBC_C1, CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                             CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6))

VlnPlot(alldata_after1, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()

### Histogram
metadata2 <- alldata_after1@meta.data
metadata2$sample <- NA
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C1"))] <- "dlPFC_C1"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C2"))] <- "dlPFC_C2"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C3"))] <- "dlPFC_C3"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C4"))] <- "dlPFC_C4"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C5"))] <- "dlPFC_C5"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_C6"))] <- "dlPFC_C6"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M1"))] <- "dlPFC_M1"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M2"))] <- "dlPFC_M2"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M3"))] <- "dlPFC_M3"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M4"))] <- "dlPFC_M4"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M5"))] <- "dlPFC_M5"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^dlPFC_M6"))] <- "dlPFC_M6"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C1"))] <- "CBC_C1"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C2"))] <- "CBC_C2"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C3"))] <- "CBC_C3"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C4"))] <- "CBC_C4"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C5"))] <- "CBC_C5"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_C6"))] <- "CBC_C6"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M1"))] <- "CBC_M1"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M2"))] <- "CBC_M2"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M3"))] <- "CBC_M3"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M4"))] <- "CBC_M4"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M5"))] <- "CBC_M5"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^CBC_M6"))] <- "CBC_M6"
alldata_after1@meta.data <- metadata2

metadata2$group <- NA
metadata2$group[which(str_detect(metadata2$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata2$group[which(str_detect(metadata2$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
metadata2$group[which(str_detect(metadata2$orig.ident, "^CBC_C"))] <- "CBC_C"
metadata2$group[which(str_detect(metadata2$orig.ident, "^CBC_M"))] <- "CBC_M"
alldata_after1@meta.data <- metadata2

### Visualize the number UMIs/transcripts per cell
metadata2 %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
  #geom_vline(xintercept = 10000)

# Visualize the distribution of genes detected per cell via histogram
metadata2 %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
  #geom_vline(xintercept = c(300, 400))


### Remove Mitochondrial gene (MT-)
mito_genes <- rownames(dlPFC_C1)[grep("^MT-", rownames(dlPFC_C1))]
mito_genes
dlPFC_C1 <- dlPFC_C1[!grepl("^MT-", rownames(dlPFC_C1)), ]
dlPFC_C2 <- dlPFC_C2[!grepl("^MT-", rownames(dlPFC_C2)), ]
dlPFC_C3 <- dlPFC_C3[!grepl("^MT-", rownames(dlPFC_C3)), ]
dlPFC_C4 <- dlPFC_C4[!grepl("^MT-", rownames(dlPFC_C4)), ]
dlPFC_C5 <- dlPFC_C5[!grepl("^MT-", rownames(dlPFC_C5)), ]
dlPFC_C6 <- dlPFC_C6[!grepl("^MT-", rownames(dlPFC_C6)), ]

dlPFC_M1 <- dlPFC_M1[!grepl("^MT-", rownames(dlPFC_M1)), ]
dlPFC_M2 <- dlPFC_M2[!grepl("^MT-", rownames(dlPFC_M2)), ]
dlPFC_M3 <- dlPFC_M3[!grepl("^MT-", rownames(dlPFC_M3)), ]
dlPFC_M4 <- dlPFC_M4[!grepl("^MT-", rownames(dlPFC_M4)), ]
dlPFC_M5 <- dlPFC_M5[!grepl("^MT-", rownames(dlPFC_M5)), ]
dlPFC_M6 <- dlPFC_M6[!grepl("^MT-", rownames(dlPFC_M6)), ]

CBC_C1 <- CBC_C1[!grepl("^MT-", rownames(CBC_C1)), ]
CBC_C2 <- CBC_C2[!grepl("^MT-", rownames(CBC_C2)), ]
CBC_C3 <- CBC_C3[!grepl("^MT-", rownames(CBC_C3)), ]
CBC_C4 <- CBC_C4[!grepl("^MT-", rownames(CBC_C4)), ]
CBC_C5 <- CBC_C5[!grepl("^MT-", rownames(CBC_C5)), ]
CBC_C6 <- CBC_C6[!grepl("^MT-", rownames(CBC_C6)), ]

CBC_M1 <- CBC_M1[!grepl("^MT-", rownames(CBC_M1)), ]
CBC_M2 <- CBC_M2[!grepl("^MT-", rownames(CBC_M2)), ]
CBC_M3 <- CBC_M3[!grepl("^MT-", rownames(CBC_M3)), ]
CBC_M4 <- CBC_M4[!grepl("^MT-", rownames(CBC_M4)), ]
CBC_M5 <- CBC_M5[!grepl("^MT-", rownames(CBC_M5)), ]
CBC_M6 <- CBC_M6[!grepl("^MT-", rownames(CBC_M6)), ]


### combine list
combined.list <- list(dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                      dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6,
                      CBC_C1, CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                      CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6)
names(combined.list) <- c("dlPFC_C1", "dlPFC_C2", "dlPFC_C3", "dlPFC_C4", "dlPFC_C5", "dlPFC_C6",
                          "dlPFC_M1", "dlPFC_M2", "dlPFC_M3", "dlPFC_M4", "dlPFC_M5", "dlPFC_M6",
                          "CBC_C1", "CBC_C2", "CBC_C3", "CBC_C4", "CBC_C5", "CBC_C6",
                          "CBC_M1", "CBC_M2", "CBC_M3", "CBC_M4", "CBC_M5", "CBC_M6")

### Normalize and identify variable features for each dataset independently
combined.list <- lapply(X = combined.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#v.genes <- VariableFeatures(object=dlPFC_C1)
#LabelPoints(plot=VariableFeaturePlot(dlPFC_C1), points=v.genes[1:10], repel=T)

### DoubletFinder
combined.list <- lapply(X = combined.list, FUN = function(x) {
 x <- ScaleData(x, vars.to.regress = "nCount_RNA", verbose = F)
 x <- RunPCA(x, verbose = F, npcs = 20)
 x <- FindNeighbors(x, reduction = "pca", dims=1:10, verbose=F)
 x <- FindClusters(x, verbose=F)
 x <- RunUMAP(x, dims=1:10, verbose=F)
 })

# save
saveRDS(object=combined.list, file="Seurat.before.DF.RDS")

sweep.res <- paramSweep_v3(combined.list$CBC_M6, PCs=1:10)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")


homotypic.prop <- modelHomotypic(combined.list$CBC_M6@meta.data$RNA_snn_res.0.8)
nExp <- round(ncol(combined.list$CBC_M6) * 0.08)  # expect 8% doublets
nExp.adj <- round(nExp*(1 - homotypic.prop))
combined.list$CBC_M6 <- doubletFinder_v3(combined.list$CBC_M6, pN = 0.25, pK = 0.29, nExp = nExp,
                                           PCs = 1:10, reuse.pANN=F)
head(combined.list$CBC_M6@meta.data)[9]
combined.list$CBC_M6 <- doubletFinder_v3(combined.list$CBC_M6, pN = 0.25, pK = 0.29, nExp = nExp.adj,
                                           PCs = 1:10, reuse.pANN="pANN_0.25_0.29_748")

# save
saveRDS(object=combined.list, file="Seurat.after.DF_Homotypic.RDS")

# Change column name
head(combined.list$dlPFC_C1@meta.data)
colnames(combined.list$dlPFC_C1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_C2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_C3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_C4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_C5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_C6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$dlPFC_M1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_M2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_M3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_M4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_M5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$dlPFC_M6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$CBC_C1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_C2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_C3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_C4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_C5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_C6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$CBC_M1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_M2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_M3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_M4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_M5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$CBC_M6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")


### remove doublet
combined.list$dlPFC_C1 <- subset(combined.list$dlPFC_C1, subset = DF_class_high == "Singlet")
combined.list$dlPFC_C2 <- subset(combined.list$dlPFC_C2, subset = DF_class_high == "Singlet")
combined.list$dlPFC_C3 <- subset(combined.list$dlPFC_C3, subset = DF_class_high == "Singlet")
combined.list$dlPFC_C4 <- subset(combined.list$dlPFC_C4, subset = DF_class_high == "Singlet")
combined.list$dlPFC_C5 <- subset(combined.list$dlPFC_C5, subset = DF_class_high == "Singlet")
combined.list$dlPFC_C6 <- subset(combined.list$dlPFC_C6, subset = DF_class_high == "Singlet")

combined.list$dlPFC_M1 <- subset(combined.list$dlPFC_M1, subset = DF_class_high == "Singlet")
combined.list$dlPFC_M2 <- subset(combined.list$dlPFC_M2, subset = DF_class_high == "Singlet")
combined.list$dlPFC_M3 <- subset(combined.list$dlPFC_M3, subset = DF_class_high == "Singlet")
combined.list$dlPFC_M4 <- subset(combined.list$dlPFC_M4, subset = DF_class_high == "Singlet")
combined.list$dlPFC_M5 <- subset(combined.list$dlPFC_M5, subset = DF_class_high == "Singlet")
combined.list$dlPFC_M6 <- subset(combined.list$dlPFC_M6, subset = DF_class_high == "Singlet")

combined.list$CBC_C1 <- subset(combined.list$CBC_C1, subset = DF_class_high == "Singlet")
combined.list$CBC_C2 <- subset(combined.list$CBC_C2, subset = DF_class_high == "Singlet")
combined.list$CBC_C3 <- subset(combined.list$CBC_C3, subset = DF_class_high == "Singlet")
combined.list$CBC_C4 <- subset(combined.list$CBC_C4, subset = DF_class_high == "Singlet")
combined.list$CBC_C5 <- subset(combined.list$CBC_C5, subset = DF_class_high == "Singlet")
combined.list$CBC_C6 <- subset(combined.list$CBC_C6, subset = DF_class_high == "Singlet")

combined.list$CBC_M1 <- subset(combined.list$CBC_M1, subset = DF_class_high == "Singlet")
combined.list$CBC_M2 <- subset(combined.list$CBC_M2, subset = DF_class_high == "Singlet")
combined.list$CBC_M3 <- subset(combined.list$CBC_M3, subset = DF_class_high == "Singlet")
combined.list$CBC_M4 <- subset(combined.list$CBC_M4, subset = DF_class_high == "Singlet")
combined.list$CBC_M5 <- subset(combined.list$CBC_M5, subset = DF_class_high == "Singlet")
combined.list$CBC_M6 <- subset(combined.list$CBC_M6, subset = DF_class_high == "Singlet")

# save
saveRDS(object=combined.list, file="Seurat.remove.DF_Homotypic.RDS")

quantile(combined.list$CBC_M6$nFeature_RNA,probs = c(0.5))
quantile(combined.list$CBC_M6$nCount_RNA,probs = c(0.5))

alldata_after <- merge(combined.list$dlPFC_C1, c(combined.list$dlPFC_C2, combined.list$dlPFC_C3, combined.list$dlPFC_C4, combined.list$dlPFC_C5, combined.list$dlPFC_C6,
                                                 combined.list$dlPFC_M1, combined.list$dlPFC_M2, combined.list$dlPFC_M3, combined.list$dlPFC_M4, combined.list$dlPFC_M5, combined.list$dlPFC_M6,
                                                 combined.list$CBC_C1, combined.list$CBC_C2, combined.list$CBC_C3, combined.list$CBC_C4, combined.list$CBC_C5, combined.list$CBC_C6,
                                                 combined.list$CBC_M1, combined.list$CBC_M2, combined.list$CBC_M3, combined.list$CBC_M4, combined.list$CBC_M5, combined.list$CBC_M6))

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
VlnPlot(alldata_after, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()

### Histogram
metadata3 <- alldata_after@meta.data
metadata3$sample <- NA
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C1"))] <- "dlPFC_C1"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C2"))] <- "dlPFC_C2"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C3"))] <- "dlPFC_C3"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C4"))] <- "dlPFC_C4"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C5"))] <- "dlPFC_C5"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_C6"))] <- "dlPFC_C6"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M1"))] <- "dlPFC_M1"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M2"))] <- "dlPFC_M2"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M3"))] <- "dlPFC_M3"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M4"))] <- "dlPFC_M4"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M5"))] <- "dlPFC_M5"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^dlPFC_M6"))] <- "dlPFC_M6"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C1"))] <- "CBC_C1"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C2"))] <- "CBC_C2"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C3"))] <- "CBC_C3"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C4"))] <- "CBC_C4"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C5"))] <- "CBC_C5"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_C6"))] <- "CBC_C6"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M1"))] <- "CBC_M1"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M2"))] <- "CBC_M2"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M3"))] <- "CBC_M3"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M4"))] <- "CBC_M4"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M5"))] <- "CBC_M5"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^CBC_M6"))] <- "CBC_M6"
alldata_after@meta.data <- metadata3

metadata3$group <- NA
metadata3$group[which(str_detect(metadata3$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata3$group[which(str_detect(metadata3$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
metadata3$group[which(str_detect(metadata3$orig.ident, "^CBC_C"))] <- "CBC_C"
metadata3$group[which(str_detect(metadata3$orig.ident, "^CBC_M"))] <- "CBC_M"
alldata_after@meta.data <- metadata3

### Visualize the number UMIs/transcripts per cell
metadata3 %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
#geom_vline(xintercept = 10000)

# Visualize the distribution of genes detected per cell via histogram
metadata3 %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
#geom_vline(xintercept = c(300, 400))

VlnPlot(alldata_after, group.by = "group", features = feats, pt.size = 0, ncol = 5) +
  NoLegend()

# save (alldata before, after)
saveRDS(object=alldata_after, file="Alldata_after_QC&DF.RDS")
#alldata_after <- readRDS("Alldata_after_QC&DF.RDS")

### select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = combined.list)

### Seurat Integration (too long)
combined.anchors <- FindIntegrationAnchors(object.list = combined.list, anchor.features = features)
full.combined <- IntegrateData(anchorset = combined.anchors)

View(full.combined@meta.data)

# save
saveRDS(object=full.combined, file="Seurat.Integration.RDS")


metadata <- full.combined@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C1_"))] <- "dlPFC_C1"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C2_"))] <- "dlPFC_C2"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C3_"))] <- "dlPFC_C3"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C4_"))] <- "dlPFC_C4"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C5_"))] <- "dlPFC_C5"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C6_"))] <- "dlPFC_C6"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M1_"))] <- "dlPFC_M1"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M2_"))] <- "dlPFC_M2"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M3_"))] <- "dlPFC_M3"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M4_"))] <- "dlPFC_M4"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M5_"))] <- "dlPFC_M5"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M6_"))] <- "dlPFC_M6"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C1_"))] <- "CBC_C1"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C2_"))] <- "CBC_C2"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C3_"))] <- "CBC_C3"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C4_"))] <- "CBC_C4"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C5_"))] <- "CBC_C5"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C6_"))] <- "CBC_C6"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M1_"))] <- "CBC_M1"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M2_"))] <- "CBC_M2"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M3_"))] <- "CBC_M3"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M4_"))] <- "CBC_M4"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M5_"))] <- "CBC_M5"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M6_"))] <- "CBC_M6"
full.combined@meta.data <- metadata

metadata$group <- NA
metadata$group[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "dlPFC_C"
metadata$group[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "dlPFC_M"
metadata$group[which(str_detect(metadata$cells, "^CBC_C"))] <- "CBC_C"
metadata$group[which(str_detect(metadata$cells, "^CBC_M"))] <- "CBC_M"
full.combined@meta.data <- metadata

metadata$region <- NA
metadata$region[which(str_detect(metadata$cells, "^dlPFC"))] <- "dlPFC"
metadata$region[which(str_detect(metadata$cells, "^CBC"))] <- "CBC"
full.combined@meta.data <- metadata

metadata$state <- NA
metadata$state[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "Con"
metadata$state[which(str_detect(metadata$cells, "^CBC_C"))] <- "Con"
metadata$state[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "MDD"
metadata$state[which(str_detect(metadata$cells, "^CBC_M"))] <- "MDD"
full.combined@meta.data <- metadata

metadata$case <- NA
metadata$case[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "Con"
metadata$case[which(str_detect(metadata$cells, "^CBC_C"))] <- "Con"
metadata$case[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "MDD"
metadata$case[which(str_detect(metadata$cells, "^CBC_M"))] <- "MDD"
full.combined@meta.data <- metadata

metadata$donor <- NA
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C1"))] <- "donor1"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C1"))] <- "donor1"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C2_"))] <- "donor2"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C2_"))] <- "donor2"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C3_"))] <- "donor3"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C3_"))] <- "donor3"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C4_"))] <- "donor4"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C4_"))] <- "donor4"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C5_"))] <- "donor4"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C5_"))] <- "donor5"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_C6_"))] <- "donor6"
metadata$donor[which(str_detect(metadata$cells, "^CBC_C6_"))] <- "donor6"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M1_"))] <- "donor7"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M1_"))] <- "donor7"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M2_"))] <- "donor8"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M2_"))] <- "donor8"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M3_"))] <- "donor9"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M3_"))] <- "donor9"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M4_"))] <- "donor10"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M4_"))] <- "donor10"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M5_"))] <- "donor11"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M5_"))] <- "donor11"
metadata$donor[which(str_detect(metadata$cells, "^dlPFC_M6_"))] <- "donor12"
metadata$donor[which(str_detect(metadata$cells, "^CBC_M6_"))] <- "donor12"
full.combined@meta.data <- metadata

View(full.combined@meta.data)

### Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")
  #geom_vline(xintercept = 30000)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
  #geom_vline(xintercept = c(300, 10000))

### 
DefaultAssay(full.combined) <- "integrated"

### Run the standard workflow for visualization and clustering
full.combined <- ScaleData(full.combined, verbose = FALSE, vars.to.regress="nCount_RNA")
full.combined <- RunPCA(full.combined, npcs = 100, verbose = FALSE)

ElbowPlot(object=full.combined, ndims=100)

saveRDS(object=full.combined, file="Seurat.before_Jack.RDS")

### Too long (about 1 day)
full.combined <- JackStraw(object=full.combined, num.replicate=100, dims=100)
full.combined <- ScoreJackStraw(object=full.combined, dims=1:100)

saveRDS(object=full.combined, file="Seurat.after_Jack.RDS")
full.combined <- readRDS("Seurat.after_Jack.RDS")

# PC
JackStrawPlot(object=full.combined, dims=1:41) + NoLegend() # PC41 
JackStrawPlot(object=full.combined, dims=1:100)

DefaultAssay(full.combined) <- "integrated"
full.combined <- RunUMAP(full.combined, reduction = "pca", dims = 1:41)
full.combined <- FindNeighbors(full.combined, reduction = "pca", dims = 1:41)
full.combined <- FindClusters(full.combined, resolution = 0.9)
DimPlot(full.combined, reduction = "umap", label = TRUE, repel=F)

DefaultAssay(full.combined) <- "RNA"

Marker <- c("RBFOX3", "MAP2", "GFAP", "SLC1A2", "CX3CR1", "MRC1", "CD44", "MBP", "CLDN5")
FeaturePlot(object=full.combined, reduction="umap", features=Marker, label=F, min.cutoff = "q9", ncol=3)

levels(full.combined)
AllMarkers <- c("SLC1A2", "GLUL", "AQP4", "GFAP",
                "PDGFRA", "PCDH15", "OLIG2", "OLIG1",
                "MBP", "MOBP", "MAG", "MOG",
                "LAPTM5", "CD74", "SPI1", "CX3CR1",
                "MRC1", "CD44", "CD96", "CYBB",
                "CLDN5", "FLT1",
                "RBFOX3", "MAP2", "STMN2", "SNAP25",
                "SATB2", "SLC17A7", "SLC17A6",
                "GAD1", "GAD2", "SLC32A1",
                "RELN", "FAT2", 
                "RYR1", "HOMER3")
my_cols = c("#d9252f", "#d9252f", "#d9252f", "#d9252f",
            "#2e75d1", "#2e75d1", "#2e75d1", "#2e75d1",
            "#d4c719", "#d4c719", "#d4c719", "#d4c719",
            "#b5287f", "#b5287f", "#b5287f", "#b5287f",
            "#967341", "#967341", "#967341", "#967341",
            "#3aa312", "#3aa312",
            "#232422", "#232422", "#232422", "#232422",
            "#e897c7", "#e897c7", "#e897c7",
            "#FFB500", "#FFB500", "#FFB500",
            "#6c43bf", "#6c43bf",
            "#04cfcb", "#04cfcb")
DotPlot(full.combined,features = AllMarkers) +
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),) +
  labs(x=NULL, y=NULL)


### Cell proportion
count_table <- table(full.combined@meta.data$integrated_snn_res.0.9, full.combined@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell_count_table_PC41_res0.9.csv")

# save
saveRDS(object=full.combined, file="Seurat.Clustering_PC41_res0.9.RDS")
full.combined <- readRDS(file="Seurat.Clustering_PC41_res0.9.RDS")

### Clustering distribution, Gender-specific gene
DefaultAssay(full.combined) <- "RNA"
VlnPlot(object=full.combined, features=c("UTY", "USP9Y", "XIST", "TSIX"), group.by="sample", pt.size=0.1, ncol = 2)


### Calculate a cell cycle score for each cell
cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes
cell.cycle.tirosh[15,2] <- "CENPU"
cell.cycle.tirosh[15,3] <- "ENSG00000151725"
cell.cycle.tirosh[59,2]<-"PIMREG"
cell.cycle.tirosh[59,3]<-"ENSG00000129195"
cell.cycle.tirosh[73,2]<-"JPT1"
cell.cycle.tirosh[73,3]<-"ENSG00000189159"
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes
full.combined <- CellCycleScoring(object=full.combined, s.features=s.genes, g2m.features=g2m.genes)
View(full.combined@meta.data)

DimPlot(full.combined,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(full.combined, reduction="umap", features=metrics, min.cutoff = "q9", ncol=3)

# Save one more time
saveRDS(object=full.combined, file="Seurat.Clustering_PC41_res0.9.RDS")
full.combined <- readRDS(file="Seurat.Clustering_PC41_res1.0.RDS")

### tSNE (long)
DefaultAssay(full.combined) <- "integrated"
full.combined <- RunTSNE(full.combined, dims=1:41)
DimPlot(full.combined, reduction="tsne", label=T, repel=F)

### Save object
saveRDS(object=full.combined, file="Seurat.UMAP_tSNE.RDS")
full.combined <- readRDS(file="Seurat.UMAP_tSNE.RDS")

### Remove cluster (after tSNE)
full.combined <- subset(full.combined, idents=c(20, 21,29,36), invert=T)
DimPlot(full.combined, reduction="umap", label=T, repel=F)
DimPlot(full.combined, reduction="tsne", label=T, repel=F)

saveRDS(object=full.combined, file="Seurat.UMAP_rm.RDS")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(full.combined, reduction="umap", features=metrics, min.cutoff = "q9", ncol=3)


### UMAP (add region, state .....)
DimPlot(object=full.combined, reduction="umap", group.by="group")
DimPlot(object=full.combined, reduction="umap", group.by="sample")
DimPlot(object=full.combined, reduction="umap", group.by="region", raster=T)
DimPlot(object=full.combined, reduction="umap", group.by="state", raster=T)
DimPlot(object=full.combined, reduction="umap", group.by="case", raster=T)
DimPlot(object=full.combined, reduction="umap", group.by="donor")

DimPlot(object=full.combined, reduction="umap", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="umap", split.by="sample", ncol=6)
DimPlot(object=full.combined, reduction="umap", split.by="region", ncol=2)
DimPlot(object=full.combined, reduction="umap", split.by="state", ncol=2)
DimPlot(object=full.combined, reduction="umap", split.by="donor", ncol=4)

DimPlot(object=full.combined, reduction="umap", split.by="integrated_snn_res.0.9", ncol=6)

DimPlot(object=full.combined, reduction="tsne", group.by="group")
DimPlot(object=full.combined, reduction="tsne", group.by="sample")
DimPlot(object=full.combined, reduction="tsne", group.by="region")
DimPlot(object=full.combined, reduction="tsne", group.by="state")
DimPlot(object=full.combined, reduction="tsne", group.by="donor")

DimPlot(object=full.combined, reduction="tsne", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="tsne", split.by="sample", ncol=6)
DimPlot(object=full.combined, reduction="tsne", split.by="region", ncol=2)
DimPlot(object=full.combined, reduction="tsne", split.by="state", ncol=2)
DimPlot(object=full.combined, reduction="tsne", split.by="donor", ncol=4)

DimPlot(object=full.combined, reduction="tsne", split.by="integrated_snn_res.0.9", ncol=6)

### DEG & Heatmap (too long) ### Default RNA!!!!
DefaultAssay(full.combined) <- "RNA"
seurat.markers <- FindAllMarkers(object=full.combined, only.pos=T, min.pct=0.1, logfc.threshold=0.25)

seurat.top10.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
seurat.top20.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

DefaultAssay(full.combined) <- "integrated"
DoHeatmap(object=full.combined, features=seurat.top10.genes$gene, size=4) + NoLegend()

### Save
write.csv(seurat.markers, "Seurat.FindAllMarkers_PCT0.1.csv")
write.csv(seurat.top10.genes, "Seurat.top10markers_PCT0.1.csv")
write.csv(seurat.top20.genes, "Seurat.top20markers_PCT0.1.csv")
saveRDS(object=seurat.markers, file="Seurat.FindAllMarker_PCT0.1.RDS")


### Visualizing marker expression (Neuron, Astrocyte, Microglia, Oligo, Endo)
# Neuron: SMAP25, STMN2, RBFOX3, MAP2
# Ex: SATB2, SLC17A7, SLC17A6
# In: GAD1, GAD2, SLC32A1
# Granule: RBFOX3, RELN, GABRA6, FAT2
# Purkinje: ITPR1, RORA, SORCS3, RYR1, GAD1, GAD2
# Molecular layer interneuron (basket, stellate): PTPRK, PRKCD, LYPD6
# Purkinje layer interneuron: KLHL1
# Interneuron: PAX2, GAD1/2
# GPC: ATOH1 X, DCC, ROBO1, PAX6, PTCH1, NEUROD1, ZIC1, ZIC2, SEMA6A, PLXNA2
# iCN: MEIS2
# Astrocyte: SLC1A2, GLUL, SOX9, AQP4, GJA1, NDRG2, GFAP, ALDH1A1, ALDH1L1, VIM
# Microglia: LAPTM5, CD74, SPI1, TMEM119, CX3CR1, HEXB
# BAM: MRC1
# Macro: CD96, CD44, MRC1, CYBB, MS4A7, AXL
# OPC: PTGDS, PDGFRA, PCDH15, OLIG2, OLIG1
# Oligo: PLP1, MAG, MOG, MOBP, MBP
# Endo: CLDN5, VTN, FLT1, RAMP2, KDR
### Cell type marker (Neuron, Ex, In, Astro, Micro, Endo, Oligo, OPC)
Neuron <- c("MAP2", "STMN2",
            "SATB2", "SLC17A7",
            "GAD1",
            "RBFOX3", "RELN",
            "RYR1", "SORCS3")
NonNeuron <- c("SLC1A2", "GFAP",
               "PDGFRA", "PCDH15",
               "MBP", "MOG",
               "CX3CR1", "CD44",
               "CLDN5")
AllMarkers <- c("MAP2", "STMN2", "SNAP25",
                "SATB2", "SLC17A7", "SLC17A6",
                "GAD1", "GAD2", "SLC32A1",
                "RBFOX3", "RELN", "FAT2", 
                "ITPR1", "RYR1", "SORCS3", "NXPH1",
                "SLC1A2", "GLUL", "GFAP",
                "PDGFRA", "PCDH15",
                "MBP", "MOG", "PLP1",
                "CD74", "CX3CR1",
                "MRC1", "CD44", "CD96",
                "CLDN5", "FLT1")
my_cols = c("#232422", "#232422", "#232422",
            "#6c43bf", "#6c43bf", "#6c43bf",
            "#e897c7", "#e897c7", "#e897c7",
            "#c8bdf0", "#c8bdf0", "#c8bdf0",                 ##13417d
            "#f7c7a1", "#f7c7a1", "#f7c7a1", "#f7c7a1",      ##04cfcb
            "#d9252f", "#d9252f", "#d9252f",
            "#FFB500", "#FFB500",
            "#967341", "#967341", "#967341",
            "#3aa312", "#3aa312",
            "#3aa312", "#3aa312", "#3aa312",
            "#2e75d1", "#2e75d1")

DefaultAssay(full.combined) <- "RNA"
FeaturePlot(object=full.combined, reduction="umap", features=Neuron, label=F, min.cutoff = "q9", ncol=3)

pdf("DotPlot_All.pdf", width = 14, height = 9)
DotPlot(full.combined,features = AllMarkers, cols=c("#e8e8e8", "#751665")) +  ##0b488a #a10d63 #751665 #0a4d54
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nexpression"),
         size = guide_legend("Percent\nexpressed"),)+
  labs(x=NULL, y=NULL)
dev.off()

#group 5/6.5 sample 8/6.5
tiff(filename = "Cell proportion/Cell_proportion_region_labeling.tif", width = 10, height = 3, units = "in", res=300)
dittoBarPlot(
  object = full.combined,
  var = "region",
  group.by = "Newcluster",
  main="",
  xlab="",
  var.labels.reorder=c(1,2,13,12,11,3,4,5,6,7,8,9,10,14,15))+
  theme(axis.text=element_text(size=8, color="black"),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.5, 'cm'))
dev.off()

# Cluster Dendrogram
DefaultAssay(full.combined) <- "integrated"
full.combined <- BuildClusterTree(full.combined,assay="integrated")
myPhyTree <- Tool(object=full.combined, slot = "BuildClusterTree")
#tiff(filename = "Marker/Dendrogram4.tif", width = 5, height = 11, units = "in", res=300)
pdf("Marker/Dendrogram.pdf", width = 5, height = 11)
ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1,
                tip.col=c("#3aa312", "#967341", "#967341", "#967341", "#3aa312", "#c8bdf0", "#c8bdf0", "#c8bdf0",  #1~8
                          "#c8bdf0", "#c8bdf0", "#c8bdf0", "#f7c7a1", "#a3a3a3", "#a3a3a3", "#a3a3a3", "#2e75d1",  #9~16
                          "#3aa312", "#d9252f", "#d9252f", "#e897c7", "#e897c7", "#e897c7", "#e897c7", "#e897c7",  #17~24
                          "#e897c7", "#FFB500", "#FFB500", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf",  #25~32
                          "#6c43bf", "#6c43bf")) +                                                                 #33~34
  geom_tiplab()
dev.off()

levels(full.combined)

levels(full.combined) <- c(0,1,2,3,4,5,6,7,8,9,10,
                           11,12,13,14,15,16,17,18,19,
                           22,23,24,25,26,27,28,30,
                           31,32,33,34,35,37)


### Convert Cluster ID
levels(full.combined)
new.cluster.ids <- c("Micro_Phago", "Oligo1", "Oligo2", "Oligo3", "Micro", "Gran2", "Gran4", "Gran1", "Gran3",
                     "Gran6", "Gran5", "Purk", "PQ3", "PQ1", "PQ2", "Endo", "Macro", "Astro1", "Astro2",
                     "InN4_MGE", "InN5_MGE_PVALB", "InN6_MGE_NOS1", "InN1_MGE_SST", "InN2_CGE_VIP", "InN3_CGE_SST",
                     "OPC1", "OPC2",
                     "ExN5_L4-6", "ExN7_L6", "ExN2_L2-4", "ExN3_L4-5", "ExN4_L4-5", "ExN1_L2-4", "ExN6_L2-5")
names(new.cluster.ids) <- levels(full.combined)
full.combined <- RenameIdents(full.combined, new.cluster.ids)
full.combined@meta.data$"Newcluster" <- as.factor(full.combined@active.ident)
pdf("UMAP_label.pdf", width = 10.367, height = 7)
DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=4, repel=T, raster=T)
dev.off()

levels(full.combined)
levels(full.combined) <- c("Micro_Phago", "Oligo1", "Oligo2", "Oligo3", "Micro", "Gran2", "Gran4", "Gran1", "Gran3",
                           "Gran6", "Gran5", "Purk", "PQ3", "PQ1", "PQ2", "Endo", "Macro", "Astro1", "Astro2",
                           "InN4", "InN5", "InN6", "InN1", "InN2", "InN3", "OPC1", "OPC2", "ExN5", "ExN7",
                           "ExN2", "ExN3", "ExN4", "ExN1", "ExN6")
levels(full.combined) <- c("Astro1", "Astro2", "OPC1", "OPC2", "Oligo1", "Oligo2", "Oligo3",
                           "Micro", "Micro_Phago", "Macro", "Endo", "PQ1", "PQ2", "PQ3",
                           "ExN1_L2-6", "ExN2_L2-4", "ExN3_L4-5", "ExN4_L4-5", "ExN5_L4-6", "ExN6_L2-5", "ExN7_L6",
                           "InN1_MGE_SST", "InN2_CGE_VIP", "InN3_CGE_SST", "InN4_MGE", "InN5_MGE_PVALB", "InN6_MGE_NOS1", 
                           "Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Purk")


# DotPlot
AllMarkers <- c("MAP2", "STMN2", "SNAP25",
                "SATB2", "SLC17A7", "SLC17A6",
                "GAD1", "GAD2", "SLC32A1",
                "RBFOX3", "RELN", "FAT2", 
                "HOMER3", "RYR1", "SORCS3",
                "SLC1A2", "GLUL", "GFAP",
                "PDGFRA", "PCDH15",
                "MBP", "MOG", "PLP1",
                "CD74", "CX3CR1",
                "MRC1", "CD44", "CD96",
                "CLDN5", "FLT1")
my_cols = c("#232422", "#232422", "#232422",
            "#6c43bf", "#6c43bf", "#6c43bf",
            "#e897c7", "#e897c7", "#e897c7",
            "#13417d", "#13417d", "#13417d",
            "#04cfcb", "#04cfcb", "#04cfcb",
            "#d9252f", "#d9252f", "#d9252f",
            "#FFB500", "#FFB500",
            "#967341", "#967341", "#967341",
            "#3aa312", "#3aa312",
            "#3aa312", "#3aa312", "#3aa312",
            "#2e75d1", "#2e75d1")

DefaultAssay(full.combined) <- "RNA"

tiff(filename = "DotPlot_All_Large.tif", width = 13, height = 7, units = "in", res=300)
DotPlot(full.combined,features = AllMarkers, cols=c("#e8e8e8", "#751665")) +  ##0b488a #a10d63 #751665 #0a4d54
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),)+
  labs(x=NULL, y=NULL)
dev.off()

tiff(filename = "VlnPlot_All1.tif", width = 13, height = 9, units = "in", res=300)
VlnPlot(full.combined, features=AllMarkers, stack=T, sort=F, flip=F) +
  theme(legend.position = "none")+
  theme(strip.text.x = element_text(face="italic", size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=270))+
  labs(x=NULL, y=NULL)
dev.off()

saveRDS(object=full.combined, file="Seurat.Final.RDS")
full.combined <- readRDS(file="Seurat.Final.RDS")

pdf("UMAP_region.pdf", width = 7.695, height = 7)
DimPlot(object=full.combined, reduction="umap", group.by="region", raster=T)
dev.off()

pdf("UMAP_case.pdf", width = 7.611, height = 7)
DimPlot(object=full.combined, reduction="umap", group.by="state", raster=T)
dev.off()

### Large cluster
levels(full.combined)
new.cluster.ids <- c("Micro", "Oligo", "Oligo", "Oligo", "Micro", "Gran", "Gran", "Gran", "Gran",
                     "Gran", "Gran", "Purk", "PQ", "PQ", "PQ", "Endo", "Micro", "Astro", "Astro",
                     "InN", "InN", "InN", "InN", "InN", "InN",
                     "OPC", "OPC",
                     "ExN", "ExN", "ExN", "ExN", "ExN", "ExN", "ExN")
names(new.cluster.ids) <- levels(full.combined)
full.combined <- RenameIdents(full.combined, new.cluster.ids)
full.combined@meta.data$"Largecluster" <- as.factor(full.combined@active.ident)
DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=4, repel=T)
levels(full.combined) <- c("ExN", "InN", "Gran", "Purk", "PQ", "Astro", "OPC", "Oligo", "Micro", "Endo")


saveRDS(object=full.combined, file="Seurat.Large.RDS")
full.combined_Large <- readRDS("Seurat.Large.RDS")
levels(full.combined_Large)
levels(full.combined_Large) <- c("ExN", "InN", "Gran", "Purk", "PQ", "Astro", "OPC", "Oligo", "Micro", "Endo")

### color
my_cols <- c('ExN'='#918fe3', 'InN'='#ed91bf','Gran'='#c8bdf0', 'Purk'='#f7c7a1', 'PQ'='#d4d2d2',
             'Astro'='#ed8b28', 'OPC'='#FFB500', 'Oligo'='#967341', 'Micro'='#A4E804', 
             'Endo'='#2e75d1')

#'Astro1'='#f2d50f', 'Astro2'='#edab19', 'Astro3'='#e4622a', 'Astro4'='#dd223a',
#'OPC'='#f0b659', 'Oligo1'='#A05837', 'Oligo2'='#ba8a1a', 'Oligo3'='#d1860d',
#'Micro_Homeo'='#0d914d', 'Micro_Phago'='#25ba46', 'BAM'='#b1f70a', 'Endo'='#31332e',
#'Doublet'='#44374a', 'ExN1_L5'='#a4eaed', 'ExN2_L2-4'='#0dbfef', 'ExN3_L2-4'='#0895e5',
#'ExN4_L4-6'='#0470dc', 'ExN5_L4-6'='#2d9cc1', 'ExN6_L4-6'='#3a79ac', 'ExN7_L5-6'='#475798',
#'ExN8_L2-6'='#2e49b3', 'ExN9_L4-6'='#1f46de', 'ExN10_L5-6'='#3e4c80', 'ExN11_L6'='#35338f',
#'InN1_SST'='#ecaded', 'InN2_VIP'='#d865db', 'InN3_PVALB'='#8c76f3', 'InN4_LAMP5'='#9c58cb',
#'InN5'='#915c91', 'InN6_SST'='#6d31d4', 'InN7_PVALB'='#b631d4', 'InN8_NOS1'='#bf049a',
#'Gran1'='#c4c0b7', 'Gran2'='#d9cca9', 'Gran3'='#9c9070', 'Gran4'='#6e6c5f', 'Purkinje'='#696440'
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)

#tiff(filename = "UMAP_Large.tif", width = 8, height = 7, units = "in", res=300)
pdf("UMAP_Large.pdf", width = 7.945, height = 7)
DimPlot(full.combined_Large, reduction = "umap", cols=my_cols2, label = TRUE, label.size=3, repel=T, raster=T)
dev.off()


DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = T, label.size=4, repel=T)
DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = T, label.size=4, repel=T)
DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = F, label.size=4, repel=T)

DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=4, repel=T)

DimPlot(object=full.combined, reduction="umap", split.by="group", ncol=2, cols=my_cols2)
DimPlot(object=full.combined, reduction="umap", split.by="sample", ncol=6, cols=my_cols2)
DimPlot(object=full.combined, reduction="umap", split.by="state", ncol=2, cols=my_cols2, label=T, repel=T)+ NoLegend()

saveRDS(object=full.combined, file="Seurat.Clustering_newlabel.RDS")
full.combined1 <- readRDS(file="Seurat.Clustering_newlabel.RDS")
