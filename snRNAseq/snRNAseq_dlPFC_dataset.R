library(ggtree)
library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(Seurat)
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
k3 <- Read10X(data.dir = "kang3/outs/filtered_feature_bc_matrix/")
k17 <- Read10X(data.dir = "kang17/outs/k17_filtered_feature_bc_matrix/")
k19 <- Read10X(data.dir = "kang19/outs/filtered_feature_bc_matrix/")
k23 <- Read10X(data.dir = "kang23/outs/filtered_feature_bc_matrix/")
k29 <- Read10X(data.dir = "kang29/outs/filtered_feature_bc_matrix/")
k33 <- Read10X(data.dir = "kang33/outs/filtered_feature_bc_matrix/")
k35 <- Read10X(data.dir = "kang35/outs/filtered_feature_bc_matrix/")
k37 <- Read10X(data.dir = "kang37/outs/filtered_feature_bc_matrix/")
k39 <- Read10X(data.dir = "kang39/outs/filtered_feature_bc_matrix/")
k41 <- Read10X(data.dir = "kang41/outs/filtered_feature_bc_matrix/")
k43 <- Read10X(data.dir = "kang43/outs/filtered_feature_bc_matrix/")

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

### remove matrix
rm (k1, k3, k17, k19, k23, k29, k33, k35, k37, k39, k41, k43)

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

### merge datasets for plot
alldata_before <- merge(dlPFC_C1, c(dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                             dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6))
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
pdf("Before_QC_nodot.pdf", width = 11, height = 6)
VlnPlot(alldata_before, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()
dev.off()
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
alldata_before@meta.data <- metadata1

metadata1$group <- NA
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
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

VlnPlot(alldata_before, group.by = "group", features = feats, pt.size = 0.1, ncol = 5) +
  NoLegend()


# save
saveRDS(object=alldata_before, file="Alldata_before_QC.RDS")
alldata_before <- readRDS("Alldata_before_QC.RDS")

### Check top 5% nCount by quantile function
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

alldata_after1 <- merge(dlPFC_C1, c(dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                             dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6))

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
alldata_after1@meta.data <- metadata2

metadata2$group <- NA
metadata2$group[which(str_detect(metadata2$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata2$group[which(str_detect(metadata2$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
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

### combine list
combined.list <- list(dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                      dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6)
names(combined.list) <- c("dlPFC_C1", "dlPFC_C2", "dlPFC_C3", "dlPFC_C4", "dlPFC_C5", "dlPFC_C6",
                          "dlPFC_M1", "dlPFC_M2", "dlPFC_M3", "dlPFC_M4", "dlPFC_M5", "dlPFC_M6")

### Normalize and identify variable features for each dataset independently
combined.list <- lapply(X = combined.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

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

sweep.res <- paramSweep_v3(combined.list$dlPFC_M6, PCs=1:10)
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


homotypic.prop <- modelHomotypic(combined.list$dlPFC_M6@meta.data$RNA_snn_res.0.8)
nExp <- round(ncol(combined.list$dlPFC_M6) * 0.08)  # expect 8% doublets
nExp.adj <- round(nExp*(1 - homotypic.prop))
combined.list$dlPFC_M6 <- doubletFinder_v3(combined.list$dlPFC_M6, pN = 0.25, pK = 0.08, nExp = nExp,
                                           PCs = 1:10, reuse.pANN=F)
head(combined.list$dlPFC_M6@meta.data)[9]
combined.list$dlPFC_M6 <- doubletFinder_v3(combined.list$dlPFC_M6, pN = 0.25, pK = 0.08, nExp = nExp.adj,
                                           PCs = 1:10, reuse.pANN="pANN_0.25_0.08_333")

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

# save
saveRDS(object=combined.list, file="Seurat.remove.DF_Homotypic.RDS")

quantile(combined.list$dlPFC_M6$nFeature_RNA,probs = c(0.5))
quantile(combined.list$dlPFC_M6$nCount_RNA,probs = c(0.5))

alldata_after <- merge(combined.list$dlPFC_C1, c(combined.list$dlPFC_C2, combined.list$dlPFC_C3, combined.list$dlPFC_C4, combined.list$dlPFC_C5, combined.list$dlPFC_C6,
                                                 combined.list$dlPFC_M1, combined.list$dlPFC_M2, combined.list$dlPFC_M3, combined.list$dlPFC_M4, combined.list$dlPFC_M5, combined.list$dlPFC_M6))

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
pdf("After_QCDF_nodot.pdf", width = 11, height = 6)
VlnPlot(alldata_after, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()
dev.off()

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
alldata_after@meta.data <- metadata3

metadata3$group <- NA
metadata3$group[which(str_detect(metadata3$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata3$group[which(str_detect(metadata3$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
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

VlnPlot(alldata_after, group.by = "group", features = feats, pt.size = 0.1, ncol = 5) +
  NoLegend()


# save (alldata before, after)
saveRDS(object=alldata_after, file="Alldata_after_QC&DF.RDS")
alldata_after <- readRDS("Alldata_after_QC&DF.RDS")

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
full.combined@meta.data <- metadata

metadata$group <- NA
metadata$group[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "dlPFC_C"
metadata$group[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "dlPFC_M"
full.combined@meta.data <- metadata

metadata$case <- NA
metadata$case[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "Con"
metadata$case[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "MDD"
full.combined@meta.data <- metadata

VlnPlot(full.combined, group.by = "group", features = feats, pt.size = 0.1, ncol = 5) +
  NoLegend()


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

rm(alldata_after, alldata_after1, alldata_before, dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
   dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6)

### Too long
full.combined <- JackStraw(object=full.combined, num.replicate=100, dims=100)
full.combined <- ScoreJackStraw(object=full.combined, dims=1:100)###

saveRDS(object=full.combined, file="Seurat.after_Jack.RDS")
full.combined <- readRDS(file="Seurat.after_Jack.RDS")

# PC
JackStrawPlot(object=full.combined, dims=1:44) + NoLegend() # PC44 or PC50
JackStrawPlot(object=full.combined, dims=1:100)

DefaultAssay(full.combined) <- "integrated"
full.combined <- RunUMAP(full.combined, reduction = "pca", dims = 1:44)
full.combined <- FindNeighbors(full.combined, reduction = "pca", dims = 1:44)
full.combined <- FindClusters(full.combined, resolution = 0.8)
DimPlot(full.combined, reduction = "umap", label = TRUE, repel=F, raster=F)

count_table <- table(full.combined@meta.data$Newcluster, full.combined@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell_count_table_labeling.csv")


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
        split.by="Phase")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(full.combined, reduction="umap", features=metrics, min.cutoff = "q9", ncol=3)



# save
saveRDS(object=full.combined, file="Seurat.Clustering_PC44_res0.8.RDS")

### tSNE (long)
DefaultAssay(full.combined) <- "integrated"
full.combined <- RunTSNE(full.combined, dims=1:44)
DimPlot(full.combined, reduction="tsne", label=T, repel=F, raster=F)

### Save object
saveRDS(object=full.combined, file="Seurat.UMAP_tSNE.RDS")
full.combined <- readRDS("Seurat.UMAP_tSNE.RDS")

### Remove cluster (after tSNE)
full.combined <- subset(full.combined, idents=c(0, 23, 18, 24, 26), invert=T)
DimPlot(full.combined, reduction="umap", label=T, repel=F, raster=F)
DimPlot(full.combined, reduction="tsne", label=T, repel=F, raster=F)

saveRDS(object=full.combined, file="Seurat.UMAP_rm.RDS")
full.combined <- readRDS("Seurat.UMAP_rm.RDS")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(full.combined, reduction="umap", features=metrics, min.cutoff = "q9", ncol=3)

### Clustering distribution, Gender-specific gene
DefaultAssay(full.combined) <- "RNA"
VlnPlot(object=full.combined, features=c("UTY", "USP9Y", "XIST", "TSIX"), group.by="sample", ncol=2)



### UMAP
DimPlot(object=full.combined, reduction="umap", group.by="group")
DimPlot(object=full.combined, reduction="umap", group.by="case", raster=T)
DimPlot(object=full.combined, reduction="umap", group.by="sample")

DimPlot(object=full.combined, reduction="umap", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="umap", split.by="sample", ncol=6)

DimPlot(object=full.combined, reduction="umap", split.by="integrated_snn_res.0.8", ncol=6)

DimPlot(object=full.combined, reduction="tsne", group.by="group")
DimPlot(object=full.combined, reduction="tsne", group.by="sample")

DimPlot(object=full.combined, reduction="tsne", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="tsne", split.by="sample", ncol=6)



### Visualizing marker expression (Neuron, Astrocyte, Microglia, Oligo, Endo)
# Neuron: SMAP25, STMN2, RBFOX3, MAP2
# Astrocyte: SLC1A2, GLUL, SOX9, AQP4, GJA1, NDRG2, GFAP, ALDH1A1, ALDH1L1, VIM
# Microglia: CD74, SPI1, MRC1, TMEM119, CX3CR1
# OPC: PTGDS, PDGFRA, PCDH15, OLIG2, OLIG1
# Oligo: PLP1, MAG, MOG, MOBP, MBP
# Ex: SATB2, SLC17A7, SLC17A6
# In: GAD1, GAD2, SLC32A1
# Endo: CLDN5, VTN
### Cell type marker (Neuron, Ex, In, Astro, Micro, Endo, Oligo, OPC)
Neuron <- c("RBFOX3", "MAP2", "STMN2",
            "SATB2", "SLC17A7", "SLC17A6",
            "GAD1", "GAD2", "SLC32A1")
NonNeuron <- c("SLC1A2", "GFAP",
               "PDGFRA", "PCDH15",
               "MBP", "MOG",
               "CX3CR1", "CD44",
               "CLDN5")
Astro <- c("SLC1A2", "GLUL", "SOX9", "AQP4", "GJA1", "NDRG2", "ALDH1A1", "ALDH1L1", "VIM")
Micro <- c("LAPTM5", "CD74", "SPI1", "TMEM119", "CX3CR1", "HEXB")
Macro <- c("CD96", "CD44", "MRC1", "CYBB", "MS4A7", "AXL")
Endo <- c("CLDN5", "VTN", "FLT1", "RAMP2", "KDR")
Oligo <- c("MBP", "MOBP", "MAG", "PLP1", "MOG",
                   "PDGFRA", "PCDH15", "OLIG2", "OLIG1")
AllMarkers <- c("RBFOX3", "MAP2", "STMN2", "SNAP25",
                "SATB2", "SLC17A7", "SLC17A6",
                "GAD1", "GAD2", "SLC32A1",
                "SLC1A2", "GFAP",
                "CD74", "CX3CR1",
                "MRC1", "CD44", "CD96", "CYBB",
                "MBP", "MOG",
                "PDGFRA", "PCDH15", 
                "CLDN5", "FLT1")

DefaultAssay(full.combined) <- "RNA"

FeaturePlot(object=full.combined, reduction="umap", features=NonNeuron, label=F, min.cutoff = "q9", ncol=3)
VlnPlot(object=full.combined, features=Oligo, pt.size=0, stack = T, flip = T) + NoLegend()
DotPlot(full.combined,features = AllMarkers) +
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),)+
  labs(x=NULL, y=NULL)


### Cell proportion
count_table <- table(full.combined@meta.data$integrated_snn_res.0.9, full.combined@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell_count_table.csv")
dittoBarPlot(
  object = full.combined,
  var = "Newcluster",
  group.by = "sample")

dittoBarPlot(
  object = full.combined,
  var = "sample",
  group.by = "Newcluster")


### DEG & Heatmap (too long)
DefaultAssay(full.combined) <- "RNA"
seurat.markers <- FindAllMarkers(object=full.combined, only.pos=T, min.pct=0.1, logfc.threshold=0.25) # 0.1 / 0.25
head(seurat.markers)

seurat.top10.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
seurat.top20.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

DefaultAssay(full.combined) <- "integrated"
DoHeatmap(object=full.combined, features=seurat.top10.genes$gene, size=4) + NoLegend()

### Save
write.csv(seurat.markers, "Seurat.FindAllMarkers_PCT0.1.csv")
write.csv(seurat.top10.genes, "Seurat.top10markers_PCT0.1.csv")
write.csv(seurat.top20.genes, "Seurat.top20markers_PCT0.1.csv")
saveRDS(object=seurat.markers, file="Seurat.FindAllMarker_PCT0.1.RDS")


### Cluster Dendrogram
DefaultAssay(full.combined) <- "integrated"
full.combined <- BuildClusterTree(full.combined, assay="integrated")
myPhyTree <- Tool(object=full.combined, slot = "BuildClusterTree")
pdf("Dendrogram.pdf", width = 4.5, height = 11)
ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1,
                tip.col=c("#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf",    #1~8
                          "#e897c7", "#e897c7", "#e897c7", "#e897c7", "#e897c7", "#e897c7", "#e897c7", "#d9252f",    #9~16
                          "#FFB500", "#FFB500", "#d9252f", "#d9252f", "#3aa312", "#3aa312", "#967341", "#967341",    #17~28
                          "#967341", "#967341", "#a3a3a3", "#a3a3a3", "#a3a3a3", "#3aa312", "#2e75d1", "#3aa312")) + #29~36
  geom_tiplab()
dev.off()

levels(full.combined)


### Convert Cluster ID
levels(full.combined)

new.cluster.ids <- c("ExN1", "Oligo1", "ExN2", "Oligo2", "ExN3", "Oligo3", "Micro", "Astro1", "OPC1", "Astro2",
                     "Oligo4", "InN1", "InN2", "ExN4", "InN3", "InN4", "Astro3", "ExN5", "ExN6",
                     "InN5", "InN6", "Endo", "ExN7", "InN7", "Micro_Phago", "ExN8",
                     "InN8", "ExN9", "Macro", "BAM", "ExN10", "OPC2")
names(new.cluster.ids) <- levels(full.combined)
full.combined <- RenameIdents(full.combined, new.cluster.ids)
full.combined@meta.data$"Newcluster" <- as.factor(full.combined@active.ident)
DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=4, repel=F) + NoLegend()


levels(full.combined) <- c("ExN7", "ExN9", "ExN10", "ExN5", "ExN6", "ExN8", "ExN3", "ExN4",
                           "InN5", "InN7", "InN1", "InN8", "InN3", "InN2", "InN4",
                           "Astro3", "OPC1", "OPC2", "Astro1", "Astro2", "Micro", "Micro_Phago",
                           "Oligo1", "Oligo4", "Oligo2", "Oligo3", "InN6", "ExN1", "ExN2",
                           "BAM", "Endo", "Macro")

levels(full.combined)

new.cluster.ids <- c("ExN5_L4-6", "ExN7_L5-6", "ExN8_L6", "ExN3_L4-5", "ExN4_L4-5", "ExN6_L4-5", "ExN1_L2-4", "ExN2_L2-4",
                     "InN5_MGE", "InN6_MGE_PVALB", "InN1_MGE_SST", "InN7_MGE_NOS1", "InN3_CGE_VIP", "InN2_CGE", "InN4_CGE_LAMP5",
                     "Astro3", "OPC1", "OPC2", "Astro1", "Astro2", "Micro", "Micro_Phago",
                     "Oligo1", "Oligo4", "Oligo2", "Oligo3", "PQ3", "PQ1", "PQ2",
                     "BAM", "Endo", "Macro")
names(new.cluster.ids) <- levels(full.combined)
full.combined <- RenameIdents(full.combined, new.cluster.ids)
full.combined@meta.data$"Newcluster" <- as.factor(full.combined@active.ident)
DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=3, repel=F) + NoLegend()


DefaultAssay(full.combined) <- "RNA"
AllMarkers <- c("RBFOX3", "MAP2", "STMN2","SNAP25",
                "SATB2", "SLC17A7", "SLC17A6",
                "GAD1", "GAD2", "SLC32A1",
                "SLC1A2", "GLUL", "GFAP",
                "PDGFRA", "PCDH15",
                "MBP", "MOG", "PLP1",
                "CD74", "CX3CR1",
                "MRC1", "CD44", "CD96",
                "CLDN5", "FLT1")
my_cols = c("#232422", "#232422", "#232422", "#232422",
            "#6c43bf", "#6c43bf", "#6c43bf",
            "#e897c7", "#e897c7", "#e897c7",
            "#d9252f", "#d9252f", "#d9252f",
            "#FFB500", "#FFB500",
            "#967341", "#967341", "#967341",
            "#3aa312", "#3aa312",
            "#3aa312", "#3aa312", "#3aa312",
            "#2e75d1", "#2e75d1")

DefaultAssay(full.combined) <- "RNA"
pdf("DotPlot_All.pdf", width = 9, height = 9)
DotPlot(full.combined,features = AllMarkers, cols=c("#e8e8e8", "#751665"), dot.scale = 7) +  ##0b488a #a10d63 #751665 #0a4d54
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=9))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),)+
  labs(x=NULL, y=NULL)
dev.off()

saveRDS(object=full.combined, file="Seurat.Final.RDS")
full.combined <- readRDS ("Seurat.Final.RDS")

levels(full.combined) <- c('Astro1', 'Astro2', 'Astro3',
                           'OPC1', 'OPC2',
                           'Oligo1', 'Oligo2', 'Oligo3', 'Oligo4',
                           'Micro', 'Micro_Phago', 'BAM', 'Macro',
                           'Endo',
                           'PQ1', 'PQ2', 'PQ3',
                           'ExN1_L2-4', 'ExN2_L2-4', 'ExN3_L4-5', 'ExN4_L4-5',
                           'ExN5_L4-6', 'ExN6_L4-5', 'ExN7_L5-6', 'ExN8_L6',
                           'InN1_MGE_SST', 'InN2_CGE', 'InN3_CGE_VIP', 'InN4_CGE_LAMP5',
                           'InN5_MGE', 'InN6_MGE_PVALB', 'InN7_MGE_NOS1')

levels(full.combined)

### UMAP color
levels(full.combined)

my_cols <- c('Astro1'='#ed8b28', 'Astro2'='#e4622a', 'Astro3'='#e33b50',
             'OPC1'='#FFB500', 'OPC2'='#cf8104',
             'Oligo1'='#c2a442', 'Oligo2'='#ba8a1a', 'Oligo3'='#967341', 'Oligo4'='#A05837',
             'Micro'='#A4E804', 'Micro_Phago'='#3aa312', 'BAM'='#5b8a5a', 'Macro'='#02705f',
             'Endo'='#2e75d1',
             'PQ1'='#ccc2c2', 'PQ2'='#a3a3a3', 'PQ3'='#8c8484',
             'ExN1_L2-4'='#c8bdf0', 'ExN2_L2-4'='#918fe3', 'ExN3_L4-5'='#7e7cd6', 'ExN4_L4-5'='#b06ae6',
             'ExN5_L4-6'='#8e6cd4', 'ExN6_L4-5'='#5d58e0', 'ExN7_L5-6'='#6f7fa6', 'ExN8_L6'='#5463a8',
             'InN1_MGE_SST'='#e07b91', 'InN2_CGE'='#ed91bf', 'InN3_CGE_VIP'='#FFAA92', 'InN4_CGE_LAMP5'='#f7c7a1',
             'InN5_MGE'='#C8A1A1', 'InN6_MGE_PVALB'='#c9908b', 'InN7_MGE_NOS1'='#7B4F4B')
my_cols <- c('ExN'='#918fe3', 'InN'='#ed91bf', 'PQ'='#d4d2d2',
             'Astro'='#ed8b28', 'OPC'='#FFB500', 'Oligo'='#967341', 'Micro'='#A4E804', 
             'Endo'='#2e75d1')
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]

pdf("UMAP_Final2.pdf", width = 7.125, height = 7)
DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = TRUE, label.size=3, repel=T, raster=T) + NoLegend()
dev.off()

pdf("tSNE_Final.pdf", width = 10.463, height = 7)
DimPlot(full.combined, reduction="tsne", label=T, cols=my_cols2, label.size=3, repel=F, raster=T)
dev.off()

scales::show_col(my_cols2)

### Cell proportion
full.combined <- subset(full.combined, idents=c("PQ1", "PQ2", "PQ3"), invert=T)
pdf("Cell_proportion_sample.pdf", width = 5.5, height = 6)
dittoBarPlot(
  object = full.combined,
  var = "Newcluster",
  group.by = "sample",
  main="",
  xlab="",
  color.panel = my_cols2,
  var.labels.reorder=c(1,2,3,28,29,24,25,26,27,22,23,4,21,5,30,31,32,
                       6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))+
  theme(axis.text=element_text(size=7, color="black"),
        legend.text=element_text(size=9),
        legend.key.size = unit(0.5, 'cm'))
dev.off()
