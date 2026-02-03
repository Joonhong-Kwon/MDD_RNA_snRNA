remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(Seurat)
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
library(fields)
library(ggdendro)
library(metap)
library(ggtree)
suppressMessages(require(DoubletFinder))

### Load data
k2 <- Read10X(data.dir = "kang2/outs/filtered_feature_bc_matrix/")
k4 <- Read10X(data.dir = "kang4/outs/filtered_feature_bc_matrix/")
k18 <- Read10X(data.dir = "kang18/outs/k18_filtered_feature_bc_matrix/")
k20 <- Read10X(data.dir = "kang20/outs/k20_filtered_feature_bc_matrix/")
k24 <- Read10X(data.dir = "kang24/outs/filtered_feature_bc_matrix/")
k30 <- Read10X(data.dir = "kang30/outs/filtered_feature_bc_matrix/")
k34 <- Read10X(data.dir = "kang34/outs/filtered_feature_bc_matrix/")
k36 <- Read10X(data.dir = "kang36/outs/filtered_feature_bc_matrix/")
k38 <- Read10X(data.dir = "kang38/outs/filtered_feature_bc_matrix/")
k40 <- Read10X(data.dir = "kang40/outs/filtered_feature_bc_matrix/")
k42 <- Read10X(data.dir = "kang42/outs/filtered_feature_bc_matrix/")
k44 <- Read10X(data.dir = "kang44/outs/filtered_feature_bc_matrix/")

### Create Seurat object
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
rm (k2, k4, k18, k20, k24, k30, k34, k36, k38, k40, k42, k44)

### QC metrics and each sample cell filtering
### Mitochondrial genes
mito_genes <- rownames(CBC_C1)[grep("^MT-", rownames(CBC_C1))]
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
ribo_genes <- rownames(CBC_C1)[grep("^RP[SL]", rownames(CBC_C1))]
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
hemo_genes <- rownames(CBC_C1)[grep("^HB[^(P)]", rownames(CBC_C1))]
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
alldata_before <- merge(CBC_C1, c(CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                                  CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6))
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
metadata1$group[which(str_detect(metadata1$orig.ident, "^CBC_C"))] <- "CBC_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^CBC_M"))] <- "CBC_M"
alldata_before@meta.data <- metadata1

### Visualize the number UMIs/transcripts per cell
metadata1 %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 9000)

# Visualize the distribution of genes detected per cell via histogram
metadata1 %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(400))

VlnPlot(alldata_before, group.by = "group", features = feats, pt.size = 0, ncol = 5) +
  NoLegend()


# save
saveRDS(object=alldata_before, file="Alldata_before_QC.RDS")
alldata_before <- readRDS("Alldata_before_QC.RDS")

### Check top 5% nCount by quantile function
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

alldata_after1 <- merge(CBC_C1, c(CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                                  CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6))

VlnPlot(alldata_after1, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

### Histogram
metadata2 <- alldata_after1@meta.data
metadata2$sample <- NA
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
mito_genes <- rownames(CBC_C1)[grep("^MT-", rownames(CBC_C1))]
mito_genes
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
combined.list <- list(CBC_C1, CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
                      CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6)
names(combined.list) <- c("CBC_C1", "CBC_C2", "CBC_C3", "CBC_C4", "CBC_C5", "CBC_C6",
                          "CBC_M1", "CBC_M2", "CBC_M3", "CBC_M4", "CBC_M5", "CBC_M6")

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
head(combined.list$CBC_C1@meta.data)
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

alldata_after <- merge(combined.list$CBC_C1, c(combined.list$CBC_C2, combined.list$CBC_C3, combined.list$CBC_C4, combined.list$CBC_C5, combined.list$CBC_C6,
                                                 combined.list$CBC_M1, combined.list$CBC_M2, combined.list$CBC_M3, combined.list$CBC_M4, combined.list$CBC_M5, combined.list$CBC_M6))

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
pdf("After_QCDF_nodot.pdf", width = 11, height = 6)
VlnPlot(alldata_after, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()
dev.off()

### Histogram
metadata3 <- alldata_after@meta.data
metadata3$sample <- NA
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
metadata3$group[which(str_detect(metadata3$orig.ident, "^CBC_C"))] <- "CBC_C"
metadata3$group[which(str_detect(metadata3$orig.ident, "^CBC_M"))] <- "CBC_M"
alldata_after@meta.data <- metadata3

### Visualize the number UMIs/transcripts per cell
metadata3 %>% 
  ggplot(aes(color=group, x=nCount_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
#geom_vline(xintercept = 10000)

# Visualize the distribution of genes detected per cell via histogram
metadata3 %>% 
  ggplot(aes(color=group, x=nFeature_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
#geom_vline(xintercept = c(300, 400))


VlnPlot(alldata_after, group.by = "group", features = feats, pt.size = 0, ncol = 5) +
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
full.combined <- readRDS("Seurat.Integration.RDS")

metadata <- full.combined@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
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
metadata$group[which(str_detect(metadata$cells, "^CBC_C"))] <- "CBC_C"
metadata$group[which(str_detect(metadata$cells, "^CBC_M"))] <- "CBC_M"
full.combined@meta.data <- metadata

metadata$case <- NA
metadata$case[which(str_detect(metadata$cells, "^CBC_C"))] <- "Con"
metadata$case[which(str_detect(metadata$cells, "^CBC_M"))] <- "MDD"
full.combined@meta.data <- metadata

VlnPlot(full.combined, group.by = "group", features = feats, pt.size = 0.1, ncol = 5) +
  NoLegend()


View(full.combined@meta.data)

### Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=group, x=nCount_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")
  #geom_vline(xintercept = 30000)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=group, x=nFeature_RNA, fill= group)) + 
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

saveRDS(object=full.combined, file="Seurat.before_Jack_PC100.RDS")

rm(alldata_after, alldata_after1, alldata_before, CBC_C1, CBC_C2, CBC_C3, CBC_C4, CBC_C5, CBC_C6,
   CBC_M1, CBC_M2, CBC_M3, CBC_M4, CBC_M5, CBC_M6)

### Too long
full.combined <- JackStraw(object=full.combined, num.replicate=100, dims=100)
full.combined <- ScoreJackStraw(object=full.combined, dims=1:100)

saveRDS(object=full.combined, file="Seurat.after_Jack_PC100.RDS")
full.combined <- readRDS(file="Seurat.after_Jack.RDS")

# PC13
JackStrawPlot(object=full.combined, dims=1:13) + NoLegend() # PC13
JackStrawPlot(object=full.combined, dims=1:100)

DefaultAssay(full.combined) <- "integrated"
full.combined <- RunUMAP(full.combined, reduction = "pca", dims = 1:13)
full.combined <- FindNeighbors(full.combined, reduction = "pca", dims = 1:13)
full.combined <- FindClusters(full.combined, resolution = 0.5)
DimPlot(full.combined, reduction = "umap", label = TRUE, repel=F, raster=F)


count_table <- table(full.combined@meta.data$integrated_snn_res.0.5, full.combined@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell_proportion.csv")


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
        split.by= "Phase")


### tSNE (long)
full.combined <- RunTSNE(full.combined, dims=1:13)
DimPlot(full.combined, reduction="tsne", label=T, repel=T, raster=T)

# save
saveRDS(object=full.combined, file="Seurat.UMAP&tSNE(before_rm).RDS")
full.combined <- readRDS(file="Seurat.UMAP&tSNE(before_rm).RDS")

### Remove cluster (after tSNE)
full.combined <- subset(full.combined, idents=10, invert=T)
DimPlot(full.combined, reduction="umap", label=T, repel=F, raster=T)
DimPlot(full.combined, reduction="tsne", label=T, repel=F, raster=T)

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(full.combined, reduction="umap", features=metrics, min.cutoff = "q9", ncol=3)

### Clustering distribution, Gender-specific gene
DefaultAssay(full.combined) <- "RNA"
VlnPlot(object=full.combined, features=c("UTY", "USP9Y", "XIST", "TSIX"), group.by="sample", ncol=2)

saveRDS(object=full.combined, file="Seurat.UMAP_rm10_seurat.RDS")
full.combined <- readRDS("Seurat.UMAP_rm10.RDS")

### UMAP
DimPlot(object=full.combined, reduction="umap", group.by="group")
DimPlot(object=full.combined, reduction="umap", group.by="case", raster=T)
DimPlot(object=full.combined, reduction="umap", group.by="sample")

DimPlot(object=full.combined, reduction="umap", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="umap", split.by="sample", ncol=6)

DimPlot(full.combined, reduction = "tsne", group.by = "group")
DimPlot(full.combined, reduction="tsne", group.by="sample")

DimPlot(object=full.combined, reduction="tsne", split.by="group", ncol=2)
DimPlot(object=full.combined, reduction="tsne", split.by="sample", ncol=6)


# Purkinje: ITPR1, RORA, SORCS3, RYR1, GAD1/2
# Molecular layer interneuron (basket, stellate): PTPRK, PRKCD, LYPD6
# Purkinje layer interneuron: KLHL1
# Interneuron: PAX2, GAD1/2
# GPC: ATOH1 X, DCC, ROBO1, PAX6, PTCH1, NEUROD1, ZIC1, ZIC2, SEMA6A, PLXNA2
# Granule: RBFOX3, RELN, GABRA6, FAT2
# iCN: MEIS2

Neuron.Marker <- c("RBFOX3", "RELN", "FAT2", "ITPR1", "RYR1", "SORCS3", "PAX2", "GAD2", "PTPRK", "PRKCD", "LYPD6")
Non.Neuron.Marker <- c("SLC1A2", "GFAP", "AQP4", "MRC1", "CD44", "CX3CR1", "CD74", "MBP", "MOG", "PDGFRA", "CLDN5", "FLT1")
Other <- c("EOMES", "LGI2", "SLC6A5", "GRM2", "OTX2", "SLIT2", "KLHL1", "PAX2", "SORCS3", "ITPR1")
All.Marker <- c("RBFOX3", "RELN", "FAT2", "ITPR1", "RYR1", "SORCS3",
                "SLC1A2", "GFAP", "CX3CR1", "CD74", "MRC1", "CD44", "MBP", "MOG", "PDGFRA", "CLDN5", "FLT1")

### Cell type marker (Neuron, Ex, In, Astro, Micro, Endo, Oligo, OPC)
Neuron <- c("RBFOX3", "RELN", "FAT2", "ITPR1", "RYR1", "SORCS3", "PAX2", "GAD2", "PTPRK", "PRKCD", "LYPD6")
Astro <- c("SLC1A2", "GLUL", "SOX9", "AQP4", "GJA1", "NDRG2", "ALDH1A1", "ALDH1L1", "VIM")
Micro <- c("LAPTM5", "CD74", "SPI1", "TMEM119", "CX3CR1", "HEXB")
Macro <- c("CD96", "CD44", "MRC1", "CYBB", "MS4A7", "AXL")
Endo <- c("CLDN5", "VTN", "FLT1", "RAMP2", "KDR")
Oligo <- c("MBP", "MOBP", "MAG", "PLP1", "MOG",
                   "PDGFRA", "PCDH15", "OLIG2", "OLIG1")
Mix <- c("MRC1", "RYR1", "GLCE", "GRM4") #Macro, Purkinjce 2, Granule


DefaultAssay(full.combined) <- "RNA"
FeaturePlot(object=full.combined, reduction="umap", features=Neuron.Marker, label=F, min.cutoff = "q9", ncol=4)
DotPlot(full.combined, features = Purkinje.markers)+
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),)+
  labs(x=NULL, y=NULL)
VlnPlot(full.combined, features=All.Marker, stack=T, sort=F, flip=F) +
  theme(legend.position = "none")+
  theme(strip.text.x = element_text(face="italic", size=11),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(angle=270))+
  labs(x=NULL, y=NULL)



### Cell proportion
count_table <- table(full.combined@meta.data$integrated_snn_res.0.5, full.combined@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell proportion/Cell_count_table.csv")
dittoBarPlot(
  object = full.combined,
  var = "integrated_snn_res.0.5",
  group.by = "sample",
  var.labels.reorder = c(1,2,12,13,14,15,16,17,18,19,3,4,5,6,7,8,9,10,11))

### DEG & Heatmap (too long)
DefaultAssay(full.combined) <- "RNA"
seurat.markers <- FindAllMarkers(object=full.combined, only.pos=T, min.pct=0.1, logfc.threshold=0.25)
head(seurat.markers)

seurat.markers <- readRDS("Seurat.FindAllMarker_PCT0.1.RDS")

seurat.top10.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
seurat.top20.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
seurat.top100.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_log2FC)

DefaultAssay(full.combined) <- "integrated"
DoHeatmap(object=full.combined, features=seurat.top10.genes$gene, size=4) + NoLegend()

### Save
write.csv(seurat.markers, "Seurat.FindAllMarkers_PCT0.1.csv")
write.csv(seurat.top10.genes, "Seurat.top10markers_PCT0.1.csv")
write.csv(seurat.top20.genes, "Seurat.top20markers_PCT0.1.csv")
write.csv(seurat.top100.genes, "Seurat.top100markers_PCT0.1.csv")
saveRDS(object=seurat.markers, file="Seurat.FindAllMarker_PCT0.1.RDS")
seurat.markers <- readRDS("Seurat.FindAllMarker_rm10.RDS")

### Convert Cluster ID
levels(full.combined)
new.cluster.ids <- c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                     "Oligo", "Astro1", "Purk1", "Purk2", "Micro", "OPC", "Astro2", "Endo")
names(new.cluster.ids) <- levels(full.combined)
full.combined <- RenameIdents(full.combined, new.cluster.ids)
full.combined@meta.data$"Newcluster" <- as.factor(full.combined@active.ident)
DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=4, repel=T) + NoLegend()


### Cluster Dendrogram
DefaultAssay(full.combined) <- "integrated"
full.combined <- BuildClusterTree(full.combined, assay="integrated")
myPhyTree <- Tool(object=full.combined, slot = "BuildClusterTree")
#tiff(filename = "Marker/Dendrogram_4.tif", width = 4, height = 7, units = "in", res=300)
pdf("Marker/Dendrogram2.pdf", width = 3, height = 8)
ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1,
                tip.col=c("#e897c7", "#e897c7",
                          "#d9252f", "#d9252f",
                          "#3aa312", "#2e75d1",
                          "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf", "#6c43bf",
                          "#967341", "#FFB500")) + geom_tiplab()
dev.off()

levels(full.combined)
levels(full.combined) <- c("Purk1", "Purk2", "Astro1", "Astro2", "Micro", "Endo",
                           "Gran7", "Gran2", "Gran4", "Gran5", "Gran1", "Gran3", "Gran6",
                           "Oligo", "OPC")


DefaultAssay(full.combined) <- "RNA"
AllMarkers <- c("RBFOX3", "RELN", "FAT2",
                "ITPR1", "RYR1", "SORCS3", "NXPH1", 
                "SLC1A2", "GLUL", "GFAP",
                "PDGFRA", "PCDH15",
                "MBP", "MOG", "PLP1",
                "CD74", "CX3CR1", "MRC1", "CD44", "CD96",
                "CLDN5", "FLT1")
Test.Marker <- c("ALDOC", "KCTD16", "ERVMER61-1", "PHACTR1", "LRRC4C", "SCN1A")
my_cols <- c("#6c43bf", "#6c43bf", "#6c43bf", 
             "#e897c7", "#e897c7", "#e897c7", "#e897c7",
             "#d9252f", "#d9252f", "#d9252f",
             "#FFB500", "#FFB500",
             "#967341", "#967341", "#967341",
             "#3aa312", "#3aa312", "#3aa312", "#3aa312", "#3aa312",
             "#2e75d1", "#2e75d1")

AllMarkers <- c("RBFOX3", "RELN",
                "ITPR1", "RYR1", "SORCS3", "NXPH1", 
                "SLC1A2", "GFAP",
                "PDGFRA", "PCDH15",
                "MBP", "MOG",
                "CX3CR1", "MRC1", "CD44",
                "CLDN5", "FLT1")
my_cols <- c("#6c43bf", "#6c43bf",
             "#e897c7", "#e897c7", "#e897c7", "#e897c7",
             "#d9252f", "#d9252f",
             "#FFB500", "#FFB500",
             "#967341", "#967341",
             "#3aa312", "#3aa312", "#3aa312",
             "#2e75d1", "#2e75d1")

pdf("DotPlot_All.pdf", width = 8.017, height = 7)
DotPlot(full.combined,features = AllMarkers, cols=c("#e8e8e8", "#751665"), dot.scale = 7) +  ##0b488a #a10d63 #751665 #0a4d54
  theme(axis.text.x = element_text(angle = 90, size=11),
        axis.text.x.bottom = element_text(colour = my_cols, face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nexpression"),
         size = guide_legend("Percent\nexpressed"),)+
  labs(x=NULL, y=NULL)
dev.off()

tiff(filename = "VlnPlot_All1.tif", width = 10, height = 7, units = "in", res=300)
VlnPlot(full.combined, features=AllMarkers, stack=T, sort=F, flip=F) +
  theme(legend.position = "none")+
  theme(strip.text.x = element_text(face="italic", size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=270))+
  labs(x=NULL, y=NULL)
dev.off()


saveRDS(object=full.combined, file="Seurat.Final.RDS")
full.combined <- readRDS(file="Seurat.Final.RDS")
full.combined_Large <- readRDS(file="Seurat.Final.RDS")

### Large cluster
levels(full.combined_Large)
new.cluster.ids <- c("Purk", "Purk", "Astro", "Astro", "Micro", "Endo",
                     "Gran", "Gran", "Gran", "Gran", "Gran", "Gran", "Gran", "Oligo", "OPC")
names(new.cluster.ids) <- levels(full.combined_Large)
full.combined_Large <- RenameIdents(full.combined_Large, new.cluster.ids)
full.combined_Large@meta.data$"Largecluster" <- as.factor(full.combined_Large@active.ident)
View(full.combined_Large@meta.data)

DimPlot(full.combined, reduction = "umap", label = TRUE, label.size=3, repel=F, raster=T) + NoLegend()

### UMAP color
levels(full.combined)
my_cols <- c('Astro1'='#ed8b28', 'Astro2'='#e33b50',
             'OPC'='#FFB500', 'Oligo'='#967341',
             'Micro' = '#A4E804',
             'Endo' = '#2e75d1',
             'Gran1'='#c8bdf0', 'Gran2'='#918fe3', 'Gran3'='#7e7cd6', 'Gran4'='#b06ae6',
             'Gran5'='#8e6cd4', 'Gran6'='#5d58e0', 'Gran7'='#6f7fa6',
             'Purk1'='#e07b91', 'Purk2'='#ed91bf')
my_cols <- c('Gran'='#918fe3', 'Purk'='#ed91bf',
             'Astro'='#ed8b28', 'OPC'='#FFB500', 'Oligo'='#967341', 'Micro'='#A4E804', 
             'Endo'='#2e75d1')
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]


pdf("UMAP_Large.pdf", width = 7.945, height = 7)
DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = TRUE, label.size=3, repel=T, raster=T)
dev.off()

pdf("UMAP_Final.pdf", width = 7.034, height = 7)
DimPlot(full.combined, reduction = "umap", cols=my_cols2, label = TRUE, label.size=3, repel=T, raster=T) + NoLegend()
dev.off()

pdf("tSNE_Final.pdf", width = 8.11, height = 7)
DimPlot(full.combined, reduction = "tsne", cols=my_cols2, label = TRUE, label.size=3, repel=T, raster=T)
dev.off()


### Cell proportion
pdf("Cell_proportion_sample.pdf", width = 3.692, height = 5.954)
dittoBarPlot(
  object = full.combined,
  var = "Newcluster",
  group.by = "sample",
  main="",
  xlab="",
  color.panel = my_cols2,
  var.labels.reorder=c(1,2,13,12,11,3,4,5,6,7,8,9,10,14,15))+
  theme(axis.text=element_text(size=7, color="black"),
        legend.text=element_text(size=9),
        legend.key.size = unit(0.5, 'cm'))
dev.off()