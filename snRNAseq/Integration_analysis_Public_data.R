library(ggtree)
library(GO.db)
library(scRNAseq)
library(SingleCellExperiment)
library(scater)
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
library(Seurat)
library(SeuratData)
suppressMessages(require(DoubletFinder))


### Load Male data
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  assign(paste0("MC", i), Read10X(data.dir = paste0("C",i,"/outs/filtered_feature_bc_matrix/")))
}

MC13_1 <- Read10X(data.dir = "C13_1/outs/filtered_feature_bc_matrix/")
MC13_2 <- Read10X(data.dir = "C13_2/outs/filtered_feature_bc_matrix/")

for(i in 1:17){
  assign(paste0("MM", i), Read10X(data.dir = paste0("M",i,"/outs/filtered_feature_bc_matrix/")))
}

### Create Seurat object
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  assign(paste0("Male_C", i), CreateSeuratObject(counts=get(paste0("MC", i)), min.cells=3, project=paste0("Male_C", i)))
}

Male_C13_1 <- CreateSeuratObject(counts=MC13_1, min.cells=3, project="Male_C13_1")
Male_C13_2 <- CreateSeuratObject(counts=MC13_2, min.cells=3, project="Male_C13_2")

for(i in 1:17){
  assign(paste0("Male_M", i), CreateSeuratObject(counts=get(paste0("MM", i)), min.cells=3, project=paste0("Male_M", i)))
}


### Load Female data
for(i in c(1:18)){
  assign(paste0("FC", i), Read10X(data.dir = paste0("C",i,"/outs/filtered_feature_bc_matrix/")))
}

for(i in 1:20){
  assign(paste0("FM", i), Read10X(data.dir = paste0("M",i,"/outs/filtered_feature_bc_matrix/")))
}

### Create Seurat object
for(i in c(1:18)){
  assign(paste0("Female_C", i), CreateSeuratObject(counts=get(paste0("FC", i)), min.cells=3, project=paste0("Female_C", i)))
}

for(i in 1:20){
  assign(paste0("Female_M", i), CreateSeuratObject(counts=get(paste0("FM", i)), min.cells=3, project=paste0("Female_M", i)))
}

### remove matrix
rm(MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8, MC10,
   MC11, MC12, MC13_1, MC13_2, MC15, MC16, MC17)
rm(MM1, MM2, MM3, MM4, MM5, MM6, MM7, MM8, MM9, MM10,
   MM11, MM12, MM13, MM14, MM15, MM16, MM17)
rm(FC1, FC2, FC3, FC4, FC5, FC6, FC7, FC8, FC9, FC10,
   FC11, FC12, FC13, FC14, FC15, FC16, FC17, FC18)
rm(FM1, FM2, FM3, FM4, FM5, FM6, FM7, FM8, FM9, FM10,
   FM11, FM12, FM13, FM14, FM15, FM16, FM17, FM18, FM19, FM20)

### Load My data
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
mito_genes <- rownames(Male_C1)[grep("^MT-", rownames(Male_C1))]
Male_C1[["percent.mt"]] <- PercentageFeatureSet(Male_C1, pattern="^MT-")
Male_C2[["percent.mt"]] <- PercentageFeatureSet(Male_C2, pattern="^MT-")
Male_C3[["percent.mt"]] <- PercentageFeatureSet(Male_C3, pattern="^MT-")
Male_C4[["percent.mt"]] <- PercentageFeatureSet(Male_C4, pattern="^MT-")
Male_C5[["percent.mt"]] <- PercentageFeatureSet(Male_C5, pattern="^MT-")
Male_C6[["percent.mt"]] <- PercentageFeatureSet(Male_C6, pattern="^MT-")
Male_C7[["percent.mt"]] <- PercentageFeatureSet(Male_C7, pattern="^MT-")
Male_C8[["percent.mt"]] <- PercentageFeatureSet(Male_C8, pattern="^MT-")
#Male_C9[["percent.mt"]] <- PercentageFeatureSet(Male_C9, pattern="^MT-")
Male_C10[["percent.mt"]] <- PercentageFeatureSet(Male_C10, pattern="^MT-")
Male_C11[["percent.mt"]] <- PercentageFeatureSet(Male_C11, pattern="^MT-")
Male_C12[["percent.mt"]] <- PercentageFeatureSet(Male_C12, pattern="^MT-")
Male_C13_1[["percent.mt"]] <- PercentageFeatureSet(Male_C13_1, pattern="^MT-")
Male_C13_2[["percent.mt"]] <- PercentageFeatureSet(Male_C13_2, pattern="^MT-")
#Male_C14_1[["percent.mt"]] <- PercentageFeatureSet(Male_C14_1, pattern="^MT-")
#Male_C14_2[["percent.mt"]] <- PercentageFeatureSet(Male_C14_2, pattern="^MT-")
Male_C15[["percent.mt"]] <- PercentageFeatureSet(Male_C15, pattern="^MT-")
Male_C16[["percent.mt"]] <- PercentageFeatureSet(Male_C16, pattern="^MT-")
Male_C17[["percent.mt"]] <- PercentageFeatureSet(Male_C17, pattern="^MT-")

Male_M1[["percent.mt"]] <- PercentageFeatureSet(Male_M1, pattern="^MT-")
Male_M2[["percent.mt"]] <- PercentageFeatureSet(Male_M2, pattern="^MT-")
Male_M3[["percent.mt"]] <- PercentageFeatureSet(Male_M3, pattern="^MT-")
Male_M4[["percent.mt"]] <- PercentageFeatureSet(Male_M4, pattern="^MT-")
Male_M5[["percent.mt"]] <- PercentageFeatureSet(Male_M5, pattern="^MT-")
Male_M6[["percent.mt"]] <- PercentageFeatureSet(Male_M6, pattern="^MT-")
Male_M7[["percent.mt"]] <- PercentageFeatureSet(Male_M7, pattern="^MT-")
Male_M8[["percent.mt"]] <- PercentageFeatureSet(Male_M8, pattern="^MT-")
Male_M9[["percent.mt"]] <- PercentageFeatureSet(Male_M9, pattern="^MT-")
Male_M10[["percent.mt"]] <- PercentageFeatureSet(Male_M10, pattern="^MT-")
Male_M11[["percent.mt"]] <- PercentageFeatureSet(Male_M11, pattern="^MT-")
Male_M12[["percent.mt"]] <- PercentageFeatureSet(Male_M12, pattern="^MT-")
Male_M13[["percent.mt"]] <- PercentageFeatureSet(Male_M13, pattern="^MT-")
Male_M14[["percent.mt"]] <- PercentageFeatureSet(Male_M14, pattern="^MT-")
Male_M15[["percent.mt"]] <- PercentageFeatureSet(Male_M15, pattern="^MT-")
Male_M16[["percent.mt"]] <- PercentageFeatureSet(Male_M16, pattern="^MT-")
Male_M17[["percent.mt"]] <- PercentageFeatureSet(Male_M17, pattern="^MT-")

Female_C1[["percent.mt"]] <- PercentageFeatureSet(Female_C1, pattern="^MT-")
Female_C2[["percent.mt"]] <- PercentageFeatureSet(Female_C2, pattern="^MT-")
Female_C3[["percent.mt"]] <- PercentageFeatureSet(Female_C3, pattern="^MT-")
Female_C4[["percent.mt"]] <- PercentageFeatureSet(Female_C4, pattern="^MT-")
Female_C5[["percent.mt"]] <- PercentageFeatureSet(Female_C5, pattern="^MT-")
Female_C6[["percent.mt"]] <- PercentageFeatureSet(Female_C6, pattern="^MT-")
Female_C7[["percent.mt"]] <- PercentageFeatureSet(Female_C7, pattern="^MT-")
Female_C8[["percent.mt"]] <- PercentageFeatureSet(Female_C8, pattern="^MT-")
Female_C9[["percent.mt"]] <- PercentageFeatureSet(Female_C9, pattern="^MT-")
Female_C10[["percent.mt"]] <- PercentageFeatureSet(Female_C10, pattern="^MT-")
Female_C11[["percent.mt"]] <- PercentageFeatureSet(Female_C11, pattern="^MT-")
Female_C12[["percent.mt"]] <- PercentageFeatureSet(Female_C12, pattern="^MT-")
Female_C13[["percent.mt"]] <- PercentageFeatureSet(Female_C13, pattern="^MT-")
Female_C14[["percent.mt"]] <- PercentageFeatureSet(Female_C14, pattern="^MT-")
Female_C15[["percent.mt"]] <- PercentageFeatureSet(Female_C15, pattern="^MT-")
Female_C16[["percent.mt"]] <- PercentageFeatureSet(Female_C16, pattern="^MT-")
Female_C17[["percent.mt"]] <- PercentageFeatureSet(Female_C17, pattern="^MT-")
Female_C18[["percent.mt"]] <- PercentageFeatureSet(Female_C18, pattern="^MT-")

Female_M1[["percent.mt"]] <- PercentageFeatureSet(Female_M1, pattern="^MT-")
Female_M2[["percent.mt"]] <- PercentageFeatureSet(Female_M2, pattern="^MT-")
Female_M3[["percent.mt"]] <- PercentageFeatureSet(Female_M3, pattern="^MT-")
Female_M4[["percent.mt"]] <- PercentageFeatureSet(Female_M4, pattern="^MT-")
Female_M5[["percent.mt"]] <- PercentageFeatureSet(Female_M5, pattern="^MT-")
Female_M6[["percent.mt"]] <- PercentageFeatureSet(Female_M6, pattern="^MT-")
Female_M7[["percent.mt"]] <- PercentageFeatureSet(Female_M7, pattern="^MT-")
Female_M8[["percent.mt"]] <- PercentageFeatureSet(Female_M8, pattern="^MT-")
Female_M9[["percent.mt"]] <- PercentageFeatureSet(Female_M9, pattern="^MT-")
Female_M10[["percent.mt"]] <- PercentageFeatureSet(Female_M10, pattern="^MT-")
Female_M11[["percent.mt"]] <- PercentageFeatureSet(Female_M11, pattern="^MT-")
Female_M12[["percent.mt"]] <- PercentageFeatureSet(Female_M12, pattern="^MT-")
Female_M13[["percent.mt"]] <- PercentageFeatureSet(Female_M13, pattern="^MT-")
Female_M14[["percent.mt"]] <- PercentageFeatureSet(Female_M14, pattern="^MT-")
Female_M15[["percent.mt"]] <- PercentageFeatureSet(Female_M15, pattern="^MT-")
Female_M16[["percent.mt"]] <- PercentageFeatureSet(Female_M16, pattern="^MT-")
Female_M17[["percent.mt"]] <- PercentageFeatureSet(Female_M17, pattern="^MT-")
Female_M18[["percent.mt"]] <- PercentageFeatureSet(Female_M18, pattern="^MT-")
Female_M19[["percent.mt"]] <- PercentageFeatureSet(Female_M19, pattern="^MT-")
Female_M20[["percent.mt"]] <- PercentageFeatureSet(Female_M20, pattern="^MT-")

### Ribosomal genes
ribo_genes <- rownames(Male_C1)[grep("^RP[SL]", rownames(Male_C1))]
Male_C1[["percent.ribo"]] <- PercentageFeatureSet(Male_C1, pattern="^RP[SL]")
Male_C2[["percent.ribo"]] <- PercentageFeatureSet(Male_C2, pattern="^RP[SL]")
Male_C3[["percent.ribo"]] <- PercentageFeatureSet(Male_C3, pattern="^RP[SL]")
Male_C4[["percent.ribo"]] <- PercentageFeatureSet(Male_C4, pattern="^RP[SL]")
Male_C5[["percent.ribo"]] <- PercentageFeatureSet(Male_C5, pattern="^RP[SL]")
Male_C6[["percent.ribo"]] <- PercentageFeatureSet(Male_C6, pattern="^RP[SL]")
Male_C7[["percent.ribo"]] <- PercentageFeatureSet(Male_C7, pattern="^RP[SL]")
Male_C8[["percent.ribo"]] <- PercentageFeatureSet(Male_C8, pattern="^RP[SL]")
#Male_C9[["percent.ribo"]] <- PercentageFeatureSet(Male_C9, pattern="^RP[SL]")
Male_C10[["percent.ribo"]] <- PercentageFeatureSet(Male_C10, pattern="^RP[SL]")
Male_C11[["percent.ribo"]] <- PercentageFeatureSet(Male_C11, pattern="^RP[SL]")
Male_C12[["percent.ribo"]] <- PercentageFeatureSet(Male_C12, pattern="^RP[SL]")
Male_C13_1[["percent.ribo"]] <- PercentageFeatureSet(Male_C13_1, pattern="^RP[SL]")
Male_C13_2[["percent.ribo"]] <- PercentageFeatureSet(Male_C13_2, pattern="^RP[SL]")
#Male_C14_1[["percent.ribo"]] <- PercentageFeatureSet(Male_C14_1, pattern="^RP[SL]")
#Male_C14_2[["percent.ribo"]] <- PercentageFeatureSet(Male_C14_2, pattern="^RP[SL]")
Male_C15[["percent.ribo"]] <- PercentageFeatureSet(Male_C15, pattern="^RP[SL]")
Male_C16[["percent.ribo"]] <- PercentageFeatureSet(Male_C16, pattern="^RP[SL]")
Male_C17[["percent.ribo"]] <- PercentageFeatureSet(Male_C17, pattern="^RP[SL]")

Male_M1[["percent.ribo"]] <- PercentageFeatureSet(Male_M1, pattern="^RP[SL]")
Male_M2[["percent.ribo"]] <- PercentageFeatureSet(Male_M2, pattern="^RP[SL]")
Male_M3[["percent.ribo"]] <- PercentageFeatureSet(Male_M3, pattern="^RP[SL]")
Male_M4[["percent.ribo"]] <- PercentageFeatureSet(Male_M4, pattern="^RP[SL]")
Male_M5[["percent.ribo"]] <- PercentageFeatureSet(Male_M5, pattern="^RP[SL]")
Male_M6[["percent.ribo"]] <- PercentageFeatureSet(Male_M6, pattern="^RP[SL]")
Male_M7[["percent.ribo"]] <- PercentageFeatureSet(Male_M7, pattern="^RP[SL]")
Male_M8[["percent.ribo"]] <- PercentageFeatureSet(Male_M8, pattern="^RP[SL]")
Male_M9[["percent.ribo"]] <- PercentageFeatureSet(Male_M9, pattern="^RP[SL]")
Male_M10[["percent.ribo"]] <- PercentageFeatureSet(Male_M10, pattern="^RP[SL]")
Male_M11[["percent.ribo"]] <- PercentageFeatureSet(Male_M11, pattern="^RP[SL]")
Male_M12[["percent.ribo"]] <- PercentageFeatureSet(Male_M12, pattern="^RP[SL]")
Male_M13[["percent.ribo"]] <- PercentageFeatureSet(Male_M13, pattern="^RP[SL]")
Male_M14[["percent.ribo"]] <- PercentageFeatureSet(Male_M14, pattern="^RP[SL]")
Male_M15[["percent.ribo"]] <- PercentageFeatureSet(Male_M15, pattern="^RP[SL]")
Male_M16[["percent.ribo"]] <- PercentageFeatureSet(Male_M16, pattern="^RP[SL]")
Male_M17[["percent.ribo"]] <- PercentageFeatureSet(Male_M17, pattern="^RP[SL]")

Female_C1[["percent.ribo"]] <- PercentageFeatureSet(Female_C1, pattern="^RP[SL]")
Female_C2[["percent.ribo"]] <- PercentageFeatureSet(Female_C2, pattern="^RP[SL]")
Female_C3[["percent.ribo"]] <- PercentageFeatureSet(Female_C3, pattern="^RP[SL]")
Female_C4[["percent.ribo"]] <- PercentageFeatureSet(Female_C4, pattern="^RP[SL]")
Female_C5[["percent.ribo"]] <- PercentageFeatureSet(Female_C5, pattern="^RP[SL]")
Female_C6[["percent.ribo"]] <- PercentageFeatureSet(Female_C6, pattern="^RP[SL]")
Female_C7[["percent.ribo"]] <- PercentageFeatureSet(Female_C7, pattern="^RP[SL]")
Female_C8[["percent.ribo"]] <- PercentageFeatureSet(Female_C8, pattern="^RP[SL]")
Female_C9[["percent.ribo"]] <- PercentageFeatureSet(Female_C9, pattern="^RP[SL]")
Female_C10[["percent.ribo"]] <- PercentageFeatureSet(Female_C10, pattern="^RP[SL]")
Female_C11[["percent.ribo"]] <- PercentageFeatureSet(Female_C11, pattern="^RP[SL]")
Female_C12[["percent.ribo"]] <- PercentageFeatureSet(Female_C12, pattern="^RP[SL]")
Female_C13[["percent.ribo"]] <- PercentageFeatureSet(Female_C13, pattern="^RP[SL]")
Female_C14[["percent.ribo"]] <- PercentageFeatureSet(Female_C14, pattern="^RP[SL]")
Female_C15[["percent.ribo"]] <- PercentageFeatureSet(Female_C15, pattern="^RP[SL]")
Female_C16[["percent.ribo"]] <- PercentageFeatureSet(Female_C16, pattern="^RP[SL]")
Female_C17[["percent.ribo"]] <- PercentageFeatureSet(Female_C17, pattern="^RP[SL]")
Female_C18[["percent.ribo"]] <- PercentageFeatureSet(Female_C18, pattern="^RP[SL]")

Female_M1[["percent.ribo"]] <- PercentageFeatureSet(Female_M1, pattern="^RP[SL]")
Female_M2[["percent.ribo"]] <- PercentageFeatureSet(Female_M2, pattern="^RP[SL]")
Female_M3[["percent.ribo"]] <- PercentageFeatureSet(Female_M3, pattern="^RP[SL]")
Female_M4[["percent.ribo"]] <- PercentageFeatureSet(Female_M4, pattern="^RP[SL]")
Female_M5[["percent.ribo"]] <- PercentageFeatureSet(Female_M5, pattern="^RP[SL]")
Female_M6[["percent.ribo"]] <- PercentageFeatureSet(Female_M6, pattern="^RP[SL]")
Female_M7[["percent.ribo"]] <- PercentageFeatureSet(Female_M7, pattern="^RP[SL]")
Female_M8[["percent.ribo"]] <- PercentageFeatureSet(Female_M8, pattern="^RP[SL]")
Female_M9[["percent.ribo"]] <- PercentageFeatureSet(Female_M9, pattern="^RP[SL]")
Female_M10[["percent.ribo"]] <- PercentageFeatureSet(Female_M10, pattern="^RP[SL]")
Female_M11[["percent.ribo"]] <- PercentageFeatureSet(Female_M11, pattern="^RP[SL]")
Female_M12[["percent.ribo"]] <- PercentageFeatureSet(Female_M12, pattern="^RP[SL]")
Female_M13[["percent.ribo"]] <- PercentageFeatureSet(Female_M13, pattern="^RP[SL]")
Female_M14[["percent.ribo"]] <- PercentageFeatureSet(Female_M14, pattern="^RP[SL]")
Female_M15[["percent.ribo"]] <- PercentageFeatureSet(Female_M15, pattern="^RP[SL]")
Female_M16[["percent.ribo"]] <- PercentageFeatureSet(Female_M16, pattern="^RP[SL]")
Female_M17[["percent.ribo"]] <- PercentageFeatureSet(Female_M17, pattern="^RP[SL]")
Female_M18[["percent.ribo"]] <- PercentageFeatureSet(Female_M18, pattern="^RP[SL]")
Female_M19[["percent.ribo"]] <- PercentageFeatureSet(Female_M19, pattern="^RP[SL]")
Female_M20[["percent.ribo"]] <- PercentageFeatureSet(Female_M20, pattern="^RP[SL]")

### Hemoglobin genes
hemo_genes <- rownames(Male_C1)[grep("^HB[^(P)]", rownames(Male_C1))]
Male_C1[["percent.hemo"]] <- PercentageFeatureSet(Male_C1, pattern="^HB[^(P)]")
Male_C2[["percent.hemo"]] <- PercentageFeatureSet(Male_C2, pattern="^HB[^(P)]")
Male_C3[["percent.hemo"]] <- PercentageFeatureSet(Male_C3, pattern="^HB[^(P)]")
Male_C4[["percent.hemo"]] <- PercentageFeatureSet(Male_C4, pattern="^HB[^(P)]")
Male_C5[["percent.hemo"]] <- PercentageFeatureSet(Male_C5, pattern="^HB[^(P)]")
Male_C6[["percent.hemo"]] <- PercentageFeatureSet(Male_C6, pattern="^HB[^(P)]")
Male_C7[["percent.hemo"]] <- PercentageFeatureSet(Male_C7, pattern="^HB[^(P)]")
Male_C8[["percent.hemo"]] <- PercentageFeatureSet(Male_C8, pattern="^HB[^(P)]")
#Male_C9[["percent.hemo"]] <- PercentageFeatureSet(Male_C9, pattern="^HB[^(P)]")
Male_C10[["percent.hemo"]] <- PercentageFeatureSet(Male_C10, pattern="^HB[^(P)]")
Male_C11[["percent.hemo"]] <- PercentageFeatureSet(Male_C11, pattern="^HB[^(P)]")
Male_C12[["percent.hemo"]] <- PercentageFeatureSet(Male_C12, pattern="^HB[^(P)]")
Male_C13_1[["percent.hemo"]] <- PercentageFeatureSet(Male_C13_1, pattern="^HB[^(P)]")
Male_C13_2[["percent.hemo"]] <- PercentageFeatureSet(Male_C13_2, pattern="^HB[^(P)]")
#Male_C14_1[["percent.hemo"]] <- PercentageFeatureSet(Male_C14_1, pattern="^HB[^(P)]")
#Male_C14_2[["percent.hemo"]] <- PercentageFeatureSet(Male_C14_2, pattern="^HB[^(P)]")
Male_C15[["percent.hemo"]] <- PercentageFeatureSet(Male_C15, pattern="^HB[^(P)]")
Male_C16[["percent.hemo"]] <- PercentageFeatureSet(Male_C16, pattern="^HB[^(P)]")
Male_C17[["percent.hemo"]] <- PercentageFeatureSet(Male_C17, pattern="^HB[^(P)]")

Male_M1[["percent.hemo"]] <- PercentageFeatureSet(Male_M1, pattern="^HB[^(P)]")
Male_M2[["percent.hemo"]] <- PercentageFeatureSet(Male_M2, pattern="^HB[^(P)]")
Male_M3[["percent.hemo"]] <- PercentageFeatureSet(Male_M3, pattern="^HB[^(P)]")
Male_M4[["percent.hemo"]] <- PercentageFeatureSet(Male_M4, pattern="^HB[^(P)]")
Male_M5[["percent.hemo"]] <- PercentageFeatureSet(Male_M5, pattern="^HB[^(P)]")
Male_M6[["percent.hemo"]] <- PercentageFeatureSet(Male_M6, pattern="^HB[^(P)]")
Male_M7[["percent.hemo"]] <- PercentageFeatureSet(Male_M7, pattern="^HB[^(P)]")
Male_M8[["percent.hemo"]] <- PercentageFeatureSet(Male_M8, pattern="^HB[^(P)]")
Male_M9[["percent.hemo"]] <- PercentageFeatureSet(Male_M9, pattern="^HB[^(P)]")
Male_M10[["percent.hemo"]] <- PercentageFeatureSet(Male_M10, pattern="^HB[^(P)]")
Male_M11[["percent.hemo"]] <- PercentageFeatureSet(Male_M11, pattern="^HB[^(P)]")
Male_M12[["percent.hemo"]] <- PercentageFeatureSet(Male_M12, pattern="^HB[^(P)]")
Male_M13[["percent.hemo"]] <- PercentageFeatureSet(Male_M13, pattern="^HB[^(P)]")
Male_M14[["percent.hemo"]] <- PercentageFeatureSet(Male_M14, pattern="^HB[^(P)]")
Male_M15[["percent.hemo"]] <- PercentageFeatureSet(Male_M15, pattern="^HB[^(P)]")
Male_M16[["percent.hemo"]] <- PercentageFeatureSet(Male_M16, pattern="^HB[^(P)]")
Male_M17[["percent.hemo"]] <- PercentageFeatureSet(Male_M17, pattern="^HB[^(P)]")

Female_C1[["percent.hemo"]] <- PercentageFeatureSet(Female_C1, pattern="^HB[^(P)]")
Female_C2[["percent.hemo"]] <- PercentageFeatureSet(Female_C2, pattern="^HB[^(P)]")
Female_C3[["percent.hemo"]] <- PercentageFeatureSet(Female_C3, pattern="^HB[^(P)]")
Female_C4[["percent.hemo"]] <- PercentageFeatureSet(Female_C4, pattern="^HB[^(P)]")
Female_C5[["percent.hemo"]] <- PercentageFeatureSet(Female_C5, pattern="^HB[^(P)]")
Female_C6[["percent.hemo"]] <- PercentageFeatureSet(Female_C6, pattern="^HB[^(P)]")
Female_C7[["percent.hemo"]] <- PercentageFeatureSet(Female_C7, pattern="^HB[^(P)]")
Female_C8[["percent.hemo"]] <- PercentageFeatureSet(Female_C8, pattern="^HB[^(P)]")
Female_C9[["percent.hemo"]] <- PercentageFeatureSet(Female_C9, pattern="^HB[^(P)]")
Female_C10[["percent.hemo"]] <- PercentageFeatureSet(Female_C10, pattern="^HB[^(P)]")
Female_C11[["percent.hemo"]] <- PercentageFeatureSet(Female_C11, pattern="^HB[^(P)]")
Female_C12[["percent.hemo"]] <- PercentageFeatureSet(Female_C12, pattern="^HB[^(P)]")
Female_C13[["percent.hemo"]] <- PercentageFeatureSet(Female_C13, pattern="^HB[^(P)]")
Female_C14[["percent.hemo"]] <- PercentageFeatureSet(Female_C14, pattern="^HB[^(P)]")
Female_C15[["percent.hemo"]] <- PercentageFeatureSet(Female_C15, pattern="^HB[^(P)]")
Female_C16[["percent.hemo"]] <- PercentageFeatureSet(Female_C16, pattern="^HB[^(P)]")
Female_C17[["percent.hemo"]] <- PercentageFeatureSet(Female_C17, pattern="^HB[^(P)]")
Female_C18[["percent.hemo"]] <- PercentageFeatureSet(Female_C18, pattern="^HB[^(P)]")

Female_M1[["percent.hemo"]] <- PercentageFeatureSet(Female_M1, pattern="^HB[^(P)]")
Female_M2[["percent.hemo"]] <- PercentageFeatureSet(Female_M2, pattern="^HB[^(P)]")
Female_M3[["percent.hemo"]] <- PercentageFeatureSet(Female_M3, pattern="^HB[^(P)]")
Female_M4[["percent.hemo"]] <- PercentageFeatureSet(Female_M4, pattern="^HB[^(P)]")
Female_M5[["percent.hemo"]] <- PercentageFeatureSet(Female_M5, pattern="^HB[^(P)]")
Female_M6[["percent.hemo"]] <- PercentageFeatureSet(Female_M6, pattern="^HB[^(P)]")
Female_M7[["percent.hemo"]] <- PercentageFeatureSet(Female_M7, pattern="^HB[^(P)]")
Female_M8[["percent.hemo"]] <- PercentageFeatureSet(Female_M8, pattern="^HB[^(P)]")
Female_M9[["percent.hemo"]] <- PercentageFeatureSet(Female_M9, pattern="^HB[^(P)]")
Female_M10[["percent.hemo"]] <- PercentageFeatureSet(Female_M10, pattern="^HB[^(P)]")
Female_M11[["percent.hemo"]] <- PercentageFeatureSet(Female_M11, pattern="^HB[^(P)]")
Female_M12[["percent.hemo"]] <- PercentageFeatureSet(Female_M12, pattern="^HB[^(P)]")
Female_M13[["percent.hemo"]] <- PercentageFeatureSet(Female_M13, pattern="^HB[^(P)]")
Female_M14[["percent.hemo"]] <- PercentageFeatureSet(Female_M14, pattern="^HB[^(P)]")
Female_M15[["percent.hemo"]] <- PercentageFeatureSet(Female_M15, pattern="^HB[^(P)]")
Female_M16[["percent.hemo"]] <- PercentageFeatureSet(Female_M16, pattern="^HB[^(P)]")
Female_M17[["percent.hemo"]] <- PercentageFeatureSet(Female_M17, pattern="^HB[^(P)]")
Female_M18[["percent.hemo"]] <- PercentageFeatureSet(Female_M18, pattern="^HB[^(P)]")
Female_M19[["percent.hemo"]] <- PercentageFeatureSet(Female_M19, pattern="^HB[^(P)]")
Female_M20[["percent.hemo"]] <- PercentageFeatureSet(Female_M20, pattern="^HB[^(P)]")

### Add cell id
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  assign(paste0("Male_C", i), RenameCells(object=get(paste0("Male_C", i)), add.cell.id=paste0("Male_C", i)))
}
Male_C13_1 <- RenameCells(object=Male_C13_1, add.cell.id="Male_C13_1")
Male_C13_2 <- RenameCells(object=Male_C13_2, add.cell.id="Male_C13_2")

for(i in 1:17){
  assign(paste0("Male_M", i), RenameCells(object=get(paste0("Male_M", i)), add.cell.id=paste0("Male_M", i)))
}

for(i in c(1:18)){
  assign(paste0("Female_C", i), RenameCells(object=get(paste0("Female_C", i)), add.cell.id=paste0("Female_C", i)))
}

for(i in 1:20){
  assign(paste0("Female_M", i), RenameCells(object=get(paste0("Female_M", i)), add.cell.id=paste0("Female_M", i)))
}

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

### merge dataset for plot
alldata_before <- merge(Male_C1, c(Male_C2, Male_C3, Male_C4, Male_C5, Male_C6, Male_C7, Male_C8, Male_C10,
                                   Male_C11, Male_C12, Male_C13_1, Male_C13_2, Male_C15, Male_C16, Male_C17,
                                   Male_M1, Male_M2, Male_M3, Male_M4, Male_M5, Male_M6, Male_M7, Male_M8, Male_M9, Male_M10,
                                   Male_M11, Male_M12, Male_M13, Male_M14, Male_M15, Male_M16, Male_M17,
                                   Female_C1, Female_C2, Female_C3, Female_C4, Female_C5, Female_C6, Female_C7, Female_C8, Female_C9, Female_C10,
                                   Female_C11, Female_C12, Female_C13, Female_C14, Female_C15, Female_C16, Female_C17, Female_C18,
                                   Female_M1, Female_M2, Female_M3, Female_M4, Female_M5, Female_M6, Female_M7, Female_M8, Female_M9, Female_M10,
                                   Female_M11, Female_M12, Female_M13, Female_M14, Female_M15, Female_M16, Female_M17, Female_M18, Female_M19, Female_M20,
                                   dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                                   dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6))
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
feats1 <- c("nFeature_RNA", "nCount_RNA")
feats2 <- c("percent.mt")
feats3 <- c( "percent.ribo", "percent.hemo")

tiff(filename = "QC/Before_QC1.tif", width = 10, height = 15, units = "in", res=200)
VlnPlot(alldata_before, group.by = "orig.ident", features = feats3, pt.size = 0.1, ncol = 1) +
  NoLegend()
dev.off()
FeatureScatter(alldata_before, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) + NoLegend()
FeatureScatter(alldata_before, "nCount_RNA", "percent.mt", group.by = "orig.ident", pt.size = 0.5) + NoLegend()

### Histogram
metadata1 <- alldata_before@meta.data
metadata1$sample <- NA
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  metadata1$sample[which(str_detect(metadata1$orig.ident, paste0("^Male_C", i)))] <- paste0("Male_C", i)
}
metadata1$sample[which(str_detect(metadata1$orig.ident, "^Male_C13_1"))] <- "Male_C13"
metadata1$sample[which(str_detect(metadata1$orig.ident, "^Male_C13_2"))] <- "Male_C13"

for(i in 1:17){
  metadata1$sample[which(str_detect(metadata1$orig.ident, paste0("^Male_M", i)))] <- paste0("Male_M", i)
}

for(i in c(1:18)){
  metadata1$sample[which(str_detect(metadata1$orig.ident, paste0("^Female_C", i)))] <- paste0("Female_C", i)
}

for(i in 1:20){
  metadata1$sample[which(str_detect(metadata1$orig.ident, paste0("^Female_M", i)))] <- paste0("Female_M", i)
}

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
metadata1$group[which(str_detect(metadata1$orig.ident, "^Male_C"))] <- "Male_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^Male_M"))] <- "Male_M"
metadata1$group[which(str_detect(metadata1$orig.ident, "^Female_C"))] <- "Female_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^Female_M"))] <- "Female_M"
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "dlPFC_C"
metadata1$group[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "dlPFC_M"
alldata_before@meta.data <- metadata1

metadata1$case <- NA
metadata1$case[which(str_detect(metadata1$orig.ident, "^Male_C"))] <- "Control"
metadata1$case[which(str_detect(metadata1$orig.ident, "^Male_M"))] <- "MDD"
metadata1$case[which(str_detect(metadata1$orig.ident, "^Female_C"))] <- "Control"
metadata1$case[which(str_detect(metadata1$orig.ident, "^Female_M"))] <- "MDD"
metadata1$case[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "Control"
metadata1$case[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "MDD"
alldata_before@meta.data <- metadata1

metadata1$batch <- NA
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C12"))] <- "MB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M9"))] <- "MB1"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C11"))] <- "MB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M2"))] <- "MB2"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C2"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C5"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C6"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C14_1"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M4"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M6"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M8"))] <- "MB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M16"))] <- "MB3"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C4"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C7"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C10"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C13_2"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M5"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M7"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M14"))] <- "MB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M15"))] <- "MB4"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C3"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C14_2"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C15"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C17"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M1"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M10"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M12"))] <- "MB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M13"))] <- "MB5"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C1"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C8"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C9"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C13_1"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_C16"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M3"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M11"))] <- "MB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Male_M17"))] <- "MB6"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C6"))] <- "FB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M8"))] <- "FB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M15"))] <- "FB1"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C7"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C8"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C9"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C12"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M16"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M17"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M19"))] <- "FB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M20"))] <- "FB2"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C2"))] <- "FB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C11"))] <- "FB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C13"))] <- "FB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M3"))] <- "FB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M4"))] <- "FB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M5"))] <- "FB3"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C1"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C3"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C4"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C5"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M2"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M12"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M13"))] <- "FB4"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M14"))] <- "FB4"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C14"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C15"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C16"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C17"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M1"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M6"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M7"))] <- "FB5"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M9"))] <- "FB5"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C10"))] <- "FB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_C18"))] <- "FB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M10"))] <- "FB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M11"))] <- "FB6"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^Female_M18"))] <- "FB6"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M2"))] <- "MyB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M3"))] <- "MyB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M4"))] <- "MyB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M5"))] <- "MyB1"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M6"))] <- "MyB1"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C2"))] <- "MyB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C3"))] <- "MyB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C4"))] <- "MyB2"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C5"))] <- "MyB2"

metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C1"))] <- "MyB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_C6"))] <- "MyB3"
metadata1$batch[which(str_detect(metadata1$orig.ident, "^dlPFC_M1"))] <- "MyB3"
alldata_before@meta.data <- metadata1

metadata1$chemistry <- NA
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Male_C"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Male_M"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "v3"

metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C1"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C2"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C3"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C4"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C5"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C10"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C11"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C13"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C14"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C15"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C16"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C17"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C18"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M1"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M2"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M3"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M4"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M5"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M6"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M7"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M9"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M10"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M11"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M12"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M13"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M14"))] <- "v3"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M18"))] <- "v3"

metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C6"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C7"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C8"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C9"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_C12"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M8"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M15"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M16"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M17"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M19"))] <- "v2"
metadata1$chemistry[which(str_detect(metadata1$orig.ident, "^Female_M20"))] <- "v2"
alldata_before@meta.data <- metadata1

metadata1$dataset <- NA
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^Male_C"))] <- "Male"
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^Male_M"))] <- "Male"
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^Female_C"))] <- "Female"
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^Female_M"))] <- "Female"
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^dlPFC_C"))] <- "My"
metadata1$dataset[which(str_detect(metadata1$orig.ident, "^dlPFC_M"))] <- "My"
alldata_before@meta.data <- metadata1


### Visualize the number UMIs/transcripts per cell
metadata1 %>% 
  ggplot(aes(color=dataset, x=nCount_RNA, fill= dataset)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = c(9000, 12000))

# Visualize the distribution of genes detected per cell via histogram
metadata1 %>% 
  ggplot(aes(color=dataset, x=nFeature_RNA, fill= dataset)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(250, 350))

plot1 <- RidgePlot(alldata_before, features = "nFeature_RNA", group.by = "sample") +
  ylab(NULL) +
  geom_vline(xintercept = c(250, 350)) +
  NoLegend()
plot1

plot1 <- RidgePlot(alldata_before, features = "nCount_RNA", group.by = "sample") +
  ylab(NULL) +
  geom_vline(xintercept = c(9000,12000)) +
  NoLegend()
plot1

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")

VlnPlot(alldata_before, group.by = "sample", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()
VlnPlot(alldata_before, group.by = "sample", features = "percent.mt", pt.size = 0) +
  geom_hline(yintercept = 5) +
  NoLegend()

# save
saveRDS(object=alldata_before, file="Alldata_before_QC.RDS")

### Check top 5% nCount by quantile function
#quantile(dlPFC_C2$nCount_RNA,probs = c(0.75, 0.95, 0.995))
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  assign(paste0("Male_C", i), subset(get(paste0("Male_C", i)), subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5))
}

Male_C13_1 <- subset(Male_C13_1, subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5)
Male_C13_2 <- subset(Male_C13_2, subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5)

for(i in 1:17){
  assign(paste0("Male_M", i), subset(get(paste0("Male_M", i)), subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5))
}

for(i in c(1:18)){
  assign(paste0("Female_C", i), subset(get(paste0("Female_C", i)), subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5))
}

for(i in 1:20){
  assign(paste0("Female_M", i), subset(get(paste0("Female_M", i)), subset=nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5))
}

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

#rm Nagy_C9, C14
alldata_after1 <- merge(Male_C1, c(Male_C2, Male_C3, Male_C4, Male_C5, Male_C6, Male_C7, Male_C8, Male_C10,
                                   Male_C11, Male_C12, Male_C13_1, Male_C13_2, Male_C15, Male_C16, Male_C17,
                                   Male_M1, Male_M2, Male_M3, Male_M4, Male_M5, Male_M6, Male_M7, Male_M8, Male_M9, Male_M10,
                                   Male_M11, Male_M12, Male_M13, Male_M14, Male_M15, Male_M16, Male_M17,
                                   Female_C1, Female_C2, Female_C3, Female_C4, Female_C5, Female_C6, Female_C7, Female_C8, Female_C9, Female_C10,
                                   Female_C11, Female_C12, Female_C13, Female_C14, Female_C15, Female_C16, Female_C17, Female_C18,
                                   Female_M1, Female_M2, Female_M3, Female_M4, Female_M5, Female_M6, Female_M7, Female_M8, Female_M9, Female_M10,
                                   Female_M11, Female_M12, Female_M13, Female_M14, Female_M15, Female_M16, Female_M17, Female_M18, Female_M19, Female_M20,
                                   dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                                   dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6))
VlnPlot(alldata_after1, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 1) +
  NoLegend()

### Histogram
metadata2 <- alldata_after1@meta.data
metadata2$sample <- NA
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  metadata2$sample[which(str_detect(metadata2$orig.ident, paste0("^Male_C", i)))] <- paste0("Male_C", i)
}
metadata2$sample[which(str_detect(metadata2$orig.ident, "^Male_C13_1"))] <- "Male_C13"
metadata2$sample[which(str_detect(metadata2$orig.ident, "^Male_C13_2"))] <- "Male_C13"

for(i in 1:17){
  metadata2$sample[which(str_detect(metadata2$orig.ident, paste0("^Male_M", i)))] <- paste0("Male_M", i)
}

for(i in c(1:18)){
  metadata2$sample[which(str_detect(metadata2$orig.ident, paste0("^Female_C", i)))] <- paste0("Female_C", i)
}

for(i in 1:20){
  metadata2$sample[which(str_detect(metadata2$orig.ident, paste0("^Female_M", i)))] <- paste0("Female_M", i)
}


##### Generate metadata using the same pattern-based approach as above. #####


### Visualize the number UMIs/transcripts per cell
metadata2 %>% 
  ggplot(aes(color=group, x=nCount_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
  #geom_vline(xintercept = 10000)

# Visualize the distribution of genes detected per cell via histogram
metadata2 %>% 
  ggplot(aes(color=group, x=nFeature_RNA, fill= group)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
  #geom_vline(xintercept = c(300, 400))


### Remove Mitochondrial gene (MT-)
mito_genes <- rownames(Male_C1)[grep("^MT-", rownames(Male_C1))]
mito_genes
Male_C1 <- Male_C1[!grepl("^MT-", rownames(Male_C1)), ]
Male_C2 <- Male_C2[!grepl("^MT-", rownames(Male_C2)), ]
Male_C3 <- Male_C3[!grepl("^MT-", rownames(Male_C3)), ]
Male_C4 <- Male_C4[!grepl("^MT-", rownames(Male_C4)), ]
Male_C5 <- Male_C5[!grepl("^MT-", rownames(Male_C5)), ]
Male_C6 <- Male_C6[!grepl("^MT-", rownames(Male_C6)), ]
Male_C7 <- Male_C7[!grepl("^MT-", rownames(Male_C7)), ]
Male_C8 <- Male_C8[!grepl("^MT-", rownames(Male_C8)), ]
#Male_C9 <- Male_C9[!grepl("^MT-", rownames(Male_C9)), ]
Male_C10 <- Male_C10[!grepl("^MT-", rownames(Male_C10)), ]
Male_C11 <- Male_C11[!grepl("^MT-", rownames(Male_C11)), ]
Male_C12 <- Male_C12[!grepl("^MT-", rownames(Male_C12)), ]
Male_C13_1 <- Male_C13_1[!grepl("^MT-", rownames(Male_C13_1)), ]
Male_C13_2 <- Male_C13_2[!grepl("^MT-", rownames(Male_C13_2)), ]
Male_C15 <- Male_C15[!grepl("^MT-", rownames(Male_C15)), ]
Male_C16 <- Male_C16[!grepl("^MT-", rownames(Male_C16)), ]
Male_C17 <- Male_C17[!grepl("^MT-", rownames(Male_C17)), ]

Male_M1 <- Male_M1[!grepl("^MT-", rownames(Male_M1)), ]
Male_M2 <- Male_M2[!grepl("^MT-", rownames(Male_M2)), ]
Male_M3 <- Male_M3[!grepl("^MT-", rownames(Male_M3)), ]
Male_M4 <- Male_M4[!grepl("^MT-", rownames(Male_M4)), ]
Male_M5 <- Male_M5[!grepl("^MT-", rownames(Male_M5)), ]
Male_M6 <- Male_M6[!grepl("^MT-", rownames(Male_M6)), ]
Male_M7 <- Male_M7[!grepl("^MT-", rownames(Male_M7)), ]
Male_M8 <- Male_M8[!grepl("^MT-", rownames(Male_M8)), ]
Male_M9 <- Male_M9[!grepl("^MT-", rownames(Male_M9)), ]
Male_M10 <- Male_M10[!grepl("^MT-", rownames(Male_M10)), ]
Male_M11 <- Male_M11[!grepl("^MT-", rownames(Male_M11)), ]
Male_M12 <- Male_M12[!grepl("^MT-", rownames(Male_M12)), ]
Male_M13 <- Male_M13[!grepl("^MT-", rownames(Male_M13)), ]
Male_M14 <- Male_M14[!grepl("^MT-", rownames(Male_M14)), ]
Male_M15 <- Male_M15[!grepl("^MT-", rownames(Male_M15)), ]
Male_M16 <- Male_M16[!grepl("^MT-", rownames(Male_M16)), ]
Male_M17 <- Male_M17[!grepl("^MT-", rownames(Male_M17)), ]

Female_C1 <- Female_C1[!grepl("^MT-", rownames(Female_C1)), ]
Female_C2 <- Female_C2[!grepl("^MT-", rownames(Female_C2)), ]
Female_C3 <- Female_C3[!grepl("^MT-", rownames(Female_C3)), ]
Female_C4 <- Female_C4[!grepl("^MT-", rownames(Female_C4)), ]
Female_C5 <- Female_C5[!grepl("^MT-", rownames(Female_C5)), ]
Female_C6 <- Female_C6[!grepl("^MT-", rownames(Female_C6)), ]
Female_C7 <- Female_C7[!grepl("^MT-", rownames(Female_C7)), ]
Female_C8 <- Female_C8[!grepl("^MT-", rownames(Female_C8)), ]
Female_C9 <- Female_C9[!grepl("^MT-", rownames(Female_C9)), ]
Female_C10 <- Female_C10[!grepl("^MT-", rownames(Female_C10)), ]
Female_C11 <- Female_C11[!grepl("^MT-", rownames(Female_C11)), ]
Female_C12 <- Female_C12[!grepl("^MT-", rownames(Female_C12)), ]
Female_C13 <- Female_C13[!grepl("^MT-", rownames(Female_C13)), ]
Female_C14 <- Female_C14[!grepl("^MT-", rownames(Female_C14)), ]
Female_C15 <- Female_C15[!grepl("^MT-", rownames(Female_C15)), ]
Female_C16 <- Female_C16[!grepl("^MT-", rownames(Female_C16)), ]
Female_C17 <- Female_C17[!grepl("^MT-", rownames(Female_C17)), ]
Female_C18 <- Female_C18[!grepl("^MT-", rownames(Female_C18)), ]

Female_M1 <- Female_M1[!grepl("^MT-", rownames(Female_M1)), ]
Female_M2 <- Female_M2[!grepl("^MT-", rownames(Female_M2)), ]
Female_M3 <- Female_M3[!grepl("^MT-", rownames(Female_M3)), ]
Female_M4 <- Female_M4[!grepl("^MT-", rownames(Female_M4)), ]
Female_M5 <- Female_M5[!grepl("^MT-", rownames(Female_M5)), ]
Female_M6 <- Female_M6[!grepl("^MT-", rownames(Female_M6)), ]
Female_M7 <- Female_M7[!grepl("^MT-", rownames(Female_M7)), ]
Female_M8 <- Female_M8[!grepl("^MT-", rownames(Female_M8)), ]
Female_M9 <- Female_M9[!grepl("^MT-", rownames(Female_M9)), ]
Female_M10 <- Female_M10[!grepl("^MT-", rownames(Female_M10)), ]
Female_M11 <- Female_M11[!grepl("^MT-", rownames(Female_M11)), ]
Female_M12 <- Female_M12[!grepl("^MT-", rownames(Female_M12)), ]
Female_M13 <- Female_M13[!grepl("^MT-", rownames(Female_M13)), ]
Female_M14 <- Female_M14[!grepl("^MT-", rownames(Female_M14)), ]
Female_M15 <- Female_M15[!grepl("^MT-", rownames(Female_M15)), ]
Female_M16 <- Female_M16[!grepl("^MT-", rownames(Female_M16)), ]
Female_M17 <- Female_M17[!grepl("^MT-", rownames(Female_M17)), ]
Female_M18 <- Female_M18[!grepl("^MT-", rownames(Female_M18)), ]
Female_M19 <- Female_M19[!grepl("^MT-", rownames(Female_M19)), ]
Female_M20 <- Female_M20[!grepl("^MT-", rownames(Female_M20)), ]

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
combined.list <- list(Male_C1, Male_C2, Male_C3, Male_C4, Male_C5, Male_C6, Male_C7, Male_C8, Male_C10,
                      Male_C11, Male_C12, Male_C13_1, Male_C13_2, Male_C15, Male_C16, Male_C17,
                      Male_M1, Male_M2, Male_M3, Male_M4, Male_M5, Male_M6, Male_M7, Male_M8, Male_M9, Male_M10,
                      Male_M11, Male_M12, Male_M13, Male_M14, Male_M15, Male_M16, Male_M17,
                      Female_C1, Female_C2, Female_C3, Female_C4, Female_C5, Female_C6, Female_C7, Female_C8, Female_C9, Female_C10,
                      Female_C11, Female_C12, Female_C13, Female_C14, Female_C15, Female_C16, Female_C17, Female_C18,
                      Female_M1, Female_M2, Female_M3, Female_M4, Female_M5, Female_M6, Female_M7, Female_M8, Female_M9, Female_M10,
                      Female_M11, Female_M12, Female_M13, Female_M14, Female_M15, Female_M16, Female_M17, Female_M18, Female_M19, Female_M20,
                      dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
                      dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6)
names(combined.list) <- c("Male_C1", "Male_C2", "Male_C3", "Male_C4", "Male_C5", "Male_C6", "Male_C7", "Male_C8", "Male_C10",
                          "Male_C11", "Male_C12", "Male_C13_1", "Male_C13_2", "Male_C15", "Male_C16", "Male_C17",
                          "Male_M1", "Male_M2", "Male_M3", "Male_M4", "Male_M5", "Male_M6", "Male_M7", "Male_M8", "Male_M9", "Male_M10",
                          "Male_M11", "Male_M12", "Male_M13", "Male_M14", "Male_M15", "Male_M16", "Male_M17",
                          "Female_C1", "Female_C2", "Female_C3", "Female_C4", "Female_C5", "Female_C6", "Female_C7", "Female_C8", "Female_C9", "Female_C10",
                          "Female_C11", "Female_C12", "Female_C13", "Female_C14", "Female_C15", "Female_C16", "Female_C17", "Female_C18",
                          "Female_M1", "Female_M2", "Female_M3", "Female_M4", "Female_M5", "Female_M6", "Female_M7", "Female_M8", "Female_M9", "Female_M10",
                          "Female_M11", "Female_M12", "Female_M13", "Female_M14", "Female_M15", "Female_M16", "Female_M17", "Female_M18", "Female_M19", "Female_M20",
                          "dlPFC_C1", "dlPFC_C2", "dlPFC_C3", "dlPFC_C4", "dlPFC_C5", "dlPFC_C6",
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
combined.list <- readRDS(file="Seurat.before.DF.RDS")


### DoubletFinder ## Nagy doublet rate 2.3%
sweep.res <- paramSweep(combined.list$Female_M20, PCs=1:10)
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


homotypic.prop <- modelHomotypic(combined.list$Female_M20@meta.data$RNA_snn_res.0.8)
nExp <- round(ncol(combined.list$Female_M20) * 0.023)  # expect 2.3% doublets (3,000 cells)
nExp.adj <- round(nExp*(1 - homotypic.prop))
combined.list$Female_M20 <- doubletFinder(combined.list$Female_M20, pN = 0.25, pK = 0.18, nExp = nExp,                                         PCs = 1:10, reuse.pANN=F)
head(combined.list$Female_M20@meta.data)[9]
combined.list$Female_M20 <- doubletFinder(combined.list$Female_M20, pN = 0.25, pK = 0.18, nExp = nExp.adj,
                                           PCs = 1:10, reuse.pANN="pANN_0.25_0.18_28")

###
### DoubletFinder ## My doublet rate 8%
sweep.res <- paramSweep(combined.list$dlPFC_M6, PCs=1:10)
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
combined.list$dlPFC_M6 <- doubletFinder(combined.list$dlPFC_M6, pN = 0.25, pK = 0.08, nExp = nExp,                                         PCs = 1:10, reuse.pANN=F)
head(combined.list$dlPFC_M6@meta.data)[9]
combined.list$dlPFC_M6 <- doubletFinder(combined.list$dlPFC_M6, pN = 0.25, pK = 0.08, nExp = nExp.adj,
                                           PCs = 1:10, reuse.pANN="pANN_0.25_0.08_333")

# save
saveRDS(object=combined.list, file="Seurat.after.DF_Homotypic.RDS")

# Change column name
head(combined.list$Male_C1@meta.data)
colnames(combined.list$Male_C1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C7@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C8@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
#colnames(combined.list$Male_C9@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C10@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C11@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C12@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C13_1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C13_2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C15@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C16@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_C17@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$Male_M1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M7@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M8@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M9@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M10@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M11@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M12@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M13@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M14@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M15@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M16@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Male_M17@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$Female_C1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C7@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C8@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C9@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C10@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C11@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C12@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C13@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C14@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C15@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C16@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C17@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_C18@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

colnames(combined.list$Female_M1@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M2@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M3@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M4@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M5@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M6@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M7@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M8@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M9@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M10@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M11@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M12@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M13@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M14@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M15@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M16@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M17@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M18@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M19@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")
colnames(combined.list$Female_M20@meta.data)[9:11] <- c("DF_score", "DF_class_low", "DF_class_high")

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
combined.list$Male_C1 <- subset(combined.list$Male_C1, subset = DF_class_high == "Singlet")
combined.list$Male_C2 <- subset(combined.list$Male_C2, subset = DF_class_high == "Singlet")
combined.list$Male_C3 <- subset(combined.list$Male_C3, subset = DF_class_high == "Singlet")
combined.list$Male_C4 <- subset(combined.list$Male_C4, subset = DF_class_high == "Singlet")
combined.list$Male_C5 <- subset(combined.list$Male_C5, subset = DF_class_high == "Singlet")
combined.list$Male_C6 <- subset(combined.list$Male_C6, subset = DF_class_high == "Singlet")
combined.list$Male_C7 <- subset(combined.list$Male_C7, subset = DF_class_high == "Singlet")
combined.list$Male_C8 <- subset(combined.list$Male_C8, subset = DF_class_high == "Singlet")
#combined.list$Male_C9 <- subset(combined.list$Male_C9, subset = DF_class_high == "Singlet")
combined.list$Male_C10 <- subset(combined.list$Male_C10, subset = DF_class_high == "Singlet")
combined.list$Male_C11 <- subset(combined.list$Male_C11, subset = DF_class_high == "Singlet")
combined.list$Male_C12 <- subset(combined.list$Male_C12, subset = DF_class_high == "Singlet")
combined.list$Male_C13_1 <- subset(combined.list$Male_C13_1, subset = DF_class_high == "Singlet")
combined.list$Male_C13_2 <- subset(combined.list$Male_C13_2, subset = DF_class_high == "Singlet")
combined.list$Male_C15 <- subset(combined.list$Male_C15, subset = DF_class_high == "Singlet")
combined.list$Male_C16 <- subset(combined.list$Male_C16, subset = DF_class_high == "Singlet")
combined.list$Male_C17 <- subset(combined.list$Male_C17, subset = DF_class_high == "Singlet")

combined.list$Male_M1 <- subset(combined.list$Male_M1, subset = DF_class_high == "Singlet")
combined.list$Male_M2 <- subset(combined.list$Male_M2, subset = DF_class_high == "Singlet")
combined.list$Male_M3 <- subset(combined.list$Male_M3, subset = DF_class_high == "Singlet")
combined.list$Male_M4 <- subset(combined.list$Male_M4, subset = DF_class_high == "Singlet")
combined.list$Male_M5 <- subset(combined.list$Male_M5, subset = DF_class_high == "Singlet")
combined.list$Male_M6 <- subset(combined.list$Male_M6, subset = DF_class_high == "Singlet")
combined.list$Male_M7 <- subset(combined.list$Male_M7, subset = DF_class_high == "Singlet")
combined.list$Male_M8 <- subset(combined.list$Male_M8, subset = DF_class_high == "Singlet")
combined.list$Male_M9 <- subset(combined.list$Male_M9, subset = DF_class_high == "Singlet")
combined.list$Male_M10 <- subset(combined.list$Male_M10, subset = DF_class_high == "Singlet")
combined.list$Male_M11 <- subset(combined.list$Male_M11, subset = DF_class_high == "Singlet")
combined.list$Male_M12 <- subset(combined.list$Male_M12, subset = DF_class_high == "Singlet")
combined.list$Male_M13 <- subset(combined.list$Male_M13, subset = DF_class_high == "Singlet")
combined.list$Male_M14 <- subset(combined.list$Male_M14, subset = DF_class_high == "Singlet")
combined.list$Male_M15 <- subset(combined.list$Male_M15, subset = DF_class_high == "Singlet")
combined.list$Male_M16 <- subset(combined.list$Male_M16, subset = DF_class_high == "Singlet")
combined.list$Male_M17 <- subset(combined.list$Male_M17, subset = DF_class_high == "Singlet")

combined.list$Female_C1 <- subset(combined.list$Female_C1, subset = DF_class_high == "Singlet")
combined.list$Female_C2 <- subset(combined.list$Female_C2, subset = DF_class_high == "Singlet")
combined.list$Female_C3 <- subset(combined.list$Female_C3, subset = DF_class_high == "Singlet")
combined.list$Female_C4 <- subset(combined.list$Female_C4, subset = DF_class_high == "Singlet")
combined.list$Female_C5 <- subset(combined.list$Female_C5, subset = DF_class_high == "Singlet")
combined.list$Female_C6 <- subset(combined.list$Female_C6, subset = DF_class_high == "Singlet")
combined.list$Female_C7 <- subset(combined.list$Female_C7, subset = DF_class_high == "Singlet")
combined.list$Female_C8 <- subset(combined.list$Female_C8, subset = DF_class_high == "Singlet")
combined.list$Female_C9 <- subset(combined.list$Female_C9, subset = DF_class_high == "Singlet")
combined.list$Female_C10 <- subset(combined.list$Female_C10, subset = DF_class_high == "Singlet")
combined.list$Female_C11 <- subset(combined.list$Female_C11, subset = DF_class_high == "Singlet")
combined.list$Female_C12 <- subset(combined.list$Female_C12, subset = DF_class_high == "Singlet")
combined.list$Female_C13 <- subset(combined.list$Female_C13, subset = DF_class_high == "Singlet")
combined.list$Female_C14 <- subset(combined.list$Female_C14, subset = DF_class_high == "Singlet")
combined.list$Female_C15 <- subset(combined.list$Female_C15, subset = DF_class_high == "Singlet")
combined.list$Female_C16 <- subset(combined.list$Female_C16, subset = DF_class_high == "Singlet")
combined.list$Female_C17 <- subset(combined.list$Female_C17, subset = DF_class_high == "Singlet")
combined.list$Female_C18 <- subset(combined.list$Female_C18, subset = DF_class_high == "Singlet")

combined.list$Female_M1 <- subset(combined.list$Female_M1, subset = DF_class_high == "Singlet")
combined.list$Female_M2 <- subset(combined.list$Female_M2, subset = DF_class_high == "Singlet")
combined.list$Female_M3 <- subset(combined.list$Female_M3, subset = DF_class_high == "Singlet")
combined.list$Female_M4 <- subset(combined.list$Female_M4, subset = DF_class_high == "Singlet")
combined.list$Female_M5 <- subset(combined.list$Female_M5, subset = DF_class_high == "Singlet")
combined.list$Female_M6 <- subset(combined.list$Female_M6, subset = DF_class_high == "Singlet")
combined.list$Female_M7 <- subset(combined.list$Female_M7, subset = DF_class_high == "Singlet")
combined.list$Female_M8 <- subset(combined.list$Female_M8, subset = DF_class_high == "Singlet")
combined.list$Female_M9 <- subset(combined.list$Female_M9, subset = DF_class_high == "Singlet")
combined.list$Female_M10 <- subset(combined.list$Female_M10, subset = DF_class_high == "Singlet")
combined.list$Female_M11 <- subset(combined.list$Female_M11, subset = DF_class_high == "Singlet")
combined.list$Female_M12 <- subset(combined.list$Female_M12, subset = DF_class_high == "Singlet")
combined.list$Female_M13 <- subset(combined.list$Female_M13, subset = DF_class_high == "Singlet")
combined.list$Female_M14 <- subset(combined.list$Female_M14, subset = DF_class_high == "Singlet")
combined.list$Female_M15 <- subset(combined.list$Female_M15, subset = DF_class_high == "Singlet")
combined.list$Female_M16 <- subset(combined.list$Female_M16, subset = DF_class_high == "Singlet")
combined.list$Female_M17 <- subset(combined.list$Female_M17, subset = DF_class_high == "Singlet")
combined.list$Female_M18 <- subset(combined.list$Female_M18, subset = DF_class_high == "Singlet")
combined.list$Female_M19 <- subset(combined.list$Female_M19, subset = DF_class_high == "Singlet")
combined.list$Female_M20 <- subset(combined.list$Female_M20, subset = DF_class_high == "Singlet")

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
combined.list <- readRDS(file="Seurat.remove.DF_Homotypic.RDS")


# Cell counting after QC
quantile(combined.list$Female_M20$nFeature_RNA,probs = c(0.5))
quantile(combined.list$Female_M20$nCount_RNA,probs = c(0.5))

alldata_after <- merge(combined.list$Male_C1, c(combined.list$Male_C2, combined.list$Male_C3, combined.list$Male_C4, combined.list$Male_C5, combined.list$Male_C6, combined.list$Male_C7, combined.list$Male_C8, combined.list$Male_C10,
                                  combined.list$Male_C11, combined.list$Male_C12, combined.list$Male_C13_1, combined.list$Male_C13_2, combined.list$Male_C15, combined.list$Male_C16, combined.list$Male_C17,
                                  combined.list$Male_M1, combined.list$Male_M2, combined.list$Male_M3, combined.list$Male_M4, combined.list$Male_M5, combined.list$Male_M6, combined.list$Male_M7, combined.list$Male_M8, combined.list$Male_M9, combined.list$Male_M10,
                                  combined.list$Male_M11, combined.list$Male_M12, combined.list$Male_M13, combined.list$Male_M14, combined.list$Male_M15, combined.list$Male_M16, combined.list$Male_M17,
                                  combined.list$Female_C1, combined.list$Female_C2, combined.list$Female_C3, combined.list$Female_C4, combined.list$Female_C5, combined.list$Female_C6, combined.list$Female_C7, combined.list$Female_C8, combined.list$Female_C9, combined.list$Female_C10,
                                  combined.list$Female_C11, combined.list$Female_C12, combined.list$Female_C13, combined.list$Female_C14, combined.list$Female_C15, combined.list$Female_C16, combined.list$Female_C17, combined.list$Female_C18,
                                  combined.list$Female_M1, combined.list$Female_M2, combined.list$Female_M3, combined.list$Female_M4, combined.list$Female_M5, combined.list$Female_M6, combined.list$Female_M7, combined.list$Female_M8, combined.list$Female_M9, combined.list$Female_M10,
                                  combined.list$Female_M11, combined.list$Female_M12, combined.list$Female_M13, combined.list$Female_M14, combined.list$Female_M15, combined.list$Female_M16, combined.list$Female_M17, combined.list$Female_M18, combined.list$Female_M19, combined.list$Female_M20,
                                  combined.list$dlPFC_C1, combined.list$dlPFC_C2, combined.list$dlPFC_C3, combined.list$dlPFC_C4, combined.list$dlPFC_C5, combined.list$dlPFC_C6,
                                  combined.list$dlPFC_M1, combined.list$dlPFC_M2, combined.list$dlPFC_M3, combined.list$dlPFC_M4, combined.list$dlPFC_M5, combined.list$dlPFC_M6))

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hemo")
VlnPlot(alldata_after, group.by = "orig.ident", features = feats3, pt.size = 0.1, ncol = 1) +
  NoLegend()

### Histogram
metadata3 <- alldata_after@meta.data
metadata3$sample <- NA
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  metadata3$sample[which(str_detect(metadata3$orig.ident, paste0("^Male_C", i)))] <- paste0("Male_C", i)
}
metadata3$sample[which(str_detect(metadata3$orig.ident, "^Male_C13_1"))] <- "Male_C13"
metadata3$sample[which(str_detect(metadata3$orig.ident, "^Male_C13_2"))] <- "Male_C13"

for(i in 1:17){
  metadata3$sample[which(str_detect(metadata3$orig.ident, paste0("^Male_M", i)))] <- paste0("Male_M", i)
}

for(i in c(1:18)){
  metadata3$sample[which(str_detect(metadata3$orig.ident, paste0("^Female_C", i)))] <- paste0("Female_C", i)
}

for(i in 1:20){
  metadata3$sample[which(str_detect(metadata3$orig.ident, paste0("^Female_M", i)))] <- paste0("Female_M", i)
}

##### Generate metadata using the same pattern-based approach as above. #####


### Visualize the number UMIs/transcripts per cell
metadata3 %>% 
  ggplot(aes(color=case, x=nCount_RNA, fill= case)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
#geom_vline(xintercept = 10000)

# Visualize the distribution of genes detected per cell via histogram
metadata3 %>% 
  ggplot(aes(color=case, x=nFeature_RNA, fill= case)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()
#geom_vline(xintercept = c(300, 400))


# save (alldata before, after)
saveRDS(object=alldata_after, file="Alldata_after_QC&DF.RDS")

rm(Male_C1, Male_C2, Male_C3, Male_C4, Male_C5, Male_C6, Male_C7, Male_C8, Male_C10,
   Male_C11, Male_C12, Male_C13_1, Male_C13_2,Male_C15, Male_C16, Male_C17,
   Male_M1, Male_M2, Male_M3, Male_M4, Male_M5, Male_M6, Male_M7, Male_M8, Male_M9, Male_M10,
   Male_M11, Male_M12, Male_M13, Male_M14, Male_M15, Male_M16, Male_M17)
rm(Female_C1, Female_C2, Female_C3, Female_C4, Female_C5, Female_C6, Female_C7, Female_C8, Female_C9, Female_C10,
   Female_C11, Female_C12, Female_C13, Female_C14,Female_C15, Female_C16, Female_C17, Female_C18,
   Female_M1, Female_M2, Female_M3, Female_M4, Female_M5, Female_M6, Female_M7, Female_M8, Female_M9, Female_M10,
   Female_M11, Female_M12, Female_M13, Female_M14, Female_M15, Female_M16, Female_M17, Female_M18, Female_M19, Female_M20)
rm(dlPFC_C1, dlPFC_C2, dlPFC_C3, dlPFC_C4, dlPFC_C5, dlPFC_C6,
   dlPFC_M1, dlPFC_M2, dlPFC_M3, dlPFC_M4, dlPFC_M5, dlPFC_M6)
rm(alldata_before, alldata_after1)


### Integration v5
combined.list <- readRDS(file="Seurat.remove.DF_Homotypic.RDS")
Integ.obj <- merge(combined.list$Male_C1, c(combined.list$Male_C2, combined.list$Male_C3, combined.list$Male_C4, combined.list$Male_C5, combined.list$Male_C6, combined.list$Male_C7, combined.list$Male_C8, combined.list$Male_C10,
                                            combined.list$Male_C11, combined.list$Male_C12, combined.list$Male_C13_1, combined.list$Male_C13_2, combined.list$Male_C15, combined.list$Male_C16, combined.list$Male_C17,
                                            combined.list$Male_M1, combined.list$Male_M2, combined.list$Male_M3, combined.list$Male_M4, combined.list$Male_M5, combined.list$Male_M6, combined.list$Male_M7, combined.list$Male_M8, combined.list$Male_M9, combined.list$Male_M10,
                                            combined.list$Male_M11, combined.list$Male_M12, combined.list$Male_M13, combined.list$Male_M14, combined.list$Male_M15, combined.list$Male_M16, combined.list$Male_M17,
                                            combined.list$Female_C1, combined.list$Female_C2, combined.list$Female_C3, combined.list$Female_C4, combined.list$Female_C5, combined.list$Female_C6, combined.list$Female_C7, combined.list$Female_C8, combined.list$Female_C9, combined.list$Female_C10,
                                            combined.list$Female_C11, combined.list$Female_C12, combined.list$Female_C13, combined.list$Female_C14, combined.list$Female_C15, combined.list$Female_C16, combined.list$Female_C17, combined.list$Female_C18,
                                            combined.list$Female_M1, combined.list$Female_M2, combined.list$Female_M3, combined.list$Female_M4, combined.list$Female_M5, combined.list$Female_M6, combined.list$Female_M7, combined.list$Female_M8, combined.list$Female_M9, combined.list$Female_M10,
                                            combined.list$Female_M11, combined.list$Female_M12, combined.list$Female_M13, combined.list$Female_M14, combined.list$Female_M15, combined.list$Female_M16, combined.list$Female_M17, combined.list$Female_M18, combined.list$Female_M19, combined.list$Female_M20,
                                            combined.list$dlPFC_C1, combined.list$dlPFC_C2, combined.list$dlPFC_C3, combined.list$dlPFC_C4, combined.list$dlPFC_C5, combined.list$dlPFC_C6,
                                            combined.list$dlPFC_M1, combined.list$dlPFC_M2, combined.list$dlPFC_M3, combined.list$dlPFC_M4, combined.list$dlPFC_M5, combined.list$dlPFC_M6))

Integ.obj[["RNA"]] <- JoinLayers(Integ.obj[["RNA"]])
Integ.obj

### Metadata
metadata <- Integ.obj@meta.data
metadata$sample <- NA
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Male_C", i)))] <- paste0("Male_C", i)
}
metadata$sample[which(str_detect(metadata$orig.ident, "^Male_C13_1"))] <- "Male_C13"
metadata$sample[which(str_detect(metadata$orig.ident, "^Male_C13_2"))] <- "Male_C13"

for(i in 1:17){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Male_M", i)))] <- paste0("Male_M", i)
}

for(i in c(1:18)){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Female_C", i)))] <- paste0("Female_C", i)
}

for(i in 1:20){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Female_M", i)))] <- paste0("Female_M", i)
}


##### Generate metadata using the same pattern-based approach as above. #####

View(Integ.obj@meta.data)

saveRDS(object=Integ.obj, file="Before_Integration_split.RDS")
Integ.obj <- readRDS("Before_Integration_split.RDS")


# split the RNA measurements batch / sample / dataset
Integ.obj[["RNA"]] <- split(Integ.obj[["RNA"]], f = Integ.obj$sample)
Integ.obj


###############################
### read object before split and integration
set.seed(1234)
library(future)
options(future.globals.maxSize = 8000 * 1024^2)
Integ.obj <- readRDS("Before_Integration_split.RDS")
Integ.obj[["RNA"]] <- JoinLayers(Integ.obj[["RNA"]])
Integ.obj

# split the RNA measurements sample
Integ.obj[["RNA"]] <- split(Integ.obj[["RNA"]], f = Integ.obj$sample)
Integ.obj

# run standard anlaysis workflow
Integ.obj <- NormalizeData(Integ.obj)
Integ.obj <- FindVariableFeatures(Integ.obj)
Integ.obj <- ScaleData(Integ.obj, vars.to.regress = "nCount_RNA", verbose = FALSE)
Integ.obj <- RunPCA(Integ.obj, npcs = 100, verbose = FALSE)
Integ.obj <- FindNeighbors(Integ.obj, dims = 1:30, reduction = "pca")
Integ.obj <- FindClusters(Integ.obj, resolution = 2, cluster.name = "unintegrated_clusters")
Integ.obj <- RunUMAP(Integ.obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(Integ.obj, reduction = "umap.unintegrated", group.by = c("batch", "dataset", "sample"))
DimPlot(Integ.obj, reduction = "umap.unintegrated", group.by = c("seurat_clusters"))

saveRDS(object=Integ.obj, file="Before_Integration_Harmony.RDS")
Integ.obj <- readRDS("Before_Integration_Harmony.RDS")

### Integration Harmony
Integ.obj <- IntegrateLayers(object = Integ.obj, method = HarmonyIntegration, orig.reduction = "pca",
                             new.reduction = "harmony", verbose = FALSE)
#epsilon.cluster = -Inf, epsilon.harmony = -Inf

saveRDS(object=Integ.obj, file="Integration_Harmony.RDS")
Integ.obj <- readRDS("Integration_Harmony.RDS")

### PCA selection (Elbow, JackStraw)
ElbowPlot(object=Integ.obj, ndims=100)

Integ.obj <- JackStraw(object=Integ.obj, num.replicate=100, dims=100)
Integ.obj <- ScoreJackStraw(object=Integ.obj, dims=1:100)

saveRDS(object=Integ.obj, file="After_Jack.RDS")
Integ.obj <- readRDS("After_Jack.RDS")

###
# re-join layers after integration
Integ.obj[["RNA"]] <- JoinLayers(Integ.obj[["RNA"]])

JackStrawPlot(object=Integ.obj, dims=1:100)
JackStrawPlot(object=Integ.obj, dims=1:72) + NoLegend() # PC72

Integ.obj <- FindNeighbors(Integ.obj, reduction = "harmony", dims = 1:72)
Integ.obj <- FindClusters(Integ.obj, resolution = 1.3)
Integ.obj <- RunUMAP(Integ.obj, dims = 1:72, reduction = "harmony", reduction.name = "umap.harmony")

DimPlot(Integ.obj, reduction = "umap.harmony", label=T, repel=F, raster=F)

# Visualization
p1 <- DimPlot(Integ.obj, reduction = "umap.unintegrated", group.by = c("batch", "dataset", "sample")) + NoLegend()
p2 <- DimPlot(Integ.obj, reduction = "umap.harmony", group.by = c("batch", "dataset", "sample")) + NoLegend()

p1 / p2

DefaultAssay(Integ.obj)
Marker <- c("RBFOX3", "MAP2", "GFAP", "SLC1A2", "CX3CR1", "CD74", "MBP", "PDGFRA", "CLDN5")
FeaturePlot(object=Integ.obj, reduction="umap.harmony", features=Marker, label=F, min.cutoff = "q9", ncol=3)
FeaturePlot(object=Integ.obj, reduction="umap.harmony", features=c("THY1", "NRGN", "BEX1"), label=F, min.cutoff = "q9", ncol=3)

saveRDS(object=Integ.obj, file="Integration_Harmony_res1.3.RDS")
Integ.obj <- readRDS("Integration_Harmony_res1.3.RDS")

View(Integ.obj@meta.data)

### Clustering distribution, Gender-specific gene
DefaultAssay(Integ.obj) <- "RNA"
VlnPlot(object=Integ.obj, features=c("UTY", "USP9Y","XIST", "TSIX"), group.by="sample", ncol=1, pt.size=0.1)
VlnPlot(object=Integ.obj, features=c("UTY", "USP9Y","XIST", "TSIX"), group.by="sample", ncol=1, pt.size=0)


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
Integ.obj <- CellCycleScoring(object=Integ.obj, s.features=s.genes, g2m.features=g2m.genes)
View(Integ.obj@meta.data)

DimPlot(Integ.obj,
        reduction = "pca",
        group.by= "Phase",
        split.by="Phase")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(Integ.obj, reduction="umap.harmony", features=metrics, min.cutoff = "q9", ncol=3)

# save
saveRDS(object=Integ.obj, file="Integration_Harmony_Final_res1.3.RDS")

### tSNE (long)
Integ.obj <- RunTSNE(Integ.obj, dims=1:72)
DimPlot(Integ.obj, reduction="tsne", label=T, repel=F, raster=F)

### Save object
saveRDS(object=Integ.obj, file="Seurat.UMAP_tSNE.RDS")


### Metadata
metadata <- Integ.obj@meta.data
metadata$sample <- NA
for(i in c(1:8, 10, 11, 12, 15, 16, 17)){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Male_C", i)))] <- paste0("Male_C", i)
}
metadata$sample[which(str_detect(metadata$orig.ident, "^Male_C13_1"))] <- "Male_C13"
metadata$sample[which(str_detect(metadata$orig.ident, "^Male_C13_2"))] <- "Male_C13"

for(i in 1:17){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Male_M", i)))] <- paste0("Male_M", i)
}

for(i in c(1:18)){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Female_C", i)))] <- paste0("Female_C", i)
}

for(i in 1:20){
  metadata$sample[which(str_detect(metadata$orig.ident, paste0("^Female_M", i)))] <- paste0("Female_M", i)
}


##### Generate metadata using the same pattern-based approach as above. #####

### Save object
saveRDS(object=Integ.obj, file="Seurat.UMAP_tSNE.RDS")
Integ.obj <- readRDS("Seurat.UMAP_tSNE.RDS")

count_table <- table(Integ.obj@meta.data$RNA_snn_res.1.3, Integ.obj@meta.data$sample)
count_table
write.csv(count_table, "Cell_proportion/Cell_count_table.csv")



### Cluster Dendrogram
Integ.obj <- BuildClusterTree(Integ.obj)
myPhyTree <- Tool(object=Integ.obj, slot = "BuildClusterTree")
tiff(filename = "Dendrogram.tif", width = 5, height = 11, units = "in", res=300)
pdf("Dendrogram.pdf", width = 5, height = 13)
ggtree(myPhyTree)+geom_tiplab()+theme_tree()+xlim(NA,200)
#ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1) + 
#  geom_tiplab(aes(label=myPhyTree$tip.label))
#PlotClusterTree(Integ.obj, direction = "rightwards")
dev.off()


levels(Integ.obj)
levels(Integ.obj) <- c("15", "26", "43", "3", "32", "41", "21", "36", "42", "2", "27",
                       "9", "40", "23", "25", "22", "39", "13", "0", "4",
                       "11", "29", "6", "38", "17", "34", "8", "19", "31", "28", "35", "10", "20",
                       "45", "30", "46", "37", "44", "12", "1", "5", "33", "7", "16", "14", "18", "24")

levels(Integ.obj)

### Remove cluster (after tSNE)
Integ.obj <- subset(Integ.obj, idents=c(32,41,45,46), invert=T)
Integ.obj <- subset(Integ.obj, idents=c("PQ4"), invert=T)
DimPlot(Integ.obj, reduction="umap.harmony", label=T, repel=T, raster=F)
DimPlot(Integ.obj, reduction="tsne", label=T, repel=T, raster=F)

saveRDS(object=Integ.obj, file="Seurat.UMAP_rm.RDS")
Integ.obj <- readRDS("Seurat.UMAP_rm.RDS")

metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hemo", "S.Score", "G2M.Score")
FeaturePlot(Integ.obj, reduction="umap.harmony", features=metrics, min.cutoff = "q9", ncol=3)

levels(Integ.obj)

### DEG & Heatmap (too long)
seurat.markers <- FindAllMarkers(object=Integ.obj, only.pos=T, min.pct=0.1, logfc.threshold=0.25) # 0.1 / 0.25
head(seurat.markers)

seurat.top10.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
seurat.top20.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

DoHeatmap(object=Integ.obj, features=seurat.top10.genes$gene, size=4) + NoLegend()

### Save
write.csv(seurat.markers, "Seurat.FindAllMarkers.csv")
write.csv(seurat.top10.genes, "Seurat.top10markers.csv")
write.csv(seurat.top20.genes, "Seurat.top20markers.csv")
saveRDS(object=seurat.markers, file="Seurat.FindAllMarker.RDS")

View(Integ.obj@meta.data)

### UMAP
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="group")
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="sample")
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="case")
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="dataset")
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="batch")

DimPlot(object=Integ.obj, reduction="umap.harmony", split.by="case", ncol=2)
DimPlot(object=Integ.obj, reduction="umap.harmony", split.by="dataset", ncol=3)

DimPlot(object=Integ.obj, reduction="umap.harmony", split.by="integrated_snn_res.1.3", ncol=6)

DimPlot(object=Integ.obj, reduction="tsne", group.by="group")
DimPlot(object=Integ.obj, reduction="tsne", group.by="sample")

DimPlot(object=Integ.obj, reduction="tsne", split.by="group", ncol=2)
DimPlot(object=Integ.obj, reduction="tsne", split.by="sample", ncol=6)

levels(Integ.obj)


### Convert Cluster ID
levels(Integ.obj) #32,41,45,46
levels(Integ.obj) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                       "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                       "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                       "30", "31", "33", "34", "35", "36", "37", "38", "39",
                       "40", "42", "43", "44")

new.cluster.ids <- c("Oligo1", "ExN1", "PQ1", "Astro1", "Oligo2", "ExN2", "OPC1", "ExN3", "InN1", "Micro1",
                     "InN2", "InN3", "ExN4", "Oligo3", "ExN5", "Astro2", "ExN6", "ExN7", "ExN8", "InN4",
                     "InN5", "PQ2", "Oligo4", "Endo1", "ExN9", "Endo2", "Astro3", "PQ3", "InN6", "InN7",
                     "Neuron", "InN8", "ExN10", "ExN11", "InN9", "PQ4", "ExN12", "OPC2", "Micro2",
                     "Micro3", "PQ5", "Astro4", "ExN13")
names(new.cluster.ids) <- levels(Integ.obj)
Integ.obj <- RenameIdents(Integ.obj, new.cluster.ids)
Integ.obj@meta.data$"Newcluster" <- as.factor(Integ.obj@active.ident)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=4, repel=F)



### Cluster Dendrogram
t1 <- myPhyTree$tip.label
Integ.obj <- BuildClusterTree(Integ.obj)
myPhyTree <- Tool(object=Integ.obj, slot = "BuildClusterTree")
#tiff(filename = "Dendrogram5.tif", width = 5, height = 11, units = "in", res=300)
pdf("Dendrogram.pdf", width = 5, height = 13)
ggtree(myPhyTree)+
  geom_tiplab(color = c("#967341", "#6c43bf", "#a3a3a3", "#d9252f", "#967341", "#6c43bf", "#FFB500", "#6c43bf", "#e897c7", "#3aa312",
                        "#e897c7", "#e897c7", "#6c43bf", "#967341", "#6c43bf", "#d9252f", "#6c43bf", "#6c43bf", "#6c43bf", "#e897c7",
                        "#e897c7", "#a3a3a3", "#967341", "#2e75d1", "#6c43bf", "#2e75d1", "#d9252f", "#a3a3a3", "#e897c7", "#e897c7",
                        "#6c43bf", "#e897c7", "#6c43bf", "#6c43bf", "#e897c7", "#6c43bf", "#FFB500", "#3aa312", "#3aa312",
                        "#a3a3a3", "#d9252f", "#6c43bf")) +
  theme_tree()+xlim(NA,200)
#ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1) + 
#  geom_tiplab(aes(label=myPhyTree$tip.label))
#PlotClusterTree(Integ.obj, direction = "rightwards")
dev.off()

pdf("Dendrogram.pdf", width = 5, height = 12)
ggtree(myPhyTree)+
  geom_tiplab(color = c("#d9252f","#d9252f","#d9252f","#d9252f","#a3a3a3","#a3a3a3","#a3a3a3","#a3a3a3",
                        "#3aa312","#3aa312","#2e75d1","#2e75d1","#967341","#3aa312","#967341","#967341","#967341",
                        "#e897c7","#e897c7","#FFB500","#FFB500","#6c43bf","#6c43bf",
                        "#e897c7","#e897c7","#e897c7","#e897c7","#e897c7","#e897c7","#e897c7",
                        "#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf","#6c43bf"),
              size=5) +
  theme_tree()+xlim(NA,350)
#ape::plot.phylo(x=myPhyTree, direction="rightwards", font=1) + 
#  geom_tiplab(aes(label=myPhyTree$tip.label))
#PlotClusterTree(Integ.obj, direction = "rightwards")
dev.off()



levels(Integ.obj) <- c("Astro2", "Astro1", "Astro3", "Astro4", "PQ2", "PQ4", "PQ5", "PQ1", "PQ3",
                       "Micro1", "Micro3", "Endo1", "Endo2", "Oligo4", "Micro2", "Oligo3", "Oligo1", "Oligo2",
                       "InN3", "InN7", "OPC1", "OPC2", "ExN7", "ExN11",
                       "InN1", "InN4", "InN8", "InN6", "InN9", "InN2", "InN5",
                       "Neuron", "ExN12", "ExN13", "ExN4", "ExN1", "ExN2", "ExN10", "ExN3",
                       "ExN6", "ExN5", "ExN8", "ExN9")

levels(Integ.obj)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=4, repel=T)

new.cluster.ids <- c("Astro2", "Astro1", "Astro3", "Astro4", "PQ2", "PQ4", "PQ1", "PQ3",
                     "Micro", "Macro", "Endo1", "Endo2", "Oligo4", "Micro_Phago", "Oligo3", "Oligo1", "Oligo2",
                     "InN3_MGE", "InN7_MGE_PVALB", "OPC1", "OPC2", "ExN7_L5-6", "ExN12_L5-6",
                     "InN1_MGE_SST", "InN4_CGE_LAMP5", "InN8_MGE_NOS1", "InN6_CGE_NDNF", "InN9_CGE", "InN2_CGE_VIP", "InN5_CGE_CALB2",
                     "ExN10_L5-6", "ExN13_L6", "ExN14_L4-5", "ExN4_L2-6", "ExN1_L2-4", "ExN2_L2-4", "ExN11_L4-5", "ExN3_L4-5",
                     "ExN6_L4-5", "ExN5_L2-4", "ExN8_L4-5", "ExN9_L4-5")
names(new.cluster.ids) <- levels(Integ.obj)
Integ.obj <- RenameIdents(Integ.obj, new.cluster.ids)
Integ.obj@meta.data$"Newcluster" <- as.factor(Integ.obj@active.ident)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=4, repel=T)

pdf("UMAP_Final.pdf", width = 7.126, height = 7)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=4, repel=T) + NoLegend()
dev.off()

pdf("UMAP_Case.pdf", width = 7.864, height = 7)
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="case")
dev.off()

pdf("UMAP_Dataset.pdf", width = 7.875, height = 7)
DimPlot(object=Integ.obj, reduction="umap.harmony", group.by="dataset")
dev.off()


### Cell proportion
count_table <- table(Integ.obj@meta.data$Newcluster, Integ.obj@meta.data$sample)
count_table
write.csv(count_table, "Cell_proportion/Cell_count_table_label.csv")

dittoBarPlot(
  object = Integ.obj,
  var = "Newcluster",
  group.by = "dataset")

dittoBarPlot(
  object = Integ.obj,
  var = "sample",
  group.by = "Newcluster")

### save
saveRDS(object=Integ.obj, file="Seurat.UMAP_label.RDS")
Integ.obj <- readRDS("Seurat.UMAP_label.RDS")
levels(Integ.obj)


DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=3, repel=F) + NoLegend()

### Large cluster
levels(Integ.obj)

new.cluster.ids <- c("Astro", "Astro", "Astro", "Astro", "PQ", "PQ", "PQ", "PQ",
                     "Micro", "Micro", "Endo", "Endo", "Oligo", "Micro", "Oligo", "Oligo", "Oligo",
                     "InN", "InN", "OPC", "OPC", "ExN", "ExN", "InN", "InN", "InN", "InN", "InN",
                     "InN", "InN", "ExN", "ExN", "ExN", "ExN", "ExN", "ExN", "ExN", "ExN", "ExN",
                     "ExN", "ExN", "ExN")

names(new.cluster.ids) <- levels(Integ.obj)
Integ.obj <- RenameIdents(Integ.obj, new.cluster.ids)
Integ.obj@meta.data$"Largecluster" <- as.factor(Integ.obj@active.ident)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=3, repel=T) + NoLegend()

my_cols <- c('Astro1'='#ed8b28', 'Astro2'='#e33b50', 'Astro3'='#e4622a', 'Astro4'='#f5a97a',
             'OPC1'='#FFB500', 'OPC2'='#cf8104',
             'Oligo1'='#c2a442', 'Oligo2'='#ba8a1a', 'Oligo3'='#967341', 'Oligo4'='#A05837',
             'Micro'='#A4E804', 'Micro_Phago'='#3aa312', 'Macro'='#02705f',
             'Endo1'='#2e75d1', 'Endo2'='#81ace3',
             'PQ1'='#ccc2c2', 'PQ2'='#a3a3a3', 'PQ3'='#8c8484', 'PQ4'='#d4d2c1',
             'ExN1_L2-4'='#c8bdf0', 'ExN2_L2-4'='#918fe3', 'ExN3_L4-5'='#7e7cd6', 'ExN4_L2-6'='#b06ae6',
             'ExN5_L2-4'='#8e6cd4', 'ExN6_L4-5'='#5d58e0', 'ExN7_L5-6'='#6f7fa6', 'ExN8_L4-5'='#5463a8',
             'ExN9_L4-5'='#b19bde', 'ExN10_L5-6'='#574999', 'ExN11_L4-5'='#8988b8', 'ExN12_L5-6'='#8d4dbf',
             'ExN13_L6'='#4a1b6e',
             'InN1_MGE_SST'='#e07b91', 'InN2_CGE_VIP'='#FFAA92', 'InN3_MGE'='#C8A1A1', 'InN4_CGE_LAMP5'='#f7c7a1',
             'InN5_CGE_CALB2'='#c486a5', 'InN6_CGE_NDNF'='#b56475', 'InN7_MGE_PVALB'='#c9908b', 'InN8_MGE_NOS1'='#7B4F4B',
             'InN9_CGE'='#ed91bf')

my_cols <- c('ExN'='#918fe3', 'InN'='#ed91bf', 'PQ'='#d4d2d2',
             'Astro'='#ed8b28', 'OPC'='#FFB500', 'Oligo'='#967341', 'Micro'='#A4E804', 
             'Endo'='#2e75d1')
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]

pdf("UMAP_Final2.pdf", width = 7.126, height = 7)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=4, repel=T, cols=my_cols2) + NoLegend()
dev.off()


pdf("UMAP_Large.pdf", width = 7.125, height = 7)
DimPlot(Integ.obj, reduction = "umap.harmony", label = TRUE, label.size=3, repel=T, cols=my_cols2) + NoLegend()
dev.off()
#levels(ExN.sub) <- c("ExN10_L6", "ExN9_L5-6", "ExN7_L4-6", "ExN1_L5", "ExN5_L4-5", "ExN6_L4-5",
#                     "ExN8_L4-5", "ExN4_L2-6", "ExN2_L2-5", "ExN3_L2-4")
#new.cluster.ids <- c("InN5", "InN7_PVALB", "InN1_SST", "InN8_NOS1",
#                     "InN2", "InN3_VIP", "InN4_LAMP5", "InN6_SST")


