library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(tidyr)
library(RColorBrewer)

full.combined <- readRDS(file="Seurat.Final.RDS")
levels(full.combined)
View(full.combined@meta.data)

metadata <- full.combined@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C1_"))] <- "dlPFCC1"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C2_"))] <- "dlPFCC2"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C3_"))] <- "dlPFCC3"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C4_"))] <- "dlPFCC4"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C5_"))] <- "dlPFCC5"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_C6_"))] <- "dlPFCC6"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M1_"))] <- "dlPFCM1"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M2_"))] <- "dlPFCM2"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M3_"))] <- "dlPFCM3"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M4_"))] <- "dlPFCM4"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M5_"))] <- "dlPFCM5"
metadata$sample[which(str_detect(metadata$cells, "^dlPFC_M6_"))] <- "dlPFCM6"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C1_"))] <- "CBCC1"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C2_"))] <- "CBCC2"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C3_"))] <- "CBCC3"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C4_"))] <- "CBCC4"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C5_"))] <- "CBCC5"
metadata$sample[which(str_detect(metadata$cells, "^CBC_C6_"))] <- "CBCC6"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M1_"))] <- "CBCM1"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M2_"))] <- "CBCM2"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M3_"))] <- "CBCM3"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M4_"))] <- "CBCM4"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M5_"))] <- "CBCM5"
metadata$sample[which(str_detect(metadata$cells, "^CBC_M6_"))] <- "CBCM6"
full.combined@meta.data <- metadata

metadata$group <- NA
metadata$group[which(str_detect(metadata$cells, "^dlPFC_C"))] <- "dlPFCC"
metadata$group[which(str_detect(metadata$cells, "^dlPFC_M"))] <- "dlPFCM"
metadata$group[which(str_detect(metadata$cells, "^CBC_C"))] <- "CBCC"
metadata$group[which(str_detect(metadata$cells, "^CBC_M"))] <- "CBCM"
full.combined@meta.data <- metadata

View(full.combined@meta.data)

DefaultAssay(full.combined) <- "RNA"
cts <- AggregateExpression(full.combined,
                           group.by = "sample",
                           slot="counts",
                           return.seurat = FALSE)

cts <- cts$RNA
dim(cts)

cts[1:10, 1:10]
cts <- as.data.frame(cts)
write.csv(cts, "Pseudo-bulk_counts.csv")

cts <- as.matrix(cts)

A <- rep(c("CBCC", "CBCM", "dlPFCC", "dlPFCM"), each=6)
group_id <- factor(A, levels=c("CBCC", "CBCM", "dlPFCC", "dlPFCM"))

y <- DGEList(cts, remove.zeros = TRUE)
keep <- filterByExpr(y, group=group_id, min.count=1, min.total.count=25)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
CPM <- cpm(y)

write.csv(CPM, "Pseudo-bulk_normalized_counts.csv")

### Pearson's correlation

RNA_cts <- read.csv("mRNA_FPKM_quan.csv")
snRNA_cts <- read.csv("Pseudo-bulk_normalized_counts_ID.csv")
head(RNA_cts)
dim(RNA_cts)
head(snRNA_cts)
dim(snRNA_cts)

RNA <- pivot_longer(RNA_cts,
                    cols = -ID,
                    names_to = "Sample",
                    values_to = "RNA_exp")
snRNA <- pivot_longer(snRNA_cts,
                      cols = -ID,
                      names_to = "Sample",
                      values_to = "snRNA_exp")

A <- inner_join(RNA_cts, snRNA_cts,
                by= c("ID"))

B <- inner_join(RNA, snRNA,
                by= c("ID", "Sample")) %>%
    group_by(ID) %>%
    summarise(rho = cor(RNA_exp, snRNA_exp),
              pvalue = cor.test(RNA_exp, snRNA_exp)$p.value)

write.csv(B, "RNA_snRNA_correlation.csv", row.names = FALSE)
