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
library(RColorBrewer)
library(limma)

full.combined <- readRDS(file="Seurat.Final.RDS")
levels(full.combined)

### Pseudobulk_Counts_Extract
# Preprocessing
metadata <- full.combined@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
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
metadata$group[which(str_detect(metadata$cells, "^CBC_C"))] <- "CBCC"
metadata$group[which(str_detect(metadata$cells, "^CBC_M"))] <- "CBCM"
full.combined@meta.data <- metadata

View(full.combined@meta.data)

######### get Pseudo-bulk normalized counts #############
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

A <- rep(c("CBCC", "CBCM"), each=6)
group_id <- factor(A, levels=c("CBCC", "CBCM"))

y <- DGEList(cts, remove.zeros = TRUE)
keep <- filterByExpr(y, group=group_id, min.count=1, min.total.count=13)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
CPM <- cpm(y)

write.csv(CPM, "Pseudo-bulk_normalized_counts_edgeR1.csv")



### Pseudobulk_DEG_EdgeR
DefaultAssay(full.combined) <- "RNA"
counts <- full.combined@assays$RNA@counts
metadata <- full.combined@meta.data[, c(13,14,19)]
#View(metadata)

metadata$sample_id <- factor(full.combined@meta.data$sample)
metadata$group_id <- factor(full.combined@meta.data$group)
metadata$cluster_id <- factor(full.combined@active.ident)

metadata <- metadata[, 4:6]
#View(metadata)

metadata$Paired <- NA
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC1"))] <- "P1"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC2"))] <- "P2"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC3"))] <- "P3"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC4"))] <- "P4"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC5"))] <- "P5"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCC6"))] <- "P6"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM1"))] <- "P1"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM2"))] <- "P2"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM3"))] <- "P3"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM4"))] <- "P4"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM5"))] <- "P5"
metadata$Paired[which(str_detect(metadata$sample_id, "^CBCM6"))] <- "P6"
metadata$Paired <- factor(metadata$Paired)


sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)
sce
groups <- colData(sce)[, c("cluster_id", "sample_id")]

assays(sce)
dim(counts(sce))
counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))

sce <- sce[rowSums(counts(sce)) > 12, ]

kids <- purrr::set_names(levels(sce$cluster_id))
kids

nk <- length(kids)
nk

sids <- purrr::set_names(levels(sce$sample_id))
sids

ns <- length(sids)
ns

table(sce$sample_id)

n_cells <- as.numeric(table(sce$sample_id))

m <- match(sids, sce$sample_id)

ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) %>%
  select(-"cluster_id")
ei


# Aggregate the counts per sample_id and cluster_id
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)
dim(pb)
pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
kable(table(sce$cluster_id, sce$sample_id))



# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "Paired")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id, Paired)

metadata        
metadata$cluster_id <- factor(metadata$cluster_id)
metadata$sample_id <- factor(metadata$sample_id)
metadata$group_id <- factor(metadata$group_id)
metadata$Paired <- factor(metadata$Paired)
metadata        

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters

clusters[1]


# construct SCE of pseudo-bulk counts for only select clusters
# If you are interested in all clusters AND you have the same samples represented in each cluster you can just use pb

# Create a character vector of the clusters to use for DE
keepClusters <-as.character(c(1:15))

# Subset the sce object
(interestingClusters <- SingleCellExperiment(assays = pb[1:15]))


res <- lapply(1:15, function(k) {
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[k]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  Paired <- factor(cluster_metadata$Paired)
  group_id <- factor(cluster_metadata$group_id, levels=c("CBCC", "CBCM"))
  design <- model.matrix(~ Paired + group_id)
  y <- SingleCellExperiment(assays = pb[k])
  y <- assays(y)[[1]]
  y <- DGEList(y, remove.zeros = TRUE)
  #keep <- rowSums(y$counts) > 12
  keep <- filterByExpr(y, group=group_id, min.count=1, min.total.count=13)
  y <- y[keep, ,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=T)
  fit <- glmQLFTest(fit)
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})


# Results filtering & overview

# filter pvalue < 0.05, |logFC| > 0.3785 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_val < 0.05, abs(logFC) > 0.3785) %>% 
                    dplyr::arrange(p_val))

# filter FDR<0.1
res_FDR <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.1) %>% 
                    dplyr::arrange(p_adj))


## Count the number of differential findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=clusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(interestingClusters) * 100, digits =2)) %>%  kable()

for(cluster in 1:length(keepClusters)){
  # Full results
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0("Pairwise/edgeR_", clusters[cluster], "_Paired_allGenes.csv"), quote=F, row.names = F)
  
  # Sig genes
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0("Pairwise/edgeR_", clusters[cluster], "_Paired_sigGenes.csv"), quote=F, row.names = F)
  
  # FDR < 0.1
  out <- res_FDR[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0("Pairwise/FDR/edgeR_", clusters[cluster], "_Paired_FDRGenes.csv"), quote=F, row.names = F)
}