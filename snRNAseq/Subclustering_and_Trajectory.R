library(dittoSeq)
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

set.seed(1234)
full.combined <- readRDS("Seurat.Final.RDS")

levels(full.combined)
DimPlot(full.combined, reduction = "umap", label = TRUE, repel=F) + NoLegend()

DefaultAssay(full.combined) <- "integrated"

Oligo.sub <- subset(x=full.combined, idents=c("OPC1", "OPC2", "Oligo1", "Oligo2", "Oligo3", "Oligo4"))
levels(Oligo.sub)

### Subclustering
DefaultAssay(Oligo.sub) <- "RNA"
Oligo.sub <- FindVariableFeatures(Oligo.sub, selection.method = "vst", nfeatures = 2000)

DefaultAssay(Oligo.sub) <- "integrated"
Oligo.sub <- ScaleData(Oligo.sub, verbose = FALSE, vars.to.regress=c("nCount_RNA"))

Oligo.sub <- RunPCA(Oligo.sub, npcs = 30, verbose = FALSE)
ElbowPlot(object=Oligo.sub, ndims=30)
Oligo.sub <- JackStraw(object=Oligo.sub, num.replicate=100, dims=30)
Oligo.sub <- ScoreJackStraw(object=Oligo.sub, dims=1:30)
JackStrawPlot(object=Oligo.sub, dims=1:9)

Oligo.sub <- RunUMAP(Oligo.sub, reduction = "pca", dims = 1:9)
Oligo.sub <- FindNeighbors(Oligo.sub, reduction = "pca", dims = 1:9)
Oligo.sub <- FindClusters(Oligo.sub, resolution = 0.3)

col = c('#FFB500','#cf8104', '#c2a442', '#ba8a1a', '#967341', '#A05837')
levels(Oligo.sub)
Oligo.sub <- SetIdent(Oligo.sub, value = "Newcluster")
levels(Oligo.sub) <- c("OPC1", "OPC2", "Oligo1", "Oligo2", "Oligo3", "Oligo4")

pdf("UMAP_Oligo_sub.pdf", width = 6, height = 4.981)
DimPlot(Oligo.sub, reduction = "umap", label = TRUE, repel=F, cols = col)
dev.off()

pdf("UMAP_Oligo_sub_split_group.pdf", width = 9.787, height = 5)
DimPlot(object=Oligo.sub, reduction="umap", split.by="group", label = TRUE, repel=T, cols = col)
dev.off()

dittoBarPlot(
  object = Oligo.sub,
  var = "integrated_snn_res.0.3",
  group.by = "sample")


DimPlot(object=Oligo.sub3, reduction="umap", split.by="group", ncol=2)
DimPlot(object=Oligo.sub1, reduction="umap", split.by="sample", ncol=6)

DefaultAssay(Oligo.sub) <- "RNA"
Oligo <- c("MAG", "MOG", "MBP", "QKI", "PCDH9", "PLP1", "CNP", "CLDN11", "CNTNAP2",
           "LPAR6", "ADAM28", "ARHGAP24", "CTSS", "CD74", "C3", "P2RY12",
           "OLIG1", "OLIG2", "CSPG4", "PDGFRA", "SOX6", "VCAN", "DSCAM", "PCDH15")
Oligo1 <- c("PDGFRA", "PCDH15", "GPR17", "SOX4", "MOG", "CDH7", "FRY", "OPALIN", "GSN",
            "APOE", "CD74", "P2RY12", "CDH20", "RBFOX1", "KLK6", "GJB1", "LURAP1L-AS1", "CDH19", "MAG")

DotPlot(Oligo.sub, features = Oligo1) +
  theme(axis.text.x = element_text(size=11),
        axis.text.x.bottom = element_text(angle=90, colour = "black", face=3),
        legend.title=element_text(size=10))+
  guides(color = guide_colorbar(title = "Average\nExpression"),
         size = guide_legend("Percent\nExpressed"),) +
  labs(x=NULL, y=NULL)

VlnPlot(Oligo.sub, features=Oligo1, stack=T, sort=F, flip=T) +
  theme(legend.position = "none")+
  theme(strip.text.y = element_text(face="italic", size=11),
        axis.text.x = element_text(size=11, angle=0, hjust=0.5)) +
  labs(x=NULL, y=NULL) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0)


### Cell proportion
count_table <- table(Oligo.sub1@meta.data$integrated_snn_res.0.2, Oligo.sub@meta.data$orig.ident)
count_table
write.csv(count_table, "Cell_count_table_PC9_res0.3.csv")


### DEG Oligo sub
DefaultAssay(Oligo.sub) <- "RNA"
seurat.markers <- FindAllMarkers(object=Oligo.sub, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
seurat.top10.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
seurat.top20.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

DefaultAssay(Oligo.sub) <- "integrated"
DoHeatmap(object=Oligo.sub, features=seurat.top10.genes$gene, size=4) + NoLegend()

### Save
write.csv(seurat.markers, "Seurat.FindAllMarkers_Oligo_sub_PC9_res0.3.csv")
write.csv(seurat.top10.genes, "Seurat.top10markers_Oligo_sub_PC9_res0.3.csv")
write.csv(seurat.top20.genes, "Seurat.top20markers_Oligo_sub_PC9_res0.3.csv")

saveRDS(object=Oligo.sub, file="Oligo_sub_PC9_res0.3.RDS")
Oligo.sub <- readRDS(file="Oligo_sub_PC9_res0.3.RDS")


### Trajectory

Oligo.sub3 <- readRDS(file="Oligo_sub_PC9_res0.3.RDS")

### Monocle3
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3
# ...1. Convert to cell_data_set object -----
DefaultAssay(Oligo.sub3) <- "RNA"
cds <- as.cell_data_set(Oligo.sub3)
cds

# to get cell metadata
colData(cds)
# to gene metadata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

# to get counts
counts(cds)

#...2. Cluster cells (using clustering infor from seurat's UMAP) -----
# let's use the clustering information have

# assign partitions
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

# Assign the cluster info
list_cluster <- Oligo.sub3@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- Oligo.sub3@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'Newcluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right") +
  scale_color_manual(values=c('#FFB500','#cf8104', '#c2a442', '#ba8a1a', '#967341', '#A05837'))

cluster.names <- plot_cells(cds,
                            color_cells_by = "integrated_snn_res.0.3",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  theme(legend.position = "right")


cluster.before.trajectory
cluster.before.trajectory | cluster.names
DimPlot(object=Oligo.sub3, reduction="umap", split.by="group", ncol=2)

# ...3. Learn trajectory graph -----
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = "integrated_snn_res.0.3",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# ...4. Order the cells in pseudotime -----
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) == 3]))

pdf("UMAP_trajectory_Newcluster.pdf", width = 6, height = 4.869)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           trajectory_graph_segment_size = 0.35,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)+
  #group_label_size = 5)+
  theme(axis.line.x=element_line(size=0.5))+
  theme(axis.line.y=element_line(size=0.5))+
  theme(axis.ticks=element_line(size=0.5,
                                colour = "black"))+
  theme(axis.text=element_text(colour = "black",
                               size=13))+
  theme(axis.title=element_text(size=15))
dev.off()

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle_pseudotime, reorder(integrated_snn_res.0.3, monocle_pseudotime, median), fill=integrated_snn_res.0.3)) +
  geom_boxplot()
ggplot(data.pseudo, aes(monocle_pseudotime, reorder(Newcluster, monocle_pseudotime, median), fill=Newcluster)) +
  geom_boxplot()

cds_subset_Con <- cds[,pData(cds)$group == "dlPFC_C"]
cds_subset_MDD <- cds[,pData(cds)$group == "dlPFC_M"]

# ...5. Finding genes that change as a function of pseudotime
deg_Oligo <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 8)
deg_Oligo_Con <- graph_test(cds_subset_Con, neighbor_graph = 'principal_graph', cores = 8)
deg_Oligo_MDD <- graph_test(cds_subset_MDD, neighbor_graph = 'principal_graph', cores = 8)

deg_Oligo1 <- deg_Oligo[deg_Oligo$status == "OK" & deg_Oligo$q_value < 0.01,]
deg_Oligo_Con1 <- deg_Oligo_Con[deg_Oligo_Con$status == "OK" & deg_Oligo_Con$q_value < 0.01,]
deg_Oligo_MDD1 <- deg_Oligo_MDD[deg_Oligo_MDD$status == "OK" & deg_Oligo_MDD$q_value < 0.01,]

write.csv(deg_Oligo, "PC9_Oligo_pseudogene_not_filter.csv")
write.csv(deg_Oligo_Con1, "PC9_Oligo_Con_pseudogene_filter.csv")
write.csv(deg_Oligo_MDD1, "PC9_Oligo_MDD_pseudogene_filter.csv")
