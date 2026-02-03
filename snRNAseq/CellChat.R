library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(tidyverse)
library(Matrix)
library(WGCNA)
library(ggplot2)
library(CellChat)
library(patchwork)
library(parallel)
library(future)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

set.seed(1234)

# Data pre-processing
full.combined <- readRDS(file="Seurat.Final.RDS")
levels(full.combined)

order <- c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
           "Purk1", "Purk2", "Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo")
levels(full.combined)
full.combined$Newcluster <- factor(full.combined$Newcluster, levels=order)
levels(full.combined)

### Comparison analysis
CBC.Con <- subset(x=full.combined, subset=group=="CBC_C")
CBC.MDD <- subset(x=full.combined, subset=group=="CBC_M")

cellchat.Con <- createCellChat(object = CBC.Con, group.by = "Newcluster", assay = "RNA")
cellchat.MDD <- createCellChat(object = CBC.MDD, group.by = "Newcluster", assay = "RNA")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat.Con@DB <- CellChatDB.use
cellchat.MDD@DB <- CellChatDB.use

cellchat.Con <- subsetData(cellchat.Con) # subset the expression data of signaling genes for saving computation cost
cellchat.MDD <- subsetData(cellchat.MDD)
future::plan("multicore", workers = 32) # do parallel

cellchat.Con <- identifyOverExpressedGenes(cellchat.Con)
cellchat.Con <- identifyOverExpressedInteractions(cellchat.Con)
cellchat.Con <- projectData(cellchat.Con, PPI.human)

cellchat.MDD <- identifyOverExpressedGenes(cellchat.MDD)
cellchat.MDD <- identifyOverExpressedInteractions(cellchat.MDD)
cellchat.MDD <- projectData(cellchat.MDD, PPI.human)


options(future.globals.maxSize = 8000 * 1024^2)
cellchat.Con <- computeCommunProb(cellchat.Con, raw.use = TRUE, population.size = TRUE)
cellchat.MDD <- computeCommunProb(cellchat.MDD, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.Con <- filterCommunication(cellchat.Con, min.cells = 10)
cellchat.MDD <- filterCommunication(cellchat.MDD, min.cells = 10)

cellchat.Con <- computeCommunProbPathway(cellchat.Con)
cellchat.MDD <- computeCommunProbPathway(cellchat.MDD)

cellchat.Con <- aggregateNet(cellchat.Con)
cellchat.MDD <- aggregateNet(cellchat.MDD)

groupSize <- as.numeric(table(cellchat.Con@idents))
groupSize <- as.numeric(table(cellchat.MDD@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.Con@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions - CON")
netVisual_circle(cellchat.MDD@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions - MDD")
netVisual_circle(cellchat.Con@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength - CON")
netVisual_circle(cellchat.MDD@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength - MDD")

cellchat.Con <- netAnalysis_computeCentrality(cellchat.Con, slot.name = "netP") 
cellchat.MDD <- netAnalysis_computeCentrality(cellchat.MDD, slot.name = "netP") 

object.list <- list(CBC_C = cellchat.Con, CBC_M = cellchat.MDD)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


my_cols <- c('Gran1'='#c8bdf0', 'Gran2'='#918fe3', 'Gran3'='#7e7cd6', 'Gran4'='#b06ae6',
             'Gran5'='#8e6cd4', 'Gran6'='#5d58e0', 'Gran7'='#6f7fa6',
             'Purk1'='#e07b91', 'Purk2'='#ed91bf',
             'Astro1'='#ed8b28', 'Astro2'='#e33b50',
             'OPC'='#FFB500', 'Oligo'='#967341',
             'Micro' = '#A4E804',
             'Endo' = '#2e75d1')


gg1 <- netVisual_heatmap(cellchat, color.use = my_cols,
                         col.show=c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                                    "Purk1", "Purk2", "Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo"),
                         row.show=c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                                    "Purk1", "Purk2", "Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo"))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = my_cols,
                         col.show=c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                                    "Purk1", "Purk2", "Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo"),
                         row.show=c("Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                                    "Purk1", "Purk2", "Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo"))
#> Do heatmap based on a merged object
gg1 + gg2


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
object.list[[2]]

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

saveRDS(cellchat, file = "cellchat_comparison_MDD_consider_proportion.rds")
#cellchat <- readRDS("cellchat_comparison_MDD_broadv2_cluster.rds")
cellchat <- readRDS("cellchat_comparison_MDD_allcluster.rds")


############################### specific cluster ##########################################
# Data pre-processing
cellchat <- readRDS("cellchat_comparison_MDD_Micro_Purk2.rds")
sort(unique(cellchat@idents$joint))
cellchat.sub <- subsetCellChat(object=cellchat, idents.use=c("Micro", "Purk2"))
sort(unique(cellchat.sub@idents$joint))

tiff(filename = "Bargraph_Micro_Purk2.tif", width = 6, height = 4, units = "in", res=300)
pdf("Bargraph_Micro_Purk2.pdf", width = 5, height = 4)
gg1 <- compareInteractions(cellchat.sub, show.legend = F, group = c(1,2), color.use = c("#c5ddde", "#E2A7C1"),
                           title.name = "Micro-Purk2 interaction")
gg2 <- compareInteractions(cellchat.sub, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#c5ddde", "#E2A7C1"))
gg1 + gg2
dev.off()

tiff(filename = "Bargraph2_Micro_Purk2.tif", width = 7.5, height = 4, units = "in", res=300)
pdf("Bargraph2_broad_Micro_Purk.pdf", width = 3.5, height = 4)
gg1 <- rankNet(cellchat.sub, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#c5ddde", "#E2A7C1"))
gg2 <- rankNet(cellchat.sub, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#c5ddde", "#E2A7C1"))
gg1 + gg2
dev.off()

pdf("Dotplot_Micro_Purk.pdf", width=8, height=7)
netVisual_bubble(cellchat.sub, sources.use = c(1:2), targets.use = c(1:2), comparison = c(1, 2), angle.x = 45)
dev.off()

pairLR.use_Con <- extractEnrichedLR(cellchat.sub, signaling = c("APP", "BMP", "COLLAGEN", "FN1", "LAMININ", "PTPRM",
                                                                "SPP1"))
pairLR.use_MDD <- extractEnrichedLR(cellchat.sub, signaling = c("SPP1"))


pdf("DotPlot2_Micro_Purk.pdf", width = 9, height = 6)
gg1 <- netVisual_bubble(cellchat.sub, sources.use = c(1:2), targets.use = c(1:2),  comparison = c(1, 2),
                        max.dataset = 2, title.name = "Increased signaling in MDD", angle.x = 45, remove.isolate = F,
                        color.text = c("#65A4A7", "#BD3E77"), pairLR.use = pairLR.use_Con)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.sub, sources.use = c(1:2), targets.use = c(1:2),  comparison = c(1, 2),
                        max.dataset = 1, title.name = "Decreased signaling in MDD", angle.x = 45, remove.isolate = F,
                        color.text = c("#65A4A7", "#BD3E77"))
#> Comparing communications on a merged object
gg1 #+ gg2
dev.off()

saveRDS(cellchat.sub, file = "cellchat_comparison_MDD_Micro_Purk2.rds")
