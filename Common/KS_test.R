library(Seurat)
library(dplyr)

full.combined <- readRDS(file="Seurat.Final.RDS")
levels(full.combined)
levels(full.combined) <- c('Astro1', 'Astro2', 'OPC', 'Oligo', 'Micro', 'Endo',
                           'Gran1', 'Gran2', 'Gran3', 'Gran4', 'Gran5', 'Gran6', 'Gran7',
                           'Purk1', 'Purk2')

full.combined$cell_type <- Idents(full.combined)

# PEM
cluster_markers <- AverageExpression(full.combined, group.by = "cell_type", assay = "RNA")$RNA
cluster_totals <- colSums(cluster_markers)

calculate_pem <- function(gene, cluster_markers, cluster_totals) {
  gene_expression <- cluster_markers[gene, ]
  pem <- gene_expression / cluster_totals
  return(pem / sum(pem))
}

all_genes <- rownames(cluster_markers)
pem_matrix <- sapply(all_genes, calculate_pem, cluster_markers, cluster_totals)
pem_matrix <- t(pem_matrix)
pem_matrix <- pem_matrix[complete.cases(pem_matrix), ]
#write.csv(pem_matrix, "PEM_CBC.csv")


module_genes <- read.csv("kME0.7E_M40.csv")
module_genes <- module_genes$Gene # Module gene set
selected_genes <- module_genes

## Background genes
# 1. othermodule genes
othermodule_gene <- read.csv("All_module_gene_kME0.7.csv")
background_genes <- setdiff(othermodule_gene$Gene, module_genes) # Background gene set

# 2. Filtering
selected_genes_filtered <- selected_genes[selected_genes %in% rownames(pem_matrix)]
background_genes_filtered <- background_genes[background_genes %in% rownames(pem_matrix)]

# 3. KS test
ks_results <- list()

for (cluster in colnames(pem_matrix)) {
  selected_pem <- pem_matrix[selected_genes_filtered, cluster]
  background_pem <- pem_matrix[background_genes_filtered, cluster]
  ks_test <- ks.test(selected_pem, background_pem, exact = TRUE, alternative = "less")
  ks_results[[cluster]] <- list(
    raw_p_value = ks_test$p.value,
    statistic = ks_test$statistic
  )
}

# 4. BH adj
raw_p_values <- sapply(ks_results, function(x) x$raw_p_value)
adjusted_p_values <- p.adjust(raw_p_values, method = "BH")


for (cluster in colnames(pem_matrix)) {
  cat("\nCluster:", cluster, "\n")
  cat("Raw p-value:", ks_results[[cluster]]$raw_p_value, "\n")
  cat("Adjusted p-value (BH):", adjusted_p_values[cluster], "\n")
  cat("KS Statistic:", ks_results[[cluster]]$statistic, "\n")
}

KStest <- data.frame(
  Cluster = names(raw_p_values), 
  D_value = sapply(ks_results, function(x) x$statistic),
  P_value = raw_p_values)
KStest$BH_adjusted <- p.adjust(KStest$P_value, method = "BH")
format(KStest$P_value, scientific = TRUE)
write.csv(KStest, "M40_KS.csv", row.names = FALSE)


### Dot Plot
library(ggplot2)
library(forcats)

C <- read.csv ("CBC_each_KS.csv")

pdf("CBC_each_KS.pdf", width = 6.5, height = 4)
C$module <- factor(C$module, levels = c("M49", "M40", "M32", "M24", "M10", "M9"))
C$Cluster <- factor(C$Cluster, levels = c("Astro1", "Astro2", "OPC", "Oligo", "Micro", "Endo",
                                          "Gran1", "Gran2", "Gran3", "Gran4", "Gran5", "Gran6", "Gran7",
                                          "Purk1", "Purk2"))
sp2<-ggplot(C, aes(x=Cluster, y=module, size=-log10(pvalue))) +
  scale_size_continuous(range = c(1,11)) + 
  geom_point(aes(fill = ifelse(pvalue < 0.01, "True", "False")), 
             colour = "black", pch = 21) +
  scale_fill_manual(values = c("True" = "red", "False" = "gray"),
                    limits = c("True", "False")) +
  labs(title="CBC each clusters (KS test)", x=NULL, y=NULL, fill="Significant", size = "-Log10 (FDR)")+
  guides(colour=guide_colourbar(order=1),
         size=guide_legend(order=2))
sp2+theme_bw()+
  theme(axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1))+
  theme(axis.text.y = element_text(colour="black", size=13))+
  theme(plot.title = element_text(size=17, colour="black", face="bold", hjust=0.5))+
  theme(panel.border = element_blank(),
        legend.position = "bottom")

dev.off()
