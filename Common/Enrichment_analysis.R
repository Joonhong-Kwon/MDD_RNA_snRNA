library(WGCNA)
### DEG Module enrichment
res_CA1_UP = c()
for (i in 1:49) {
  a = read.csv(spaste("kME0.7E_M",i,".csv"))
  b = read.csv("Limma_trend_CA1_DEG_paired_0.05_1.3_up.csv")
  c = subset(a$ID, is.element(a$ID, b$ID) == T)
  aa = length(c)
  bb = length(a$ID) - aa
  cc = 687 - aa
  dd = 15738 - (aa + bb + cc)
  ee = matrix(c(aa, bb, cc, dd), nrow = 2)
  res_CA1_UP[i] = fisher.test(ee, alternative = "greater")$p.value
  BH = p.adjust(res_CA1_UP, "BH")
}
res_CA1_UP
log_CA1_UP = -log10(BH)
log_CA1_UP
write.csv(log_CA1_UP, "enrich_CA1_1.3_UP.csv", row.names = F)


Total = cbind(res_PFC_UP, res_PFC_DOWN, res_CA1_UP, res_CA1_DOWN, res_DG_UP, res_DG_DOWN, res_CBC_UP, res_CBC_DOWN)
Total_BH = p.adjust(Total, "BH")
Total_BH
log_Total = -log10(Total_BH)
log_Total
write.csv(Total_BH, "enrich_DEG_1.3_Total_pvalue.csv", row.names = F)
write.csv(log_Total, "enrich_DEG_1.3_Total.csv", row.names = F)


### Cell type Module enrichment
H <- read.csv("Cell type_H.csv")

# Cell type enrichment analysis for whole genes (kME > 0.7)-------------------------------------------------
a = read.csv("All_module_gene_kME0.7_UPPER.csv")
b = merge(x=a, y=H, by='Gene', all.x=T)
colnames(b) = c("Gene", "ID", "ModuleLable", "ModuleColor", "kME", "p.kMEME", "Fetal Astrocyte", "Mature Astrocyte", "Neuron", 
                "Oligodendrocyte", "Microglia/Macrophage", "Endothelial Cell", "Max", "Cell Type")
write.csv(b, file = "Cell type_H_kME0.7.csv", 
          row.names = FALSE, quote = FALSE)

# gene expression of each cell type (kME > 0.7)---------------------------------------------------------

c = read.csv("Cell type Geschwind.csv")
c = as.matrix(c)
N_Neuron = c[,1]
N_Astro = c[,2]
N_Micro = c[,3]
N_Endo = c[,4]
N_Oligo = c[,5]
N_GABA = c[,6]
N_GLUT = c[,7]

CBC = read.csv("Log2_FPKM_quan_CBC.csv")

# Ben Barres Human from Geschwind
# Astro, Neuron, Oligo, Micro, Endo, GABA, GLUT
Astro_CBC = subset(CBC, is.element(CBC$Gene, N_Astro) == T)
dim(Astro_CBC)
write.csv(Astro_CBC, "Astro_CBC_GW.csv", row.names = F)


Astro_CBC = cbind(Astro_CBC[,1:2], Astro_CBC[,9:14]-Astro_CBC[,3:8])
dim(Astro_CBC)
names(Astro_CBC) = c("Gene","ID",1,2,3,4,5,6)
View(Astro_CBC)
write.csv(Astro_CBC, "Astro_CBC_box_GW.csv", row.names = F)


# enrichment graph ----------------------------------------------------------------------
# Astro, Neuron, Oligo, Micro, Endo, GABA, GLUT
res_Astro = c()
for (i in 1:49) {
  a = read.csv(spaste("kME0.7E_M",i,".csv"))
  b = read.csv("Astro_CBC_box_GW.csv")
  c = subset(a$ID, is.element(a$ID, b$ID) == T)
  aa = length(c)
  bb = length(a$ID) - aa
  cc = 217 - aa
  dd = 15738 - (aa + bb + cc)
  ee = matrix(c(aa, bb, cc, dd), nrow = 2)
  res_Astro[i] = fisher.test(ee, alternative = "greater")$p.value
  BH = p.adjust(res_Astro, "BH")
}
res_Astro
log_Astro = -log10(BH)
log_Astro
write.csv(log_Astro, "enrich_Astro.csv", row.names = F)


Total = cbind(res_Neuron, res_GABA, res_GLUT, res_Astro, res_Oligo, res_Micro, res_Endo)
Total_BH = p.adjust(Total, "BH")
log_Total = -log10(Total_BH)
log_Total
write.csv(Total_BH, "enrich_Total_pvalue.csv", row.names = F)
write.csv(log_Total, "enrich_Total.csv", row.names = F)
