library(WGCNA)
# enrichment graph ----------------------------------------------------------------------

res_My15 = c()
for (i in 1:11) {
  a = read.csv(spaste("Lake",i,".csv"))
  b = read.csv("My15.csv")
  c = subset(a$Gene, is.element(a$Gene, b$Gene) == T)
  aa = length(c)
  bb = length(a$Gene) - aa
  cc = length(b$Gene) - aa
  dd = 1733 - (aa + bb + cc) #1733
  ee = matrix(c(aa, bb, cc, dd), nrow = 2)
  res_My15[i] = fisher.test(ee, alternative = "greater")$p.value
  BH = p.adjust(res_My15, "BH")
}
res_My15
log_My15 = -log10(BH)
log_My15
write.csv(log_My15, "Enrich_cluster_My15.csv", row.names = F)

Total = cbind(res_My1, res_My2, res_My3, res_My4, res_My5, res_My6, res_My7, res_My8, res_My9, res_My10,
              res_My11, res_My12, res_My13, res_My14, res_My15)
Total_BH = p.adjust(Total, "BH")
log_Total = -log10(Total_BH)
log_Total
write.csv(Total_BH, "Enrich_cluster_pvalue.csv", row.names = F)
write.csv(log_Total, "Enrich_cluster_Total.csv", row.names = F)
