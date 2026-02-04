library(preprocessCore)


FPKM <- read.csv("Before_Quan.csv")
FPKM_mat <- data.matrix(FPKM)

FPKM_nor <- normalize.quantiles(FPKM_mat[,-c(1:4)], copy=TRUE)

write.csv(FPKM_nor, "After_Quan.csv")
