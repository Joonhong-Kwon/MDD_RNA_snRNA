library(rgl)
library(magrittr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tibble)
library(svglite)

# 3D
data <- read.csv("Log2_FPKM_MDS.csv")
dim(data)

mRF = t(data[1:15738,-1]) 
Region = cbind(t(data[15739,-1]), t(data[1:15738,-1]))   # Region
Disease = cbind(t(data[15740,-1]), t(data[1:15738,-1]))   # State
dim(mRF)

d = dist(mRF) # euclidean distances between the rows
fit = cmdscale(d, eig = T, k = 3) # k is the number of dim
fit$eig

eig_values <- fit$eig
eig_values <- eig_values[eig_values > 0]

variance_explained <- eig_values / sum(eig_values) * 100
variance_explained

head(fit) # view results

# plot solution 
dim1 <- fit$points[,1]
dim2 <- fit$points[,2]
dim3 <- fit$points[,3]

mds <- mRF %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2", "Dim.3")

# groups
group_label = c(rep("dlPFC",12), rep("CA1",12), rep("DG",12), rep("CBC",12))
group_label1 = c(rep("Con",6), rep("MDD",6), rep("Con",6), rep("MDD",6),rep("Con",6), rep("MDD",6),rep("Con",6), rep("MDD",6))
mds <- mds %>%
  mutate(Disease = as.factor(group_label1))
mds <- mds %>%
  mutate(Region = as.factor(group_label))
mds$Region <- factor(mds$Region, levels = c("dlPFC", "CA1", "DG", "CBC"))
mds$Disease <- factor(mds$Disease, levels = c("Con", "MDD"))
ldash = rownames(mRF)

plot3d(dim1, dim2, dim3, col=as.integer(mds$Region), size= 13)
text3d(dim1, dim2, dim3+5, ldash, cex=0.8, col=as.integer(mds$Region))
legend3d("topright", legend = c("dlPFC", "CA1", "DG", "CBC"), pch=16, cex = 2, col = palette("default"))

rgl.postscript('3D Plot3.pdf', fmt = 'pdf')


# 2D
data <- read.csv("Log2_FPKM_MDS.csv")
dim(data)

mRF = t(data[1:15738,-1]) 
Region = cbind(t(data[15739,-1]), t(data[1:15738,-1]))   # Region
Disease = cbind(t(data[15740,-1]), t(data[1:15738,-1]))   # State
dim(mRF)

d = dist(mRF) # euclidean distances between the rows
fit = cmdscale(d, eig = T, k = 2) # k is the number of dim

head(fit) # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

# Cmpute MDS
mds <- mRF %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

# groups
group_label = c(rep("dlPFC",12), rep("CA1",12), rep("DG",12), rep("CBC",12))
group_label1 = c(rep("Con",6), rep("MDD",6), rep("Con",6), rep("MDD",6),rep("Con",6), rep("MDD",6),rep("Con",6), rep("MDD",6))
group_label2 = rep(c("C1", "C2", "C3", "C4", "C5", "C6", "M1", "M2", "M3", "M4", "M5", "M6"),4)
group_label2
mds <- mds %>%
  mutate(Disease = as.factor(group_label1))
mds <- mds %>%
  mutate(Sample = as.factor(group_label2))
mds <- mds %>%
  mutate(Region = as.factor(group_label))

# Plot and color by groups
pdf("HUman_MDD_2D_MDS_plot_Region.pdf", width = 5, height = 3.9)
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          color = "Region",
          #palette = c("#5993bd", "#d9346e"),
          palette = c("#b05394", "#EBB035", "forestgreen", "#06A2CB"),
          shape = "Disease",
          size = 3,
          ellipse = TRUE,
          ellipse.level = 0.95,
          ellipse.type = "confidence",
          font.label = c(8,"plain","black"),
          legend="right",
          label=NULL,  #rownames(mRF)
          repel = TRUE)
dev.off()