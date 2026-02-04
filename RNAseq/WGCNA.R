### WGCNA ###--------------------------------------------------------
library(doParallel)
library(WGCNA);

options(stringsAsFactors = FALSE);
HuMDD = read.csv("Log2_FPKM_quan.csv", as.is = T, header = T);
dim(HuMDD);
names(HuMDD);

datExpr0 = as.data.frame(t(HuMDD[,-c(1:4)]));
names(datExpr0) = HuMDD$Gene;
rownames(datExpr0) = names(HuMDD)[-c(1:4)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = read.csv("HuMDD_trait.csv");
dim(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
MDDSamples = rownames(datExpr);
traitRows = match(MDDSamples, traitData$Sample);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
Col1 = numbers2colors(datTraits$Region, signed = TRUE, centered = FALSE, colors = c("#538ADB", "#81B9B7", "#B6E0BB", "#E06966"));
Col2 = numbers2colors(datTraits$Disease, signed = TRUE, centered = FALSE, colors = c("white", "#69549B"));
Col3 = numbers2colors(datTraits$Diagnosis, signed = TRUE, centered = FALSE, colors = colorRampPalette(c("white", "#69549B"))(n = 100));
Col4 = numbers2colors(datTraits$Age, signed = FALSE, centered = TRUE, colors = colorRampPalette(c("white", "#8B650B"))(n = 100));
Col5 = numbers2colors(datTraits$Sex, signed = FALSE, centered = FALSE, colors = c("#FF6C61","#A3CAEC"));
Col6 = numbers2colors(datTraits$PMI, signed = FALSE, centered = FALSE, colors = colorRampPalette(c("white", "#4E82B6"))(n = 100));
Col7 = numbers2colors(datTraits$pH, signed = FALSE, centered = FALSE, lim=c(1,14), colors = colorRampPalette(c("white", "darkslategray4"))(n = 100));
Col8 = numbers2colors(datTraits$RIN, signed = FALSE, centered = FALSE, lim=c(0,10), colors = colorRampPalette(c("white", "grey62"))(n = 100));
Col9 = numbers2colors(datTraits$Suicide, signed = TRUE, centered = FALSE, colors = c("white","#D9418C"));
Col10 = numbers2colors(datTraits$Smoke, signed = TRUE, centered = FALSE, colors = c("white","#BDBDBD", "#353535"));
Col11 = numbers2colors(datTraits$Medication, signed = TRUE, centered = FALSE, colors = c("white", "#3bad66"));
Col12 = numbers2colors(datTraits$Toxicology, signed = TRUE, centered = FALSE, colors = c("white","#D3C64A"));

col = cbind(Col1, Col2, Col3, Col4, Col5, Col6, Col7, Col8, Col9, Col10, Col11, Col12)

# Plot the sample dendrogram and the colors underneath.
pdf("Sample dendrogram and trait heatmap.pdf", width=9, height=6)
plotDendroAndColors(sampleTree2, col, marAll = c(1,6,3,1), 
                    groupLabels = c("Region", "Disease", "Diagnosis", "Age", "Sex", "PMI", "pH", "RIN", "Suicide", "Smoke", "Medication", "Toxicology"),
                    main = "Sample dendrogram and trait heatmap",
                    cex.dendroLabels = 0.8)
dev.off()

save(datExpr, datTraits, file = "All sample-01-dataInput.RData")


##########################################################################################

powers = 5:25
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 2)
collectGarbage();

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# package default, Neuron paper & Duman paper
net = blockwiseModules(
  datExpr,
  blocks = NULL,
  power = 12,
  maxBlockSize = 30000,
  minModuleSize = 20,
  randomSeed = 12345,
  corType = "pearson",
  networkType = "signed", 
  TOMType = "signed",
  TOMDenom = "min",
  detectCutHeight = 0.995, 
  deepSplit = 2,
  maxCoreScatter = NULL, minGap = NULL,
  maxAbsCoreScatter = NULL, minAbsGap = NULL,
  pamStage = TRUE, pamRespectsDendro = TRUE,
  minCoreKME = 0.7, minCoreKMESize = 3,
  minKMEtoStay = 0.7,
  reassignThreshold = 1e-6,
  mergeCutHeight = 0.1,
  impute = TRUE, 
  getTOMs = NULL,
  saveTOMs = FALSE, 
  trapErrors = FALSE, numericLabels = FALSE,
  checkMissingData = TRUE,
  maxPOutliers = 1, 
  quickCor = 0,
  saveTOMFileBase = "All sample 12STOM",
  pearsonFallback = "individual",
  nThreads = 0,
  verbose = 0, indent = 0)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "All sample 12S-02-networkConstruction-auto.RData")

#####################################################################

lnames = load(file = "All sample-01-dataInput.RData");
lnames
lnames = load(file = "All sample 12S-02-networkConstruction-auto.RData");
lnames

# The merged module colors
mergedColors = net$colors;
# Eigengenes of the new merged modules:
mergedMEs = net$MEs;
table(mergedColors)

# Calculated intra- and total-connectivity of modules
kMEs <- signedKME(datExpr, mergedMEs)
ind <- sapply(mergedColors, function(color) {which(substr(colnames(kMEs), 4, 100)==color)})
modKMEs <- kMEs[cbind(1:nrow(kMEs), ind)]
write.table(cbind(ID=colnames(datExpr),dynamicColors=net$unmergedColors,mergedColors,modKMEs),
            "modules.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(cbind(sample=rownames(mergedMEs),mergedMEs),
            "mergedMEs.txt",sep="\t",quote=F,col.names=T,row.names=F)

#########################

rownames(mergedMEs) <- rownames(datExpr)

## plot mergedME tree
library(flashClust)
# Calculate dissimilarity of module eigengenes
mergedMEDiss = 1-cor(mergedMEs)

# Cluster module eigengenes
mergedMETree = flashClust(as.dist(mergedMEDiss), method = "average");
pdf("mergedME_tree.pdf", 12, 6)
plot(mergedMETree, main = "Clustering of module eigengenes",xlab = "", sub = "",cex=.7)
dev.off()

# reassignment 
modules0 <- read.table("modules.txt", as.is=T, header=T)
mergedMEs <- read.table("mergedMEs.txt", as.is=T, header=T, row.names=1)

ReassignMods <- function(datExpr, modColors0, cutoff=.7, max.iter=100) {
  for (i in 1:max.iter) {
    MEs <- moduleEigengenes(datExpr, colors = modColors0)$eigengenes
    kMEs <- signedKME(datExpr, MEs)
    modColors <- substr(colnames(kMEs), 4, 100)[apply(kMEs, 1, which.max)]
    modColors[apply(kMEs, 1, max) < cutoff] <- "grey"
    
    if(all(modColors0 == modColors)) {
      converge <- TRUE
      break;
    } else {
      modColors0 <- modColors
      converge <- FALSE
    }
  }
  return(list(modColors=modColors, MEs=MEs, kMEs=kMEs, converge=converge, num.iter=i))
}

modList <- ReassignMods(datExpr, modules0$mergedColors, cutoff=.7, max.iter=100)

modNames.color <- substr(colnames(modList$MEs), 3, 100)
modNames.num <- rep(0, length(modNames.color))
modNames.num[modNames.color=="grey"] <- 0
modNames.num[modNames.color!="grey"] <- 1:(length(modNames.color)-1)
newColors <- modNames.num[match(modList$modColors, modNames.color)]
newMEs <- moduleEigengenes(datExpr, colors = newColors)$eigengenes

newKMEs <- signedKME(datExpr, newMEs)
ind <- sapply(newColors, function(color) {which(substr(colnames(newKMEs), 4, 100)==color)})
newModKMEs <- newKMEs[cbind(1:nrow(newKMEs), ind)]
write.table(cbind(modules0, newColors, newModKMEs),"reassigned_modules.txt",
            sep="\t",quote=F,row.names=F,col.names=T)
write.table(cbind(sample=rownames(datExpr),newMEs),"reassignedMEs.txt",
            sep="\t",quote=F,col.names=T,row.names=F)

#######################
## plot reassignedME tree
# Calculate dissimilarity of module eigengenes
newMEDiss = 1-cor(newMEs)
# Cluster module eigengenes
newMETree = flashClust(as.dist(newMEDiss), method = "average");
pdf("reassignedME_tree.pdf", 12, 6)
plot(newMETree, main = "Clustering of module eigengenes",xlab = "", sub = "",cex=.7)
dev.off()

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(newColors)
# Plot the dendrogram and the module colors underneath
tiff(filename = "Cluster Dendrogram.tif", width = 5, height = 3, units = "in", res=500)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = newColors
moduleColors = labels2colors(newColors)
MEs = newMEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "12S-03-Reassignment.RData")

lnames = load(file = "12S-03-Reassignment.RData");
lnames

# Calculate dissimilarity of module eigengenes
newMEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
newMEDiss = 1-cor(newMEs)
# Cluster module eigengenes
newMETree = flashClust(as.dist(newMEDiss), method = "average");
pdf("reassignedME_tree_color.pdf", 12, 6)
plot(newMETree, main = "Clustering of module eigengenes",xlab = "", sub = "",cex=.7)
dev.off()

table(moduleColors)

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

tiff(filename = "Module-trait relationship.tif", width = 8, height = 12, units = "in", res=300)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 10, 2, 1));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


# Define variable weight containing the weight column of datTrait
MDD = as.data.frame(datTraits$Disease);
names(MDD) = "MDD"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")); # MM = kME
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, MDD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(MDD), sep="");
names(GSPvalue) = paste("p.GS.", names(MDD), sep="");

probes = names(datExpr)
# Create the starting data frame
geneInfo0 = data.frame(probes,moduleColor = moduleColors,geneTraitSignificance,GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, MDD, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.MDD));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_GS.csv")

# correlation |r|
sig_moduleTraitCor = signif(moduleTraitCor, 2)
sig_moduleTraitPvalue = signif(moduleTraitPvalue, 1)
write.csv(sig_moduleTraitCor, file="signif_moduleTraitCor.csv")
write.csv(sig_moduleTraitPvalue, file="signif_moduleTraitPvalue.csv")

#kME-----------------------------------------------------------------
options(stringsAsFactors = FALSE);
lnames = load(file = "All sample-01-dataInput.RData");
lnames
lnames = load(file = "12S-03-Reassignment.RData");
lnames

HuMDD = read.csv("Log2_FPKM_quan.csv", as.is = T, header = T);
# Take a quick look at what is in the data set:
dim(HuMDD);

names(HuMDD);

datExpr0 = as.data.frame(t(HuMDD[,-(1:4)]));
names(datExpr0) = HuMDD$Gene;
rownames(datExpr0) = names(HuMDD)[-(1:4)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

datExpr = datExpr0

kME = corAndPvalue(datExpr, MEs); # module eigengene(ME), correlation, P value 
kMEmat = cbind(kME$cor, kME$p);
MEnames = colnames(MEs); 
probes = names(datExpr)
ID=HuMDD$ID
rownames(kMEmat) = probes;
colnames(kMEmat) = c(spaste("kME", MEnames), spaste("p.kME", MEnames))

info = data.frame(Gene = probes, ID=ID, ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  kMEmat);
write.csv(info, file = "CombinedNetworkResults_Ensembl.csv",
          row.names = FALSE, quote = FALSE)
info = read.csv("CombinedNetworkResults_Ensembl.csv", as.is = T, header = T)

dim(info)
table(moduleLabels)

for (i in 0:50) {
  a = cbind(subset(info, ModuleLabel == i)[,c(1:4)],
            subset(info, ModuleLabel == i)[,c(spaste("kMEME",i), spaste("p.kMEME",i))])
  write.csv(subset(a,a[,4] >= 0.7),
            file = spaste("kME0.7Ensembl/kME0.7E_M",i,".csv"))
}

# Visualizing the network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
disease = as.data.frame(datTraits$Disease);
names(disease) = "MDD"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, disease))
# Plot the relationships among the eigengenes and the trait
tiff(filename = "ME relationship.tif", width = 8, height = 13, units = "in", res=300)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Module-module relationships", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

cor_eigengene = cor(MET,MET)
write.csv(cor_eigengene, file = "cor_eigengene.csv")

# Recalculate module eigengenes for eigengene raw data file
MEs = moduleEigengenes(datExpr, moduleLabels)$eigengenes
write.csv(MEs, file = "eigengenes.csv")


### Cytoscape
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

lnames = load(file = "All sample-01-dataInput.RData");
lnames
### need datExpr with gene symbols ###
lnames = load(file = "12S-03-Reassignment.RData");
lnames

TOM = TOMsimilarityFromExpr(datExpr, networkType = "signed", power = 12);

modules = "module color" #ex. blue, red...

# Select module probes  
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read for cell type modules (small size)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Cytoscape/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("Cytoscape/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.23,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])

Com = read.csv("CombinedNetworkResults_Ensembl.csv")

num = "module number" #ex. 1, 2...

# edge
edge = read.delim(paste("Cytoscape/CytoscapeInput-edges-",modules, ".txt", sep=""))
edge = edge[,1:3]
write.table(edge, paste("Cytoscape/edge_M", num, "_0.23.txt", sep=""),
            sep = "\t", row.names = F, quote = F)

# node
node = read.delim(paste("Cytoscape/CytoscapeInput-nodes-", modules, ".txt", sep=""))
node = node[, c(1,3)]
names(node) = c("Gene", "ModuleColor")

kME = subset (Com, ModuleLabel == num, select = c(Gene, ModuleLabel, kMEME39)) # M45

names(kME) = c("Gene", "ModuleLabel", "kME")
Cell = read.csv(paste("Cell type_H_M",num,".csv", sep=""))
node = merge(x=node, y=kME, by='Gene', all.x=T)
node = merge(x=node, y=Cell, by='Gene', all.x=T)
write.table(node, paste0("Cytoscape/node_M", num, "_0.23.txt"), sep = "\t", row.names = FALSE, quote = FALSE)



### Linear regression ------------------------------------------------
eigengene = read.csv("eigengenes.csv")
trait = read.csv("HuMDD_trait_LR.csv")
dim(eigengene)
dim(trait)

LRdata = merge(eigengene, trait)
dim(LRdata)
head(LRdata)
names(LRdata)

# - Diagnosis, Cause.of.Death, Toxicology, Smoker, Suicide
LR0 = lm(ME0 ~ Region + Disease + Age + Sex + PMI +
           pH + RIN + Suicide + Smoke + Medication + Toxicology,
         data = LRdata)
library(car)
ALR = Anova(LR0)
R = ALR$`Sum Sq`[-12]/sum(ALR$`Sum Sq`, na.rm = T)
P = ALR$`Pr(>F)`[-12]

for(i in 1:50) {
  a = lm(LRdata[,i+2] ~ Region + Disease + Age + Sex + PMI +
           pH + RIN + Suicide + Smoke + Medication + Toxicology, 
         data = LRdata)
  b = Anova(a)
  R = rbind(R, b$`Sum Sq`[-12]/sum(b$`Sum Sq`, na.rm = T))
  P = rbind(P, b$`Pr(>F)`[-12])
  C = cbind(R, P)
}
dim(C)
write.csv(C, "LR.csv")
