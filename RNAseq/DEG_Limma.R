library(limma)

expr <- read.csv("Log2_FPKM_quan_dlPFC.csv", row.names = "GeneID")
pheno <- read.csv("condition_dlPFC.csv")
colnames(pheno)
Disease <- pheno$Disease
Paired <- pheno$Paired
#RIN <- pheno$RIN
#Age <- pheno$Age
#Sex <- pheno$Sex
#pH <- pheno$pH
#Medi <- pheno$Medication
#Smoke <- pheno$Smoke
#Toxi <- pheno$Toxicology

design <- model.matrix(~0 + Disease + Paired)
colnames(design)
fit <- lmFit(expr,design)
cont <- makeContrasts(DiseaseMDD-DiseaseCon,levels=design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont, trend=TRUE)  # trend=TRUE 
res <- topTable(fit.cont,number=Inf)
head(res)

write.csv(res, file="Limma_trend_dlPFC_DEG.csv")

