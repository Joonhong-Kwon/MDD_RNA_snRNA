library(lm.beta)
library(car)
library(dplyr)
library(readr)
library(factoextra)

# Generalized variance inflation factor
Trait <- read.csv("HuMDD_trait_PCA.csv")
head(Trait)
Trait$Region <- as.factor(Trait$Region)
Trait$Disease <- as.factor(Trait$Disease)
Trait$Sex <- as.factor(Trait$Sex)
Trait$Medication <- as.factor(Trait$Medication)
Trait$Smoke <- as.factor(Trait$Smoke)
Trait$Toxicology <- as.factor(Trait$Toxicology)

multiple.reg <- lm(Gene ~ Region + Disease + Age + Sex + RIN + pH + Medication + Smoke + Toxicology, data=Trait)
summary(multiple.reg)
vif(multiple.reg)


# Principal variance component analysis
dataExpr <- read.csv("Log2_FPKM_quan.csv")
Expr <- dataExpr[, -c(1:4)]
dim(Expr)
names(Expr)

trait <- read.csv("HuMDD_trait_PCA_region.csv")

pca.Expr <- prcomp(t(Expr))
pcaExpr <- pca.Expr$x

pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  abline(h=0.8, col="red")
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}

pcaCharts(pca.Expr)
summary(pca.Expr)


trait$Region = as.factor(trait$Region)
trait$Disease = as.factor(trait$Disease)
trait$Sex = as.factor(trait$Sex)
trait$Smoke = as.factor(trait$Smoke)
trait$Medication = as.factor(trait$Medication)
trait$Toxicology = as.factor(trait$Toxicology)

features <- c('Region', 'Disease', 'Age', 'Sex', 'pH', 'RIN',
              'Smoke', 'Medication', 'Toxicology')
features_id <- match(features, names(trait))
npc <- 9
rf <- sapply(1:9, function(i){
  lm.pc <- lm(pcaExpr[,1:npc]~trait[,features_id[i]])
  lapply(summary(lm.pc), function(x) x$r.squared) %>% unlist
})
rp <- sapply(1:9, function(i){
  lm.pc <- lm(pcaExpr[,1:npc]~trait[,features_id[i]])
  lapply(summary(lm.pc), function(x) pf(x$fstatistic["value"], x$fstatistic["numdf"], x$fstatistic["dendf"], lower=FALSE)) %>% unlist
})
colnames(rf) <- features 
rownames(rf) <- paste0('PC',1:npc)
write.csv(rf, file="PCA_trait_PC10_R2_VIF.csv")

colnames(rp) <- features 
rownames(rp) <- paste0('PC',1:npc)
write.csv(rp, file="PCA_trait_PC10_p_value_VIF.csv")
