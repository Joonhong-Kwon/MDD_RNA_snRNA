library(ggplot2)
library(forcats)
library(muted)

# M9
C <- read.csv ("M9_hub_genes_DEG_top20.csv")
pdf("M9_hubgene_DEG_top20_FC1.3.pdf", width = 5.5, height = 7.5)
C$Gene <- factor(C$Gene, levels = c("GPSM3", "LY86", "C3", "C1QC", "CYBA",
                                    "UCP2", "UBA7", "TYROBP", "CD4", "LAIR1",
                                    "SPI1", "CD74", "CXCL16", "LAPTM5","ALOX5AP",
                                    "CD37", "HLA-DMA", "CSF1R", "C1QB", "SELPLG"))
C$Region <- factor(C$Region, levels = c("dlPFC", "CA1", "DG", "CBC"))
sp2<-ggplot(C, aes(x=Region, y=Gene, size=-log10(pvalue))) +
  scale_size_continuous(limits = c(0.0001, 4), range = c(1,13), breaks=c(0.1, 0.5, 1, 2, 4)) + 
  geom_point(aes(fill=Log2FC), colour="black", pch=21) +
  annotate("text", x=1, y=c(7.95, 8.95), label="*")+
  annotate("text", x=2, y=c(13.95, 14.95), label="*")+
  annotate("text", x=3, y=c(3.95), label="*")+
  annotate("text", x=4, y=c(0.95, 2.95, 3.95, 6.95, 7.95, 9.95, 10.95, 11.95, 13.95, 14.95, 15.95, 17.95, 18.95), label="*")+
  labs(title="M9 Hub genes", x=NULL, y=NULL, fill="Log2 (FC)", size = "-Log10 (p-value)")+
  guides(colour=guide_colourbar(order=1),
         size=guide_legend(order=2))
sp2+scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-1.5, 1.5), breaks=c(-1.5, 0, 1.5))+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black", size=15))+
  theme(axis.text.y = element_text(colour="black", size=12, face="italic"))+
  theme(plot.title = element_text(size=19, colour="black", face="bold", hjust=0.5))+
  theme(panel.border = element_blank())

dev.off()
