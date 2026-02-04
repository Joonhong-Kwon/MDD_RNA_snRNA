library(ggplot2)
library(tidyverse)
library(forcats)

### Combined Bar Plot
C <- read.csv("Total_All_Region_GO_BP_top5.csv")

pdf("All_BarPlot_GOBP.pdf", width = 10, height = 5)
Termorder <- C$Term[order(C$ratio)]
C$Term <- factor(C$Term, levels = rev(unique(C$Term)))
C$Category <- factor(C$Category, levels = unique(C$Category))
C$Color <- factor(C$Color, levels = c("dlPFC_top5", "CA1_top5", "DG_top5", "CBC_top5"))
ggplot(C, aes(x=Term, y=LogP, fill=Color)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c("#b05394",
                             "#cc924b",
                             "#8bad51",
                             "#3692ba")) +
  coord_flip() +
  theme_bw()+
  geom_hline(yintercept=1.3, linetype = 'dashed', color='red', size = 0.8) +
  facet_grid(. ~ Category, scales = "free_y", space = "free_y")+
  expand_limits(x=0) +
  labs(title="GO_BP", x=NULL, y="-Log10 (p-value)")+
  theme(plot.title = element_text(size=18, colour="black", face="bold", hjust=0.5),
        axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=11))
dev.off()



### Dot Plot
C <- read.csv ("Module_GO_total_top5.csv")
pdf("Module_GO_top5.pdf", width = 5.6, height = 6.5)
Termorder <- C$Term[order(C$Category, C$ratio)]
C$Term <- factor(C$Term, levels = Termorder)
sp2<-ggplot(C, aes(x=ratio, y=Term, colour=PValue, size=Count, shape=Category)) +
  geom_point() +
  scale_shape_manual(values=c(16, 17, 15, 18)) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")+
  expand_limits(x=0) +
  labs(title="Module", x="Gene ratio (%)", y=NULL, colour="p-value", size = "Count")+
  guides(size=guide_legend(order=1),
         shape=guide_legend(order=2),
         colour=guide_colourbar(order=3))
sp2+scale_color_gradient(low="#db1f5e", high="#7aa5de")+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.text.y = element_text(colour="black", size=10))+
  theme(plot.title = element_text(size=18, colour="black", face="bold", hjust=0.5))+
  scale_x_continuous(breaks = c(c(0:100)*10))
dev.off()



### Bar Plot
data <- read.csv ("Bar_Cluster_all.csv")

pdf("Bar_GO_Cluster.pdf", width = 6, height = 8.5)
data %>% 
  mutate(
    Cluster = fct_relevel(Cluster, "Cluster UP DEGs GO", "Cluster DOWN DEGs GO"),
    Term = fct_reorder(Term, Log10pvalue)
  ) %>% 
  ggplot(aes(x = Term, y = Log10pvalue, fill = Cluster)) +
  geom_col(alpha = 0.75, width = 0.85) +
  scale_fill_manual(values=c("#8E6CD4", "#C6B5E9")) +
  scale_y_continuous(expand = c(0, 0.1)) +
  coord_flip() +
  facet_grid(rows = vars(Cluster), scales = "free_y", switch = "y", space = "free_y") +
  labs(
    title = "Top5 GO_BP",
    y = "-Log10 (p-value)"
  ) +
  theme_minimal() +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    plot.title = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(angle = 270, face = "bold", size=13),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
    axis.title.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  geom_text(
    aes(Term, y = 0, label = Term),
    hjust = 0,
    nudge_x = 0.05,
    nudge_y = 0.03,
    colour = "black",
    size = 5
  )
dev.off()
