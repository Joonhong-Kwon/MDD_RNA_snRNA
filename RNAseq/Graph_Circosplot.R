library(circlize)
library(ComplexHeatmap)
library(grid)
#library(gridExtra)
#library(ggplot2)
#library(lemon)
#library(RColorBrewer)

### DEG circos plot
m2plot <- read.csv("Circos_DEG.csv", header=T)
MC <- t(m2plot[, 1])
Code <- t(m2plot[, 2])
rownames(Code) = "Module Color"
m2plot1 <- t(m2plot[, 3:10])
lable1=row.names(m2plot1);lable1
lable2=m2plot$Module;lable2
lable3=row.names(Code);lable3
colfunc<-colorRampPalette(c("yellow","red"))
cuts <- c(1.3, 2, 3, 4, 5, 6, 7, 8, 9, 10)
m2plot2 <- m2plot1
i <- 1
m2plot2[m2plot1 < cuts[i]] <- "#E6E6E6"
for(i in 1:9){
  m2plot2[m2plot1 >= cuts[i] & m2plot1 < cuts[i+1]] <- colfunc(10)[i]
}
i <- 10
m2plot2[m2plot1 >= cuts[i]] <- colfunc(10)[i]

m2plot3 <- rbind(Code, m2plot2)
lable4=row.names(m2plot3);lable4

nr = nrow(m2plot3);nr
nc = ncol(m2plot3);nc

col_mat = m2plot3
col_mat[,1]
col_code = m2plot3[1,]
col_code

col_fun = colorRamp2(c(0, 1.3, 2, 3, 4, 5, 6, 7, 8, 9,10),
                     c("#E6E6E6", "#FFFF00", "#FFE200", "#FFC600", "#FFAA00", "#FF8D00", "#FF7100", "#FF5500", "#FF3800", "#FF1C00", "#FF0000"))
col_leg=rev(c("#E6E6E6", "#FFFF00", "#FFE200", "#FFC600", "#FFAA00", "#FF8D00", "#FF7100", "#FF5500", "#FF3800", "#FF1C00", "#FF0000"))
value = rev(c("NS", "1.3", "2", "3", "4", "5", "6", "7", "8", "9", "≥10"))
col_value = structure(col_leg, names = unique(value))


pdf("Circos_DEG.pdf", width = 7, height = 7)
par(mar=c(0,0,0,0))
circos.clear();circos.par(canvas.xlim =c(-1.6,1.6),
                          canvas.ylim = c(-1.6,1.6),
                          cell.padding = c(0,0,0,0), 
                          gap.degree =60)
factors = "a"
circos.initialize(factors, xlim = c(0, ncol(m2plot3)))

circos.track(ylim = c(0, nr),bg.border = NA,track.height = 0.7,
             panel.fun = function(x, y) {
               for(i in 1:nr) {
                 circos.rect(xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
                             xright = 1:nc, ytop = rep(nr - i +1, nc),
                             border = "white",
                             col = col_mat[i,])
                 circos.text(x = nc,
                             y = 8.8 -i,
                             labels = lable4[i],
                             facing = "downward", niceFacing = TRUE,
                             cex = 0.6,
                             adj = c(-0.2, 0))
               }
             })

for(i in 1:nc){
  circos.text(x = i-0.4,
              y = 9.2,
              labels = lable2[i],
              facing = "clockwise", niceFacing = TRUE,
              cex = 0.8,adj = c(0, 0))
}

draw.sector(0, 60, rou1 = 0.94, rou2=0.915,  col="white", border="white", lwd = 1, lty = 1)
draw.sector(359.5, 60.5, rou1 = 0.757, rou2=0.757,  border="grey50", lwd = 1, lty = 1)
draw.sector(359.5, 60.5, rou1 = 0.601, rou2=0.601,  border="grey50", lwd = 1, lty = 1)
draw.sector(359.5, 60.5, rou1 = 0.446, rou2=0.446,  border="grey50", lwd = 1, lty = 1)
draw.sector(359.5, 60.5, rou1 = 0.291, rou2=0.291,  border="grey50", lwd = 1, lty = 1)
text(0, 0.01, "DEG\nenrichment", cex = 1)

lgd <- Legend(at = names(col_value), 
              legend_gp = gpar(fill = col_value),
              title_position = "lefttop-rot",title = "-Log10 (FDR)")
draw(lgd, x = unit(0.9, "npc"), y = unit(0.7, "npc"))
draw.sector(304.9, 311, rou1 = 1.65, rou2 = 0.27, clock.wise = FALSE, border=4, lwd=2)
dev.off()
