library(pheatmap)
library(grid)
library(RColorBrewer)

Peaks <- read.csv('openness_RUMEN.csv')
peak <- Peaks[,4:17]
peakLog <- log2(peak+1)

###################################
#heatmap of Chromatin accessibility
###################################
corp<-cor(peakLog,method = 'spearman')
annotation_row = data.frame(
  Tissue = factor(c(rep("Rumen", 14))) 
)
rownames(annotation_row) = colnames(Peaks)[4:17]
annotation_col = data.frame(
  Time = factor(c(rep("E60", 4),rep("D1", 3),rep("D7",3),rep("D28",3),rep('Y1',1)))
)
rownames(annotation_col) = colnames(Peaks)[4:17]
ann_colors = list(
  Time=c(E60="red",D1="orange",D7="cyan",D28="blue",Y1="purple")
)

pdf(file='heatmap of Chromatin accessibility.pdf')
p=pheatmap(corp,show_rownames=T,show_colnames=T,cellwidth = 10, cellheight = 10,pointsize=16,fontsize_row = 7,annotation_col = annotation_col, annotation_colors = ann_colors, border_color="black",fontsize_col = 7,gp = gpar(fontface = 'bold'))$gtable;
grid.draw(p)
dev.off()

###################################
#PCA of Chromatin accessibility
###################################
pca = prcomp(t(peakLog))
a<- summary(pca)
a$importance[,1:4]
colour = c(rep('red',4),rep('orange',3),rep('cyan',3),rep('blue',3), rep('purple',1))
pcha=c(rep(17,14))
pdf("PCA of Chromatin accessibility.pdf")
plot(pca$x[,c(1,2)], col=colour,pch=pcha, xlab = "Component1(68.459%)", ylab = "Component2(7.339%)", cex=2)
legend('bottomleft',legend = c("E60","D1","D7","D28","Y1"),fill = c("red","orange","cyan","blue","purple"), border =FALSE,pch = c(rep(NA,5)), bty="o")
title('PCA of Chromatin accessibility')
abline(h=0,v=0,lty=3)
dev.off()
