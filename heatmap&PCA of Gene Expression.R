library(limma)
library(readxl)
library(pheatmap)
library(grid)
library(RColorBrewer)
RNA <- read.csv("RNA-seq_RUMEN.csv")
csif <- read_excel("sif.xlsx")
y<-as.matrix(log2(RNA[,2:15]+1))
batch <- csif$batch
design=model.matrix(~csif$cdv)
RNA_RBE <- removeBatchEffect(y, batch=batch,design=design)

############################
#heatmap of Gene Expression
############################
corp<-cor(RNA_RBE,method = 'spearman')
annotation_row = data.frame(
  Tissue = factor(c(rep("Rumen", 14))) 
)
rownames(annotation_row) = colnames(RNA_RBE)[1:14]
annotation_col = data.frame(
  Time = factor(c(rep("E60", 4),rep("D1", 3),rep("D7",3),rep("D28",3),rep('Y1',1)))
)
rownames(annotation_col) = colnames(RNA_RBE)[1:14]
ann_colors = list(
  Time=c(E60="red",D1="orange",D7="cyan",D28="blue",Y1="purple")
)

pdf(file='heatmap of Gene expression.pdf')
p=pheatmap(corp,show_rownames=T,show_colnames=T,cellwidth = 10, cellheight = 10,pointsize=16,fontsize_row = 7,annotation_col = annotation_col,  annotation_colors = ann_colors, border_color="black",fontsize_col = 7,gp = gpar(fontface = 'bold'))$gtable;
grid.draw(p)
dev.off()

#######################
#PCA of Gene Expression
#######################
rnaremove=RNA_RBE[-which(apply(RNA_RBE,1,var)==0),]
pca = prcomp(t(rnaremove),scale. = TRUE)
a<- summary(pca)
a$importance[,1:2]
colour = c(rep('red',4),rep('orange',3),rep('cyan',3),rep('blue',3), rep('purple',1))
pcha=c(rep(17,14))
pdf("PCA of Gene expression.pdf")
plot(pca$x[,1:2], col=colour,pch=pcha, xlab = "Component1(36.272%)", ylab = "Component2(17.110%)", cex=2)
legend('bottomright',legend = c("E60","D1","D7","D28","Y1"),fill = c("red","orange","cyan","blue","purple"), border =FALSE,pch = c(rep(NA,5)), bty="o")
title('PCA of Gene expression')
abline(h=0,v=0,lty=3)
dev.off()
