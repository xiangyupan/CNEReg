JSDdata<-read.table('18TTFs-JSD.csv',sep=",",header=T)
mediandata<-read.table('18TTFs-median.csv',sep=",",header=T)
JSD<-1./JSDdata[,-1]
median<-(mediandata[,-1])^(1/3)

library(factoextra)

pdf("Phylogeny of 50 tissues.pdf",width=6,height=5)
result <- get_dist(t(sqrt(JSD*median)),stand = TRUE,method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")

fviz_dend(result_hc, k = 11, cex = 0.8, 
          k_colors = c("#A6D854","#FDAE61","#252525","#807DBA","#A6D854","#252525","#C2A5CF",
                       "#C2A5CF","#FA9FB5","#FA9FB5","#252525"),
          color_labels_by_k = TRUE, 
          rect = TRUE,main="Phylogeny of 50 tissues by 18TTF")
dev.off()