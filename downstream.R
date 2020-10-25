openness<-read.table('openness.csv',sep=",",header=T)
Filtered_RNAseq<-read.table('RNAseq.csv',sep=",",header=T)
label<-read.table('17TTF-downstream.csv',sep=",",header=T)
dim(label)
corr<-vector();p_value<-vector();q_value<-vector()

for(i in 1:4527){
  Odata<-subset(openness[,-c(1,2)],openness$region == as.character(label[i,1]))
  O<-as.numeric((Odata-min(Odata))/(max(Odata)-min(Odata)))
  Edata<-subset(Filtered_RNAseq[,-1],Filtered_RNAseq$gene == as.character(label[i,5]))
  E<-(Edata-min(Edata))/(max(Edata)-min(Edata))
  TFdata<-subset(Filtered_RNAseq[,-1],Filtered_RNAseq$gene == as.character(label[i,3]))
  TF<-as.numeric(((TFdata-min(TFdata))/(max(TFdata)-min(TFdata))))
  corr[i] <- cor(as.numeric(sqrt(O*TF)),as.numeric(E),method="spearman")
  p_value[i] <- cor.test(as.numeric(sqrt(O*TF)),as.numeric(E),method="spearman")$p.value
}

q_value <- p.adjust(p_value,method = "fdr")
write.csv(corr,"corr.csv")
write.csv(p_value,"pvalue.csv")
write.csv(q_value,"qvalue(fdr).csv")

