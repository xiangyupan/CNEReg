library(readxl)

RNA<-read_excel("RNA_average.xlsx")
ATAC<-read_excel("ATAC_average.xlsx")
BC_down<-read.csv("Binding&Correlation_down.csv")

TF<-c()
TG<-c()
OPEN<-c()
for(i in 1:dim(BC_down)[1]){
  s1 <- which(toupper(RNA$gene)==toupper(BC_down$TF[i]))
  s2 <- which(toupper(RNA$gene)==toupper(BC_down$TG[i]))
  s3 <- which(toupper(ATAC$openRSCNE)==toupper(BC_down$RSCNE[i]))
  TF<-rbind(TF,RNA[s1,])
  TG<-rbind(TG,RNA[s2,])
  OPEN<-rbind(OPEN,ATAC[s3,])
}

####################
#integrate features
####################
ex_score<-TG
for(i in 1:dim(BC_down)[1]){
  for(j in 2:dim(TF)[2]){
    ex_score[i,j]<-sqrt(TF[i,j]*TG[i,j])*(BC_down$MotifScore[i])*(OPEN[i,j])*2^(BC_down$correlation[i])
  }
}
aRSCNE<-BC_down$RSCNE
TF..TG<-paste(BC_down$TF,ex_score$gene,sep="->")
all<-as.data.frame(cbind(aRSCNE, TF..TG, ex_score[,-1]))

##########
#SUM TFTG
##########
RSCNE<-read_excel("DRSCNE.xlsx")
type1<-read.table('type1.conservativeS.txt')
type2<-read.table('type2-a.conservativeS.txt')
type <- rbind(type1, type2)

score<-c()
for(i in 1:dim(RSCNE)[1]){
  s1<-which(toupper(all$aRSCNE)==toupper(RSCNE$DRSCNE)[i])
  s2<-matrix(0,dim(RSCNE)[1],10)
  s3<-c()
  t <- which(type$V1==RSCNE$DRSCNE[i])
  for(j in s1){
    s2<-s2+all[j,3:12]
    s3<-(cbind(s3, as.character(all$TF..TG[j])))
  }
  DRSCNE<-RSCNE$DRSCNE[i]
  conservation_score<-type$V2[t] 
  TFTG<-s3
  s4<-as.data.frame(cbind(DRSCNE, s2, conservation_score,TFTG))
  score<-plyr::rbind.fill(score,s4)
}

write.csv(score, 'down-FI.csv',row.names = FALSE)

