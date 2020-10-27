library(readxl)

al<-read.csv("Binding&Correlation.csv")
RNA<-read_excel("RNA_average.xlsx")
ATAC<-read_excel("ATAC_average.xlsx")

TF<-c()
TG<-c()
OPEN<-c()
for(i in 1:dim(al)[1]){
  s1 <- which(toupper(RNA$gene)==toupper(al$TFname[i]))
  s2 <- which(toupper(RNA$gene)==toupper(al$TTF[i]))
  s3 <- which(toupper(ATAC$openRSCNE)==toupper(al$RSCNE[i]))
  TF<-rbind(TF,RNA[s1,])
  TG<-rbind(TG,RNA[s2,])
  OPEN<-rbind(OPEN,ATAC[s3,])
}

####################
#integrate features
####################
ex_score<-TG
for(i in 1:dim(al)[1]){
  for(j in 2:dim(TF)[2]){
    ex_score[i,j]<-sqrt(TF[i,j]*TG[i,j])*(al$MotifScore[i])*(OPEN[i,j])*2^(al$cor.TTF.TF.S.[i])
  }
}
write.csv(cbind(al$RSCNE, al$TFname,ex_score),'all-score-single.csv',row.names = FALSE)

##########
#SUM TFTG
##########
library(readxl)
al<-read.table("all-score-single-tftg.txt", header = TRUE)
RSCNE<-read_excel("URSCNE.xlsx")
type1<-read.table('type1.conservativeS.txt')
type2<-read.table('type2-a.conservativeS.txt')
type <- rbind(type1, type2)

score<-c()
for(i in 1:dim(RSCNE)[1]){
  s1<-which(toupper(al$RSCNE)==toupper(RSCNE$URSCNE)[i])
  s2<-matrix(0,39,10)
  s3<-c()
  t <- which(type$V1==RSCNE$URSCNE[i])
  for(j in s1){
    s2<-s2+al[j,3:12]
    s3<-(cbind(s3, as.character(al$TF..TG[j])))
  }
  s4<-cbind(RSCNE$URSCNE[i], s2, type$V2[t] ,s3)
  score<-plyr::rbind.fill(score,s4)
}

write.csv(score, 'up-FI.csv',row.names = FALSE)

