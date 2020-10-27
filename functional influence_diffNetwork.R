library(readxl)

all<-read.table("diff-single-system.txt",header = TRUE)
RSCNE<-read_excel("DiffRSCNE.xlsx")
type1<-read.table('type1.conservativeS.txt')
type2<-read.table('type2-a.conservativeS.txt')
type <- rbind(type1, type2)

score<-c()
for(i in 1:dim(RSCNE)[1]){
  s1<-which(toupper(all$RSCNE)==toupper(RSCNE$DiffRSCNE)[i])
  s2<-matrix(0,dim(RSCNE)[1],10)
  s3<-c()
  t <- which(type$V1==RSCNE$DiffRSCNE[i])
  for(j in s1){
    s2<-s2+all[j,3:12]
    s3<-(cbind(s3, as.character(all$TF.tsys...TG.gsys.[j])))
  }
  DiffRSCNE<-RSCNE$DiffRSCNE[i]
  conservation_score<-type$V2[t] 
  TFTG<-s3
  s4<-as.data.frame(cbind(DiffRSCNE, s2, conservation_score,TFTG))
  score<-plyr::rbind.fill(score,s4)
}

write.csv(score, 'diff-FI.csv',row.names = FALSE)

