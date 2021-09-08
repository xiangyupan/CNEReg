TFexpdata<-read.table('RNA-seq-TF.csv',sep=",",header=T)
opennessdata<-read.table('openness.csv',sep=",",header=T)
mapping<-read.table('19TTF-up-lm(cor0.6).csv',sep=",",header=T)
TFname<-unique(mapping[,1])

k=19       ###########different TTF#########
RSCNE<-unique(subset(mapping[,2],mapping$TTF == TFname[k]))
A<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);Ascale<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);
O<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);Oscale<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);

for(j in 1:length(RSCNE)){
  TFcomplex<-matrix(0,ncol=29);TFcomplexdata<-matrix(0,ncol=29)
  for(i in 1:nrow(mapping)){
    if(mapping[i,2] == RSCNE[j] && mapping[i,1]==TFname[k]){
        TFcomplex=as.numeric(subset(TFexpdata[,-1],TFexpdata$TF == as.character(mapping[i,3])))+TFcomplex
      TFcomplexdata<-(TFcomplex-min(TFcomplex))/(max(TFcomplex)-min(TFcomplex))
    }
  }
  A[j,]<-as.matrix(TFcomplex)
  Ascale[j,]<-as.matrix(TFcomplexdata)
  O[j,]<-as.numeric(subset(opennessdata[,-c(1,2)],opennessdata$openRSCNE == as.character(RSCNE[j])))
  Oscale[j,]<-(O[j,]-min(O[j,]))/(max(O[j,])-min(O[j,]))
}

TTF<-t(subset(TFexpdata[,-1],TFexpdata$TF == as.character(TFname[k])))
TTFscale<-(TTF-min(TTF))/(max(TTF)-min(TTF))
X=as.matrix(t(sqrt(Oscale*Ascale)))

openness=as.data.frame(X)

lm_TF = summary(lm (TTFscale ~ V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11, data = openness))  ########Choose Significant REs##########
write(RSCNE[which(lm_TF$coefficients[-1,4]<0.05)],"19TTF-up-RSCNE.csv")

