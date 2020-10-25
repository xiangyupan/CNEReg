TFexpdata<-read.table('RNA-seq-TF.csv',sep=",",header=T)
opennessdata<-read.table('openness.csv',sep=",",header=T)
TTFbinding<-read.table('19TTF-upbinding.csv',sep=",",header=T)

TFname<-unique(TTFbinding[,1])
corrsm<-vector()

k=19       ###########different TTF#########
RSCNE<-unique(subset(TTFbinding[,2],TTFbinding$TTF == TFname[k]))
A<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);Ascale<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);
O<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);Oscale<-matrix(nrow=length(RSCNE),ncol=29,byrow=T);

for(j in 1:length(RSCNE)){
    TFcomplex<-matrix(0,ncol=29);TFcomplexdata<-matrix(0,ncol=29)
  for(i in 1:nrow(TTFbinding)){
    if(TTFbinding[i,2] == RSCNE[j] && TTFbinding[i,1]==TFname[k]){
        corrsm[i]<-cor(as.numeric(subset(TFexpdata[,-1],TFexpdata$TF == TFname[k])),as.numeric(subset(TFexpdata[,-1],TFexpdata$TF == as.character(TTFbinding[i,3]))),method="spearman")
       if(corrsm[i]>0.6){
          #print(cbind(TTFbinding[i,c(2,3)],corrsm[i]))
          TFcomplex=as.numeric(subset(TFexpdata[,-1],TFexpdata$TF == as.character(TTFbinding[i,3])))+TFcomplex
        }
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
library(glmnet)
Y=as.matrix(TTFscale)
X=as.matrix(t(sqrt(Oscale*Ascale)))
X
colnames(X)<-RSCNE
set.seed(1)
cv.fit=cv.glmnet(X, Y, grouped=FALSE)
coef(cv.fit)
library(coefplot)

rownames(extract.coef(cv.fit))[-1]    ##########Features after LASSO screening#########


#######################subsets after LASSO #############################
TTFbinding<-read.table('19TTF-upbinding.csv',sep=",",header=T)
lable<-read.table('19TTF-RSCNE(LASSO).csv',sep=",",header=T)
E<-matrix(nrow=nrow(TTFbinding),ncol=ncol(TTFbinding),byrow=T)
for (n in 1:nrow(TTFbinding)) {
  for(m in 1:nrow(lable)){
    if(TTFbinding[n,2] == lable[m,2])
      E[n,]<-as.matrix(TTFbinding[n,])
    else next
  }
}
write.csv(na.omit(E),"19TTF-up-lm.csv")


