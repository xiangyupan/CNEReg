########################JSD-sheep###########################
TFexpdata<-read.table('sheep832.csv',sep=",",header=T)
library(philentropy)

p<-matrix(0,nrow=50,ncol=ncol(TFexpdata));P<-matrix(0,nrow=50,ncol=ncol(TFexpdata))
index<-as.vector(TFexpdata[1,-1])
for(i in 1:50){
  sum=0
  for(l in 1:ncol(TFexpdata) ){
    if(index[l]==i){
      p[i,l]<-1
      sum<-sum+p[i,l]
    }
    P[i,]<-p[i,]/sum
  }
}

q<-matrix(nrow=nrow(TFexpdata),ncol=ncol(TFexpdata));Q<-matrix(nrow=nrow(TFexpdata),ncol=ncol(TFexpdata));
x<-matrix(nrow=2,ncol=ncol(TFexpdata))
JSD<-matrix(nrow=nrow(TFexpdata),ncol=50)
for(m in 1:50){
  for(k in 1:nrow(TFexpdata)){
    q[k,]<-as.numeric(as.character(as.data.frame(TFexpdata[k+1,-1])))
    sum=0
    for(j in 1:ncol(TFexpdata)){
      sum<-sum+q[k,j]
    }
    Q[k,]<-q[k,]/sum;
    x<- rbind(P[m,],Q[k,]);
    JSD[k,m]<-distance(x, method = "jensen-shannon");
  }
}
write.csv(JSD,"JSD_sheep.csv")

########################JSD-median###########################
TFexpdataset<-read.table('sheep832.csv',sep=",",header=T)
a<-matrix(0,nrow=nrow(TFexpdata),ncol=50);b<-vector("numeric")
index<-as.vector(TFexpdataset[1,-1])
TFexpdata<-TFexpdataset[-1,-1]

for(k in 1:nrow(TFexpdata)){
  for(m in 1:50){
    b=c()
    for(j in 1:ncol(TFexpdata)){
      if(index[j]==m){
        b[j]<-TFexpdata[k,j]
      }
    }
    a[k,m]<-median(b, na.rm = T)
  }
}
write.csv(a,"median_sheep.csv")

##########################JMscore###########################
JMscore<-matrix(nrow=nrow(TFexpdata),ncol=50)
JSDdata<-1./JSD
mediandata<-a^(1/3)
JMscore<-sqrt(JSDdata*mediandata)

write.csv(JMscore,"JMscore_sheep.csv")