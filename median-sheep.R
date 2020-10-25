TFexpdataset<-read.table('sheep832.csv',sep=",",header=T)
p<-matrix(0,nrow=nrow(TFexpdata),ncol=50);q<-vector("numeric")
index<-as.vector(TFexpdataset[1,-1])
TFexpdata<-TFexpdataset[-1,-1]

for(k in 1:nrow(TFexpdata)){
  for(m in 1:50){
    q=c()
    for(j in 1:ncol(TFexpdata)){
      if(index[j]==m){
        q[j]<-TFexpdata[k,j]
      }
    }
    p[k,m]<-median(q, na.rm = T)
  }
}

write.csv(p,"median.csv")
