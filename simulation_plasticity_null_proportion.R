
num=10000; fp=0;
for (gene in 1:num){
  a<-rnorm(8,mean=0,sd=0.1)
  b<-rnorm(8,mean=0,sd=0.1)
  dif=c()
  for (i in 1:8){
    dif=c(dif,a[i]-b[i])
  }
  avg_dif=mean(dif)
  n=0; p=0;
  for (j in 1:8){
    if(dif[j]<0){
      n=n+1;
    }
    if(dif[j]>0){
      p=p+1;
    }
  }
  if ((n >6)&(avg_dif <0)){
    fp=fp+1;
  }
  if ((p >6)&(avg_dif >0)){
    fp=fp+1;
  }
}
fp/num
