## The true coefficient matrix data is saved in BetaCoef_m20.csv
Beta<-as.matrix(read.table("D:/BetaCoef_m20.csv",sep=",",head=F))
m<-ncol(Beta)  ## number of time series
p<-5  ## number of lag
N<-600 ## size of data set
replication<-100
Data<-NULL
for(i in 1:replication)
{
Y<-rmvnorm(5,rep(0,m),diag(1,m))

for(k in 1:N){
  Y<-rbind(Y,rmvnorm(1,c(t(Y[nrow(Y):(nrow(Y)-(p-1)),]))%*%Beta[-c(1:m),],diag(1,m)))
  Y[nrow(Y),]<-Y[nrow(Y),]%*%solve(diag(1,m)-Beta[c(1:m),])
}
Data<-cbind(Data,Y[-c(1:p),])

}

