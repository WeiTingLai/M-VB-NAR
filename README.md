# A Modified VAR-deGARCH Model for Asynchronous Multivariate Financial Time Series via Variational Bayesian Inference
This code is the function of the Variational Bayesian Inference for the Modified VAR-deGARCH (M-VAR-deGARCH) Model.

All the simulation results can be reproduced by using this code.

# Data Generation
```
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


```
# Simulation:
```
## Read Data
## The true coefficient matrix data is saved in BetaCoef_m20.csv
Beta<-as.matrix(read.table("D:/BetaCoef_m20.csv",sep=",",head=F))

## All independently 100 replications data are saved in Datam20_rep100.csv
Data<-as.matrix(read.table("D:/Datam20_rep100.csv",sep=",",head=F))

## Set all the parameters:
p=5  ## number of lag
m=ncol(Beta) ## number of time series
N=300 ## size of training set
replication=100

## No group structure
segment <- c(rep(1,m))

coef_beta=phi <- matrix(0,nrow=replication,ncol=(((p+1)*m)*m))
coef_variance <-list()

## Set initial indicator matrix
adjma <- matrix(0,ncol=m,nrow=m)
adjma[c(1:(m/2)),c((m/2+1):m)]<-1

mspe_all <-as.list(NULL)
sigma_hat<-matrix(0,nrow=replication,ncol=m*m)

mspe=Y_hat<-NULL

## Here, we can conduct 100 replications and display the average results. If we want to show average results, the code replaces "for(k in 1:replication)"
## In the next example, we only display the 100th result.
start.time<-proc.time()
for(k in 100:100)
{
  
  Y_test<-as.matrix(Data[c((N+1):nrow(Data)),c(((k-1)*m+1):(k*m))])
  Y_train<-as.matrix(Data[c(1:N),c(((k-1)*m+1):(k*m))]) 
  
  message("Info: Replication:",k)

  ## When 'current' is set to T, the model considers the latest data and employs the M-VAR-deGARCH model. 
  fit<-VB_NAR(Y_train,segment=segment,lag=p,adj=adjma,maxit=10000,tol=1e-8,current=T)

  ###################################################################################
  ##      When 'current' is set to F, it employs the VAR-deGARCH model.             #
  #  fit<-VB_NAR(Y_train,segment=segment,lag=p,maxit=10000,tol=1e-8,current=F)      #
  ###################################################################################

  coef_beta[k,]<-fit$coef_beta
  Y_hat[[k]] <-  predict.VB_NAR(fit,Y_test,step_ahead = nrow(Data)-N,current=T)$y_hat

  mspe<-rbind(mspe,(Y_test-Y_hat[[k]])^2)

  betahat<-predict.VB_NAR(fit,Y_test,step_ahead = nrow(Data)-N,current=T)$bb

  mspe_all[[k]]<-mspe

## M-VAR-deGARCH
ACC<-TPR_FPR_TCR(fit,coef_beta,Beta1,100,current=T)

###################################################################################
##                  M-VAR-deGARCH (excluding $A_0$)                               #
#  number<-c(1:m)                                                                 #
#  for(i in 1:(m-1))                                                              #
#  {                                                                              #    
#    number<-c(number,c(((i*(m*(p+1)))+1):((i*(m*(p+1)))+m)))                     #
#  }                                                                              #
#  ACC1<-TPR_FPR_TCR(fit,coef_beta[,-number],Beta1[-c(1:m),],100,current=F)       #
###################################################################################

###################################################################################
##                                 VAR-deGARCH                                    #
#ACC<-TPR_FPR_TCR(fit,coef_beta,Beta1[-c(1:m),],100,current=F)                    #
###################################################################################
}

end.time<-proc.time()-start.time



## Average time
end.time/replication

## Show TPR, FPR, and average model size
ACC

###################################################################################
##                    M-VAR-deGARCH (excluding $A_0$)                             #
#    ACC1                                                                         #
###################################################################################
## Mean square prediction error
mean(apply(sapply(1:100,function(k)apply(mspe_all[[k]],2,mean)),2,mean))

## Heatmap of coefficient matrix ($A_0$, $A_1$, $A_3$, and $A_5$)

beta <-apply(coef_beta,2,mean)
beta<-matrix(beta,ncol=m)

library("RColorBrewer")
library("lattice")
library("gridExtra")
library("ggplot2")
library("grid")
brks <- c( seq(-0.6,0.6, length.out = 22))
rgb.palette <- colorRampPalette(c("red","white","blue"))

m=20
p=6
plt<-list()
#colnames(beta1) <- paste( c(1:50) , sep=" ")
for(l in 1:5)
{
  
  plt[[l]]<-levelplot(t(beta),at=brks,col.regions = rgb.palette , colorkey =F,  xlab = NULL, ylab = NULL,
                      panel=function(...) {
                        arg <- list(...)
                        panel.levelplot(...)
                        panel.abline(a = NULL, b = 1, h =seq(m+.5,m*p,by=m), lwd = .5,col="gray")
                        panel.abline(a = NULL, b = 1, h =seq(1.5,m*p,by=1), lwd = .5,col="gray",alpha=0.5)
                        panel.abline(a = NULL, b = 1, v =seq(1.5,m,by=1), lwd = .5,col="gray",alpha=0.5)
                      }, main= textGrob(paste("m20NG"," lag-",l-1,sep=""),gp = gpar(fontsize = 10)), ylim=c(m*l+.5, m*(l-1)+0.5), xaxt = "n"
                      ,scales=list(draw=FALSE ) )
  
}
plt[[2]]<-levelplot(t(beta),at=brks,col.regions = rgb.palette , colorkey =T,  xlab = NULL, ylab = NULL,
                    panel=function(...) {
                      arg <- list(...)
                      panel.levelplot(...)
                      panel.abline(a = NULL, b = 1, h =seq(m+.5,m*p,by=m), lwd = .5,col="gray")
                      panel.abline(a = NULL, b = 1, h =seq(1.5,m*p,by=1), lwd = .5,col="gray",alpha=0.5)
                      panel.abline(a = NULL, b = 1, v =seq(1.5,m,by=1), lwd = .5,col="gray",alpha=0.5)
                    }, main= textGrob(paste("m20NG"," lag-",1,sep=""),gp = gpar(fontsize = 10)), ylim=c(m*2+.5, m*(2-1)+0.5), xaxt = "n"
                    ,scales=list(draw=FALSE ) )


plt[[6]]<-levelplot(t(beta),at=brks,col.regions = rgb.palette , colorkey =T,  xlab = NULL, ylab = NULL,
                    panel=function(...) {
                      arg <- list(...)
                      panel.levelplot(...)
                      panel.abline(a = NULL, b = 1, h =seq(m+.5,m*p,by=m), lwd = .5,col="gray")
                      panel.abline(a = NULL, b = 1, h =seq(1.5,m*p,by=1), lwd = .5,col="gray",alpha=0.5)
                      panel.abline(a = NULL, b = 1, v =seq(1.5,m,by=1), lwd = .5,col="gray",alpha=0.5)
                    }, main= textGrob(paste("m20NG"," lag-",5,sep=""),gp = gpar(fontsize = 10)), ylim=c(m*6+.5, m*(6-1)+0.5), xaxt = "n"
                    ,scales=list(draw=FALSE ) )



grid.arrange(plt[[1]],plt[[2]],plt[[4]],plt[[6]],ncol=2, top = textGrob("The estimated coefficient matrix", gp = gpar(fontsize = 13, fontface = 'bold')))



```
