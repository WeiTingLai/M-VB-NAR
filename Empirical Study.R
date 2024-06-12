
VARData  <- as.matrix( read.csv( "D:/Empirical.csv", head = F, row.names = 1 ) )
adjma    <- as.matrix( read.csv( "D:/opentime.csv", head = T, row.names = 1 ) )

rownames(VARData)[1352]  #2015/04/01
rownames(VARData)[1608]  #2016/03/31
rownames(VARData)[1868]  #2017/03/31
rownames(VARData)[2123]  #2018/03/29
rownames(VARData)[2380]  #2019/03/29
rownames(VARData)[2638]  #2020/03/31
rownames(VARData)[2890]  #2021/03/24
###############################   the example for M1 training model with pi=0.001.
p = 5;m = ncol(adjma);N = 2123 # 2123 2380 2638

segment   <- c(rep(1,m)) 
mspe_all  <- as.list(NULL)
sigma_hat <- matrix(0,nrow=m*p,ncol=m*m)
phi       <- matrix(0,nrow=m*p,ncol=length(segment)*m*(p+1)+(p+1)*m)


start.time<-proc.time()
mspe = Y_hat = mspe_VAR = Y_hat_VAR <- NULL

Y_train<-as.matrix(VARData1[c(1352:N),])  # 1352 1609 1869
Y_test<-as.matrix(VARData1[c((N+1):2890),])

########################  M-VAR-DeGARCH  #########################

fit       <- VB_NAR(Y_train,segment=segment,adj = adjma,alpha=c(0.001,0.001) ,lag=p,maxit=10000,tol=1e-8,current=T)
Y_hat     <- predict.VB_NAR(fit,Y_test,step_ahead = nrow(Y_test),current=T,rolling = T, rolling_num= 250)$y_hat
betahat   <- predict.VB_NAR(fit,Y_test,step_ahead =nrow(Y_test),current=T,rolling = T, rolling_num= 250)$bb
mspe      <- rbind(mspe,(Y_test-Y_hat)^2)
mae       <- abs((Y_test-Y_hat))
coef_beta <- fit$coef_beta

###########################  Table 3  #############################
# Here shows the AIC in M1 with pi=0.001.

AIC       <- (nrow(Y_train)*m)*log(sum((fit$Y-fit$X%*%matrix(fit$coef_beta,ncol=ncol(Y_train)))^2 )/(nrow(Y_train)*ncol(Y_train)))+2*ncol(Y_train)^2*p

###########################  VAR-DeGARCH  ##########################

fit_VAR       <-  VB_NAR(Y_train,segment=segment,adj = adjma,alpha=c(0.001,0.001) ,lag=p,maxit=10000,tol=1e-8,current=F)
Y_hat_VAR     <-  predict.VB_NAR(fit_VAR,Y_test,step_ahead = nrow(Y_test),current=F,rolling = T, rolling_num= 250)$y_hat
betahat_VAR   <-  predict.VB_NAR(fit_VAR,Y_test,step_ahead =nrow(Y_test),current=F,rolling = T, rolling_num= 250)$bb
mspe_VAR      <-  rbind(mspe_VAR,(Y_test-Y_hat_VAR)^2)
mae_VAR       <-  abs((Y_test-Y_hat_VAR))
coef_beta_VAR <- fit_VAR$coef_beta

############################  ARMA-GARCH  ###########################

window            <- 250
Y_hat_AR          <- matrix(0,nrow=nrow(Y_test),ncol=m)
for(i in 1:m)
{
  AR_model        <- stats::arima(Y_train[,i],order=c(p,0,0),include.mean = T,method="ML")
  Y_hat_AR[1,i]   <- predict(AR_model,1)$pred[1]
  for(j in 2:(nrow(Y_test)-1))
  {
    if(j<=window)
    {
      AR_model    <- stats::arima(c(Y_train[c((nrow(Y_train)-window+j):nrow(Y_train)),i],Y_test[c(0:(j-1)),i]),order=c(p,0,0),include.mean = T,method="ML")  
    }else{
      AR_model    <- stats::arima(Y_test[c(0:(j-1)),i],order=c(p,0,0),include.mean = T,method="ML")
    }
    
    Y_hat_AR[j,i] <- predict(AR_model,1)$pred[1]
  }
}
mspe_AR           <- (Y_test-Y_hat_AR)^2
mae_AR            <- abs((Y_test-Y_hat_AR))

end.time<-proc.time()-start.time

#############################  Table 4  ###############################
##   M1-P1
apply(  mspe     [c(1:257), ], 2 , mean )[11:15]    # 1:257  1:258  1:nrow(Y_test) 
apply(  mspe_VAR [c(1:257), ], 2 , mean )[11:15]
apply(  mspe_AR  [c(1:257), ], 2 , mean )[11:15]

apply(  mspe     [c(1:257), ], 2 , sd )[11:15]
apply(  mspe_VAR [c(1:257), ], 2 , sd )[11:15]
apply(  mspe_AR  [c(1:257), ], 2 , sd )[11:15]
##   M1-P2
apply(  mspe     [c(258:516), ], 2 , mean )[11:15]    # 258：516  259:nrow(Y_test)
apply(  mspe_VAR [c(258:516), ], 2 , mean )[11:15]
apply(  mspe_AR  [c(258:516), ], 2 , mean )[11:15]

apply(  mspe     [c(258:516), ], 2 , sd )[11:15] 
apply(  mspe_VAR [c(258:516), ], 2 , sd )[11:15]
apply(  mspe_AR  [c(258:516), ], 2 , sd )[11:15]

##   M1-P3
apply(  mspe     [c(517:nrow(Y_test)), ], 2 , mean )[11:15]
apply(  mspe_VAR [c(517:nrow(Y_test)), ], 2 , mean )[11:15]
apply(  mspe_AR  [c(517:nrow(Y_test)), ], 2 , mean )[11:15]

apply(  mspe     [c(517:nrow(Y_test)), ], 2 , sd )[11:15] 
apply(  mspe_VAR [c(517:nrow(Y_test)), ], 2 , sd )[11:15]
apply(  mspe_AR  [c(517:nrow(Y_test)), ], 2 , sd )[11:15]

###########################  Figure 5  ##################################

library("RColorBrewer")
library("gridExtra")
library("lattice")
library("ggplot2")
library("grid")

brks        <- c( seq(-1,1, length.out = 22))
rgb.palette <- colorRampPalette(c("red","white","blue"))

##   M1-P1
beta        <- apply( betahat[c(1:256), ], 2 , mean )
beta        <- matrix( beta , ncol = m )

plt[[1]] <- levelplot(t(beta), at = brks, col.regions = rgb.palette, colorkey = T,  xlab = NULL, ylab = NULL,
                        panel=function(...) {
                          arg <- list(...)
                          panel.levelplot(...)
                          panel.abline(a = NULL, b = 1, h = seq( m+.5, m*p, by = m ), lwd = .5, col = "gray" )
                          panel.abline(a = NULL, b = 1, h = seq( 1.5, m*(p+1), by = 1 ), lwd = .5, col = "gray", alpha = 0.5 )
                          panel.abline(a = NULL, b = 1, v = seq( 1.5, m, by = 1 ), lwd = .5, col = "gray", alpha = 0.5 )
                        }, main = textGrob( expression(A[0]), gp = gpar(fontsize = 10) ), ylim = c( m*l+.5, m*(l-1)+0.5 ), xaxt = "n"
                        , scales = list( draw = FALSE ) )
plt[[1]]

##   M1-P2
beta        <- apply( betahat[ c(257:515), ], 2 , mean )
beta        <- matrix( beta , ncol = m )

plt[[1]] <- levelplot(t(beta), at = brks, col.regions = rgb.palette, colorkey = T,  xlab = NULL, ylab = NULL,
                        panel=function(...) {
                          arg <- list(...)
                          panel.levelplot(...)
                          panel.abline(a = NULL, b = 1, h = seq( m+.5, m*p, by = m ), lwd = .5, col = "gray" )
                          panel.abline(a = NULL, b = 1, h = seq( 1.5, m*(p+1), by = 1 ), lwd = .5,col="gray",alpha = 0.5 )
                          panel.abline(a = NULL, b = 1, v = seq( 1.5, m,by=1 ), lwd = .5, col = "gray", alpha = 0.5 )
                        }, main= textGrob( expression(A[0]), gp = gpar(fontsize = 10)), ylim = c( m*l+.5, m*(l-1)+0.5 ), xaxt = "n"
                        , scales = list( draw = FALSE ) )
plt[[1]]

##   M1-P3
beta        <- apply( betahat[ c(516:nrow(betahat)), ], 2 , mean )
beta        <- matrix( beta , ncol = m )

plt[[1]] <- levelplot(t(beta), at = brks, col.regions = rgb.palette, colorkey = T,  xlab = NULL, ylab = NULL,
                        panel=function(...) {
                          arg <- list(...)
                          panel.levelplot(...)
                          panel.abline(a = NULL, b = 1, h =seq( m+.5, m*p, by = m ), lwd = .5,col = "gray" )
                          panel.abline(a = NULL, b = 1, h =seq( 1.5, m*(p+1), by = 1 ), lwd = .5,col = "gray", alpha = 0.5 )
                          panel.abline(a = NULL, b = 1, v =seq( 1.5, m, by = 1 ), lwd = .5, col = "gray", alpha = 0.5 )
                        }, main= textGrob(expression(A[0]),gp = gpar(fontsize = 10)), ylim=c(m*l+.5, m*(l-1)+0.5), xaxt = "n"
                        ,scales=list(draw=FALSE ) )
plt[[1]]

