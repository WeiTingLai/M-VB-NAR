
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
mae_AR            <- abs((Y_test-Y_hat_VAR))

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

###########################  Figure 6  ##################################

library("reshape")
library("dplyr") 

## Here, we adopt the 1% percentile of the in-sample returns as the threshold of the occurrence of high-risk events under the previous rolling window framework  

high_risk1   <- round( apply( Y_train[ c(nrow(Y_train):(nrow(Y_train)-249)) , ] , 2 , function(a) quantile( a , 0.01 ) ) , 2 )
for(i in 1:(nrow(Y_test)-1))
{
  Y_hirisk   <- rbind(Y_train , Y_test[c(1:i) , ] )
  high_risk1 <- rbind(high_risk1 , 
                    round( apply( Y_hirisk[ c(nrow(Y_hirisk):(nrow(Y_hirisk)-249)) , ] , 2 , function(a) quantile( a , 0.01 ) ) , 2 ) )
}


# Since each method is computationally expensive, we preserve the results for JTOPI, RTSI, FCHI, FTSE, and DAX across all different training models (M1, M2, and M3) for each method (VAR-deGARCH, and M-VAR-deGARCH).  
## To show the result for different trained models M, show the file path   M1:D:/Yhat/2019/.   M2:D:/Yhat/2020/.    M3: D:/Yhat/2021/.
## The example shows M1-P1, M1-P2, and M1-P3 in Figure 6
                                 
## M1-JTOPI ##
                                
Unfixed   <- read.csv("D:/Yhat/2019/南非/fixed_real return compare.csv")    # M-VAR-deGARCH
VB_NAR    <- read.csv("D:/Yhat/2019/南非/VBNAR return compare.csv")         # VAR-deGARCH
Unit_time <- read.csv("D:/Yhat/2019/南非/unit_real return compare.csv")     # ARMA-GARCH

new_dataS.Africa <- data.frame( Date = as.Date(Unfixed[,1]) , True = Unfixed[,3] , ARMA_GARCH = Unit_time[,4] , VAR_GARCH = VB_NAR[,4] , M_VAR_GARCH = Unfixed[,4] )
                                
## M1-RTSI ##
                                
Unfixed   <- read.csv("D:/Yhat/2019/俄羅斯/fixed_real return compare.csv")   # M-VAR-deGARCH
VB_NAR    <- read.csv("D:/Yhat/2019/俄羅斯/VBNAR return compare.csv")        # VAR-deGARCH
Unit_time <- read.csv("D:/Yhat/2019/俄羅斯/unit_real return compare.csv")    # ARMA-GARCH

new_dataRUS <- data.frame( Date = as.Date(Unfixed[,1]) , True = Unfixed[,3] , ARMA_GARCH = Unit_time[,4] , VAR_GARCH = VB_NAR[,4] , M_VAR_GARCH = Unfixed[,4] )

## M1-FCHI ##
                                
Unfixed   <- read.csv("D:/Yhat/2019/法國/fixed_real return compare.csv")   # M-VAR-deGARCH
VB_NAR    <- read.csv("D:/Yhat/2019/法國/VBNAR return compare.csv")        # VAR-deGARCH
Unit_time <- read.csv("D:/Yhat/2019/法國/unit_real return compare.csv")    # ARMA-GARCH

new_dataFR < -data.frame( Date = as.Date(Unfixed[,1]) , True = Unfixed[,3] , ARMA_GARCH = Unit_time[,4] , VAR_GARCH = VB_NAR[,4] , M_VAR_GARCH = Unfixed[,4] )

## M1-FTSE ##
                                
Unfixed   <- read.csv("D:/Yhat/2019/英國/fixed_real return compare.csv")   # M-VAR-deGARCH
VB_NAR    <- read.csv("D:/Yhat/2019/英國/VBNAR return compare.csv")        # VAR-deGARCH
Unit_time <- read.csv("D:/Yhat/2019/英國/unit_real return compare.csv")    # ARMA-GARCH

new_dataUK <- data.frame( Date = as.Date(Unfixed[,1]) , True = Unfixed[,3] , ARMA_GARCH = Unit_time[,4] , VAR_GARCH = VB_NAR[,4] , M_VAR_GARCH = Unfixed[,4] )

## M1-DAX ##
                                
Unfixed   <- read.csv("D:/Yhat/2019/德國/fixed_real return compare.csv")   # M-VAR-deGARCH
VB_NAR    <- read.csv("D:/Yhat/2019/德國/VBNAR return compare.csv")        # VAR-deGARCH
Unit_time <- read.csv("D:/Yhat/2019/德國/unit_real return compare.csv")    # ARMA-GARCH

new_dataBRD <- data.frame( Date = as.Date(Unfixed[,1]) , True = Unfixed[,3] , ARMA_GARCH = Unit_time[,4] , VAR_GARCH = VB_NAR[,4] , M_VAR_GARCH = Unfixed[,4] )


## To show the 1% percentile of the high-risk events


new_data_clearS.Africa  <- new_dataS.Africa[ which( new_dataS.Africa[,2] <= high_risk1[,11] ) , ]
new_data_clearRUS       <- new_dataRUS[ which( new_dataRUS[,2] <= high_risk1[,12] ) , ]
new_data_clearFR        <- new_dataFR[ which( new_dataFR[,2]   <= high_risk1[,13] ) , ]
new_data_clearUK        <- new_dataUK[ which( new_dataUK[,2]   <= high_risk1[,14] ) , ]
new_data_clearBRD       <- new_dataBRD[ which( new_dataBRD[,2] <= high_risk1[,15] ) , ]

## caculate the ratio R in equation(11)

## JTOPI                                
attach(new_data_clearS.Africa)
new_data_clearS.Africa1 <- abs( data.frame( ARMA_GARCH = ARMA_GARCH-True , VAR_GARCH = VAR_GARCH-True , M_VAR_GARCH = M_VAR_GARCH-True ) )
attach(new_data_clearS.Africa1)
new_data_clearS.Africa1 <- data.frame( ARMA_GARCH = ARMA_GARCH/M_VAR_GARCH , VAR_GARCH = VAR_GARCH/M_VAR_GARCH )
S.Africa <- new_data_clearS.Africa1[c(1:3) , ]  # M1-P1 
#S.Africa <- new_data_clearS.Africa1[c(4:9) , ] # M1-P2

S.Africa1 <- data.frame( ID = as.character( new_data_clearS.Africa$Date )[ c(1:3) ] , S.Africa[ , 2 ] )
newdataS.Africa <- data.frame( Index = "JTOPI" , melt( S.Africa1 , id = c("ID") ) )
#newdataS.Africa <-data.frame( Index = "JTOPI" , melt( data.frame( ID = "2020-06-11" , ARMA_GARCH = 0 , VAR_GARCH = 0 ) , id = c("ID") ) )   # M1-P3
                                
## RTSI
attach(new_data_clearRUS)
new_data_clearRUS1 <- abs( data.frame( ARMA_GARCH = ARMA_GARCH-True , VAR_GARCH = VAR_GARCH-True , M_VAR_GARCH = M_VAR_GARCH-True ) )
attach(new_data_clearRUS1)
new_data_clearRUS1 <- data.frame( ARMA_GARCH = ARMA_GARCH/M_VAR_GARCH , VAR_GARCH = VAR_GARCH/M_VAR_GARCH )
RUS <- new_data_clearRUS1[c(1:2) , ]     # M1-P1
#RUS <- new_data_clearRUS1[c(3:9) , ]    # M1-P2 

RUS1 <- data.frame( ID = as.character(new_data_clearRUS$Date)[ c(1:2) ] , RUS[ , 2 ] )
newdataRUS <- data.frame( Index = "RTSI" , melt( RUS1 , id = c("ID") ) )                               
#newdataRUS <- data.frame( Index = "RUS" , melt( data.frame( ID = "2020-06-11" , ARMA_GARCH = 0 , VAR_GARCH = 0 ) , id = c("ID") ) )   # M1-P3
                                
## FCHI
attach(new_data_clearFR)
new_data_clearFR1 <- abs( data.frame( ARMA_GARCH = ARMA_GARCH-True , VAR_GARCH = VAR_GARCH-True , M_VAR_GARCH = M_VAR_GARCH-True ) )
attach(new_data_clearFR1)
new_data_clearFR1 <- data.frame( ARMA_GARCH = ARMA_GARCH/M_VAR_GARCH , VAR_GARCH = VAR_GARCH/M_VAR_GARCH )
FR <- new_data_clearFR1[c(1:5) , ]        # M1-P1  
#FR <- new_data_clearFR1[c(6:9) , ]       # M1-P2                                
#FR <- new_data_clearFR1[c(10:11) , ]     # M1-P3

FR1<-data.frame( ID = as.character( new_data_clearFR$Date)[ c(1:5) ] , FR[ , 2 ] )
newdataFR<-data.frame( Index = "FCHI" , melt( FR1, id = c("ID") ) )
                                
## FTSE                               
attach(new_data_clearUK)
new_data_clearUK1 <- abs( data.frame( ARMA_GARCH = A RMA_GARCH-True , VAR_GARCH = VAR_GARCH-True , M_VAR_GARCH = M_VAR_GARCH-True ) )
attach(new_data_clearUK1)
new_data_clearUK1 <- data.frame( ARMA_GARCH = ARMA_GARCH/M_VAR_GARCH , VAR_GARCH = VAR_GARCH/M_VAR_GARCH )
UK <- new_data_clearUK1[c(1:2) , ]     # M1-P1  
#UK <- new_data_clearUK1[c(3:8) , ]    # M1-P2
#UK <- new_data_clearUK1[c(9:10) , ]   # M1-P3
                                
UK1 <- data.frame( ID = as.character( new_data_clearUK$Date)[ c(1:2) ] , UK[ , 2 ] )
newdataUK <- data.frame( Index = "FTSE" , melt( UK1 , id = c("ID") ) ) 
                                
## DAX
attach(new_data_clearBRD)
new_data_clearBRD1 <- abs( data.frame( ARMA_GARCH = ARMA_GARCH-True , VAR_GARCH = VAR_GARCH-True , M_VAR_GARCH = M_VAR_GARCH-True ) )
attach(new_data_clearBRD1)
new_data_clearBRD1 <- data.frame( ARMA_GARCH = ARMA_GARCH/M_VAR_GARCH , VAR_GARCH = VAR_GARCH/M_VAR_GARCH )
BRD <- new_data_clearBRD1[c(1:2) , ]     # M1-P1
#BRD <- new_data_clearBRD1[c(3:11) , ]   # M1-P2                               

BRD1 <- data.frame( ID = as.character( new_data_clearBRD$Date)[ c(1:2) ] , BRD[ , 2 ] )
newdataBRD <- data.frame( Index = "DAX" , melt( BRD1 , id = c("ID") ) )                                
#newdataBRD <- data.frame( Index = "DAX" , melt( data.frame( ID = "2020-06-11" , ARMA_GARCH = 0 , VAR_GARCH = 0 ) , id = c("ID") ) )   # M1-P3


new_data <- rbind( newdataS.Africa , newdataRUS , newdataFR , newdataUK , newdataBRD )
A        <- new_data[ order( as.Date( new_data$ID ) ) , ]
names(A) <- c("Index","ID","Method","R")

## P1
A$ID <- paste( A$ID , c( rep( "" , 26 ) , rep( "**" , 2 ) ) , sep="" )

## P2
#A$ID <- paste( A$ID , c( rep( "**" , 26 ) , rep( "*" , 38 ) ) , sep="" )

print(A %>% 
        ggplot( aes(x=ID, y=R) ) +
        geom_point(aes(shape = Index),alpha=0.6,size=10) +
        scale_shape_manual(values=c(20,8,10,6,7))+
        scale_y_continuous(limits = c(0.4, 1.8))+
        geom_abline(intercept = 1, slope = 0,colour="red")+
        theme_bw()+
        theme(axis.text.x = element_text(angle =-40, vjust = 1.5, hjust=-.2,size=20,face = "bold"),
              axis.text.y = element_text(size=12),
              legend.text = element_text(size=10),
              legend.title = element_text(size=15),
              axis.title.x=element_blank())
)                               
