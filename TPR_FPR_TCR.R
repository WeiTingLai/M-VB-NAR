#============= TPR_FPR_CTR ==================================================
TPR_FPR_TCR<-function(object,coef_beta,beta,Iteration,current=F)
{
  m<-ncol(beta)
  Y<-object$Y
  p<-object$Lag
  if(current==T)
  {
    p<-p+1
  }else{
    if(nrow(beta)/m!=p)
    {
      beta =beta[-c(1:m),]
    }
  }
    coef_beta<-matrix(coef_beta,nrow=Iteration)
    beta<-as.matrix(beta)
    phi_true<-unlist(sapply(1:(2*p*m),function(j)ifelse(j%%2!=0,ifelse(beta[j%/%2+1,ifelse((j%/%2+1)%%ncol(Y)==0,ncol(Y),(j%/%2+1)%%ncol(Y))]==0,0,1),list(sapply(1:length(fit$segment[[j%/%2]]),function(s)ifelse(sum(abs(beta[j%/%2,ifelse(j%/%2%%ncol(Y)==0,-ncol(Y),-(j%/%2%%ncol(Y)))][c((sum(fit$segment[[j%/%2]][c(1:s)])-fit$segment[[j%/%2]][[s]]+1):sum(fit$segment[[j%/%2]][c(1:s)]))]))==0,0,1)))) ))
    TP<-sapply(1:Iteration,function(j)sum(ifelse(abs(coef_beta[j,which(c(abs(beta[c(1:(m*p)),]))>0)])>0,1,0)))
    FN<-sapply(1:Iteration,function(j)length(which(c(abs(beta[c(1:(m*p)),]))>0))-sum(ifelse(abs(coef_beta[j,which(c(abs(beta[c(1:(m*p)),]))>0)])>0,1,0)))
    FP<-sapply(1:Iteration,function(j)sum(ifelse(abs(coef_beta[j,which(c(abs(beta[c(1:(m*p)),]))==0)])>0,1,0)))
    TN<-sapply(1:Iteration,function(j)length(which(c(abs(beta[c(1:(m*p)),]))==0))-sum(ifelse(abs(coef_beta[j,which(c(abs(beta[c(1:(m*p)),]))==0)])>0,1,0)))

  
  TPRFPRTCR<-list()
  
  TPRFPRTCR$TPR<-mean(TP/(TP+FN))
  TPRFPRTCR$FPR<-mean(FP/(FP+TN))
  TPRFPRTCR$TNR<-mean(1-TPRFPRTCR$FPR)
  TPRFPRTCR$TCR<-mean((TP+TN)/length(c(beta[c(1:(m*p)),])))
  
  TPRFPRTCR$indicator_true<-sum(phi_true)

  TPRFPRTCR$coef_est<-mean(sapply(1:Iteration,function(j)(length(which(abs(coef_beta[j,])>0)))))
    
  }
  
  return(TPRFPRTCR)
}
