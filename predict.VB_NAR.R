#============= Predict y_hat given X and object model====================
predict.VB_NAR<-function(object,Y,X=NULL, step_ahead=1,current=F,rolling=T,rolling_num=N){
 p<-object$Lag
  if(current==T)
  {
    bbhat<-matrix(0,nrow=nrow(Y)-1,ncol=m*m*(p+1))
  }else{
    bbhat<-matrix(0,nrow=nrow(Y)-1,ncol=m*m*(p))
    
  }
  Y1<-as.matrix(Y)
  m<-ncol(Y)
  if(current==T)
  {
    if(is.null(X))
    {
      Y<-rbind(as.matrix(object$Y),Y1)
      Y<-Y[c((nrow(Y)-nrow(Y1)-p+1):nrow(Y)),]
      X1<-t(sapply((p+1):nrow(Y),function(j)c(t(Y[c((j-1):(j-p)),]))))
      X <- cbind(Y1,X1)
    }
    X<-rbind(object$X,X)
    if(step_ahead==1)
    {
      y_hat<- X[nrow(X),]%*%matrix(object$coef_beta,ncol=m)
    }else{
      y_hat<- X[(nrow(X)-step_ahead+1),]%*%matrix(object$coef_beta,ncol=m)
      for(i in 2:step_ahead){
        bb<-matrix(0,nrow=m*(p+1),ncol=m)
        phi_idx<-object$phi
        for(j in 1:m)
        {
          if(length(which(phi_idx[,j]>0))>0)
          {
            if(rolling==T)
            {
              bb[c(which(phi_idx[,j]>0)),j]<-ginv(crossprod(X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))]))%*%crossprod(X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],rbind(matrix(object$Y[,j],ncol=1),matrix(Y1[c(1:(i-1)),j],ncol=1))[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),])
            }else{
              bb[c(which(phi_idx[,j]>0)),j]<-solve(crossprod(X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))]))%*%crossprod(X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],rbind(matrix(object$Y[,j],ncol=1),matrix(Y1[c(1:(i-1)),j],ncol=1)))
              
            }
          }
        }
        bbhat[i-1,]<-bb
        y_hat<-rbind(y_hat, X[(nrow(object$X)+i),]%*%bb)
      }
    }
  }else{
    if(is.null(X))
    {
      Y<-rbind(as.matrix(object$Y),Y1)
      Y<-Y[c((nrow(Y)-nrow(Y1)-p+1):nrow(Y)),]
      X<-t(sapply((p+1):nrow(Y),function(j)c(t(Y[c((j-1):(j-p)),]))))
    }
    X<-rbind(object$X,X)
    if(step_ahead==1)
    {
      y_hat<- X[nrow(X),]%*%matrix(object$coef_beta,ncol=m)
    }else{
      y_hat<- X[(nrow(X)-step_ahead+1),]%*%matrix(object$coef_beta,ncol=m)
      for(i in 2:step_ahead){
        bb<-matrix(0,nrow=m*p,ncol=m)
        phi_idx<-object$phi
        for(j in 1:m)
        {
          if(length(which(phi_idx[,j]>0))>0)
          {
            if(rolling==T)
            {
              bb[c(which(phi_idx[,j]>0)),j]<-ginv(crossprod(X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))]))%*%crossprod(X[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],rbind(matrix(object$Y[,j],ncol=1),matrix(Y1[c(1:(i-1)),j],ncol=1))[c((nrow(object$X)+i-rolling_num):(nrow(object$X)+i-1)),])
            }else{
              bb[c(which(phi_idx[,j]>0)),j]<-solve(crossprod(X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))]))%*%crossprod(X[c(1:(nrow(object$X)+i-1)),c(which(phi_idx[,j]>0))],rbind(matrix(object$Y[,j],ncol=1),matrix(Y1[c(1:(i-1)),j],ncol=1)))
              
            }
          }
        }
        bbhat<-rbind(bbhat,c(bb))
        y_hat<-rbind(y_hat, X[(nrow(object$X)+i),]%*%bb)
      }
    }
    
  }
  
  
  fit<-list()
  fit$y_hat<-y_hat
  fit$bb<-bbhat
  
  return(fit)
  
  
  
}

