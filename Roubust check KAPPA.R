### Kappa lag-0  
## Since each method is computationally expensive, we preserve the results for different numbers of time series.
## coefficient matrix

## Case 1 in Table 6
coef_1_12 <- read.csv("D:/coef unfixed 0.001 (2015-3_2018-3)_new (12).csv")
coef_2_12 <- read.csv("D:/coef unfixed 0.001 (2016-3_2019-3)_new (12).csv")
coef_3_12 <- read.csv("D:/coef unfixed 0.001 (2017-3_2020-3)_new (12).csv")

## Case 2 in Table 6
coef_1_12_2 <- read.csv("D:/coef unfixed 0.001 (2015-3_2018-3) (12-N255+STI).csv")
coef_2_12_2 <- read.csv("D:/coef unfixed 0.001 (2016-3_2019-3) (12-N255+STI).csv")
coef_3_12_2 <- read.csv("D:/coef unfixed 0.001 (2017-3_2020-3) (12-N255+STI).csv")

## Case 3 in Table 6
coef_1_12_3 <- read.csv("D:/coef unfixed 0.001 (2015-3_2018-3) (12-HIS+STI).csv")
coef_2_12_3 <- read.csv("D:/coef unfixed 0.001 (2016-3_2019-3) (12-HIS+STI).csv")
coef_3_12_3 <- read.csv("D:/coef unfixed 0.001 (2017-3_2020-3) (12-HIS+STI).csv")

## 20 time series
coef_1_20 <- read.csv("D:/coef unfixed 0.001 (2015-3_2018-3)_new (20).csv")
coef_2_20 <- read.csv("D:/coef unfixed 0.001 (2016-3_2019-3)_new (20).csv")
coef_3_20 <- read.csv("D:/coef unfixed 0.001 (2017-3_2020-3)_new (20).csv")

##  I(A0) in Table 2
adjma <- read.csv("D:/opentime.csv")


#######################################  Contingency Table  ##############################################

## Case1: N225 HSI JTOPI RSI FCHI FSTSE DAX BVSP S&P500 MXX GSPTSE MERV 

adjma_1 <- adjma[c(2,5,11:20),c(3,6,12:21)]

coef_1_20_lag0 <- coef_1_20[c(2,5,11:20),c(2,5,11:20)]
coef_2_20_lag0 <- coef_2_20[c(2,5,11:20),c(2,5,11:20)]
coef_3_20_lag0 <- coef_3_20[c(2,5,11:20),c(2,5,11:20)]

coef_1_12_lag0 <- coef_1_12[c(1:12),c(1:12)]
coef_2_12_lag0 <- coef_2_12[c(1:12),c(1:12)]
coef_3_12_lag0 <- coef_3_12[c(1:12),c(1:12)]

## Active Coefficient 

activecoef_1_12_lag0 <- (abs(coef_1_12_lag0)>0)*1
activecoef_2_12_lag0 <- (abs(coef_2_12_lag0)>0)*1
activecoef_3_12_lag0 <- (abs(coef_3_12_lag0)>0)*1

activecoef_1_20_lag0 <- (abs(coef_1_20_lag0)>0)*1
activecoef_2_20_lag0 <- (abs(coef_2_20_lag0)>0)*1
activecoef_3_20_lag0 <- (abs(coef_3_20_lag0)>0)*1

#####################################  M1  #######################################

activecoef_1_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_1_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_1_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_1_20_lag0[j,which(adjma_1[j,]>0)]))

Table_1  <- table(activecoef_1_12_lag0_new,activecoef_1_20_lag0new)
sum_1    <- sum(Table_1)
coef12_1 <- apply(Table_1,1,sum)
coef20_1 <- apply(Table_1,2,sum)

P0_1     <- sum(diag(Table_1))/sum_1
Pe_1     <- sum(sapply(1:length(coef20_1),function(k)coef20_1[k]*coef12_1[k]))/(sum_1)^2

Kappa_1  <- (P0_1-Pe_1)/(1-Pe_1);round(Kappa_1,3)

#####################################  M2  #######################################

activecoef_2_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_2_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_2_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_2_20_lag0[j,which(adjma_1[j,]>0)]))


Table_2  <- table(activecoef_2_12_lag0_new,activecoef_2_20_lag0new)
sum_2    <- sum(Table_2)
coef12_2 <- apply(Table_2,1,sum)
coef20_2 <- apply(Table_2,2,sum)

P0_2     <- sum(diag(Table_2))/sum_2
Pe_2     <- sum(sapply(1:length(coef20_2),function(k)coef20_2[k]*coef12_2[k]))/(sum_2)^2

Kappa_2  <- (P0_2-Pe_2)/(1-Pe_2);round(Kappa_2,3)

#####################################  M3  #######################################   

activecoef_3_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_3_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_3_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_3_20_lag0[j,which(adjma_1[j,]>0)]))

Table_3  <- table(activecoef_3_12_lag0_new,activecoef_3_20_lag0new)
sum_3    <- sum(Table_3)
coef12_3 <- apply(Table_3,1,sum)
coef20_3 <- apply(Table_3,2,sum)

P0_3     <- sum(diag(Table_3))/sum_3
Pe_3     <- sum(sapply(1:length(coef20_3),function(k)coef20_3[k]*coef12_3[k]))/(sum_3)^2

Kappa_3  <- (P0_3-Pe_3)/(1-Pe_3);round(Kappa_3,3)

#######################################  Contingency Table  ##############################################

## Case2: HSI STI JTOPI RSI FCHI FSTSE DAX BVSP S&P500 MXX GSPTSE MERV 

adjma_2 <- adjma[c(5,7,11:20),c(6,8,12:21)]

coef_1_12_lag0_1 <- coef_1_12_2[c(1:12),c(1:12)]
coef_2_12_lag0_1 <- coef_2_12_2[c(1:12),c(1:12)]
coef_3_12_lag0_1 <- coef_3_12_2[c(1:12),c(1:12)]

coef_1_20_lag0   <- coef_1_20[c(5,7,11:20),c(5,7,11:20)]
coef_2_20_lag0   <- coef_2_20[c(5,7,11:20),c(5,7,11:20)]
coef_3_20_lag0   <- coef_3_20[c(5,7,11:20),c(5,7,11:20)]

## Active Coefficient 

activecoef_1_12_lag0_1 <- (abs(coef_1_12_lag0_1)>0)*1
activecoef_2_12_lag0_1 <- (abs(coef_2_12_lag0_1)>0)*1
activecoef_3_12_lag0_1 <- (abs(coef_3_12_lag0_1)>0)*1

activecoef_1_20_lag0   <- (abs(coef_1_20_lag0)>0)*1
activecoef_2_20_lag0   <- (abs(coef_2_20_lag0)>0)*1
activecoef_3_20_lag0   <- (abs(coef_3_20_lag0)>0)*1

#####################################  M1  #######################################

activecoef_1_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_1_12_lag0_1[j,which(adjma_2[j,]>0)]))
activecoef_1_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_1_20_lag0[j,which(adjma_2[j,]>0)]))

Table_11  <- table(activecoef_1_12_lag0_1new,activecoef_1_20_lag0new)
sum_11    <- sum(Table_11)
coef10_11 <- apply(Table_11,1,sum)
coef20_11 <- apply(Table_11,2,sum)

P0_11     <- sum(diag(Table_11))/sum_11
Pe_11     <- sum(sapply(1:length(coef20_11),function(k)coef20_11[k]*coef10_11[k]))/(sum_11)^2

Kappa_11  <- (P0_11-Pe_11)/(1-Pe_11);round(Kappa_11,3)

#####################################  M2  #######################################

activecoef_2_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_2_12_lag0_1[j,which(adjma_2[j,]>0)]))
activecoef_2_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_2_20_lag0[j,which(adjma_2[j,]>0)]))

Table_21  <- table(activecoef_2_12_lag0_1new,activecoef_2_20_lag0new)
sum_21    <- sum(Table_21)
coef10_21 <- apply(Table_21,1,sum)
coef20_21 <- apply(Table_21,2,sum)

P0_21     <- sum(diag(Table_21))/sum_21
Pe_21     <- sum(sapply(1:length(coef20_21),function(k)coef20_21[k]*coef10_21[k]))/(sum_21)^2

Kappa_21  <- (P0_21-Pe_21)/(1-Pe_21);round(Kappa_21,3)

#####################################  M3  #######################################   

activecoef_3_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_3_12_lag0_1[j,which(adjma_2[j,]>0)]))
activecoef_3_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_3_20_lag0[j,which(adjma_2[j,]>0)]))

Table_31  <- table(activecoef_3_12_lag0_1new,activecoef_3_20_lag0new)
sum_31    <- sum(Table_31)
coef10_31 <- apply(Table_31,1,sum)
coef20_31 <- apply(Table_31,2,sum)

P0_31     <- sum(diag(Table_31))/sum_31
Pe_31     <- sum(sapply(1:length(coef20_31),function(k)coef20_31[k]*coef10_31[k]))/(sum_31)^2

Kappa_31  <- (P0_31-Pe_31)/(1-Pe_31);round(Kappa_31,3)


#######################################  Contingency Table  ##############################################

## Case1: N225 HSI JTOPI RSI FCHI FSTSE DAX BVSP S&P500 MXX GSPTSE MERV 

adjma_1 <- adjma[c(2,5,11:20),c(3,6,12:21)]

coef_1_20_lag0 <- coef_1_20[c(2,5,11:20),c(2,5,11:20)]
coef_2_20_lag0 <- coef_2_20[c(2,5,11:20),c(2,5,11:20)]
coef_3_20_lag0 <- coef_3_20[c(2,5,11:20),c(2,5,11:20)]

coef_1_12_lag0 <- coef_1_12[c(1:12),c(1:12)]
coef_2_12_lag0 <- coef_2_12[c(1:12),c(1:12)]
coef_3_12_lag0 <- coef_3_12[c(1:12),c(1:12)]

## Active Coefficient 

activecoef_1_12_lag0 <- (abs(coef_1_12_lag0)>0)*1
activecoef_2_12_lag0 <- (abs(coef_2_12_lag0)>0)*1
activecoef_3_12_lag0 <- (abs(coef_3_12_lag0)>0)*1

activecoef_1_20_lag0 <- (abs(coef_1_20_lag0)>0)*1
activecoef_2_20_lag0 <- (abs(coef_2_20_lag0)>0)*1
activecoef_3_20_lag0 <- (abs(coef_3_20_lag0)>0)*1

#####################################  M1  #######################################

activecoef_1_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_1_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_1_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_1_20_lag0[j,which(adjma_1[j,]>0)]))

Table_1  <- table(activecoef_1_12_lag0_new,activecoef_1_20_lag0new)
sum_1    <- sum(Table_1)
coef12_1 <- apply(Table_1,1,sum)
coef20_1 <- apply(Table_1,2,sum)

P0_1     <- sum(diag(Table_1))/sum_1
Pe_1     <- sum(sapply(1:length(coef20_1),function(k)coef20_1[k]*coef12_1[k]))/(sum_1)^2

Kappa_1  <- (P0_1-Pe_1)/(1-Pe_1);round(Kappa_1,3)

#####################################  M2  #######################################

activecoef_2_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_2_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_2_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_2_20_lag0[j,which(adjma_1[j,]>0)]))


Table_2  <- table(activecoef_2_12_lag0_new,activecoef_2_20_lag0new)
sum_2    <- sum(Table_2)
coef12_2 <- apply(Table_2,1,sum)
coef20_2 <- apply(Table_2,2,sum)

P0_2     <- sum(diag(Table_2))/sum_2
Pe_2     <- sum(sapply(1:length(coef20_2),function(k)coef20_2[k]*coef12_2[k]))/(sum_2)^2

Kappa_2  <- (P0_2-Pe_2)/(1-Pe_2);round(Kappa_2,3)

#####################################  M3  #######################################   

activecoef_3_12_lag0_new <- unlist(sapply(1:12,function(j)activecoef_3_12_lag0[j,which(adjma_1[j,]>0)]))
activecoef_3_20_lag0new  <- unlist(sapply(1:12,function(j)activecoef_3_20_lag0[j,which(adjma_1[j,]>0)]))

Table_3  <- table(activecoef_3_12_lag0_new,activecoef_3_20_lag0new)
sum_3    <- sum(Table_3)
coef12_3 <- apply(Table_3,1,sum)
coef20_3 <- apply(Table_3,2,sum)

P0_3     <- sum(diag(Table_3))/sum_3
Pe_3     <- sum(sapply(1:length(coef20_3),function(k)coef20_3[k]*coef12_3[k]))/(sum_3)^2

Kappa_3  <- (P0_3-Pe_3)/(1-Pe_3);round(Kappa_3,3)

#######################################  Contingency Table  ##############################################

## Case3: N255 STI JTOPI RSI FCHI FSTSE DAX BVSP S&P500 MXX GSPTSE MERV 

adjma_3 <- adjma[c(2,7,11:20),c(3,8,12:21)]

coef_1_12_lag0_1 <- coef_1_12_3[c(1:12),c(1:12)]
coef_2_12_lag0_1 <- coef_2_12_3[c(1:12),c(1:12)]
coef_3_12_lag0_1 <- coef_3_12_3[c(1:12),c(1:12)]

coef_1_20_lag0  <- coef_1_20[c(2,7,11:20),c(2,7,11:20)]
coef_2_20_lag0  <- coef_2_20[c(2,7,11:20),c(2,7,11:20)]
coef_3_20_lag0  <- coef_3_20[c(2,7,11:20),c(2,7,11:20)]

## Active Coefficient 

activecoef_1_12_lag0_1 <- (abs(coef_1_12_lag0_1)>0)*1
activecoef_2_12_lag0_1 <- (abs(coef_2_12_lag0_1)>0)*1
activecoef_3_12_lag0_1 <- (abs(coef_3_12_lag0_1)>0)*1

activecoef_1_20_lag0   <- (abs(coef_1_20_lag0)>0)*1
activecoef_2_20_lag0   <- (abs(coef_2_20_lag0)>0)*1
activecoef_3_20_lag0   <- (abs(coef_3_20_lag0)>0)*1

#####################################  M1  #######################################

activecoef_1_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_1_12_lag0_1[j,which(adjma_3[j,]>0)]))
activecoef_1_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_1_20_lag0[j,which(adjma_3[j,]>0)]))

Table_12  <- table(activecoef_1_12_lag0_1new,activecoef_1_20_lag0new)
sum_12    <- sum(Table_12)
coef10_12 <- apply(Table_12,1,sum)
coef12_12 <- apply(Table_12,2,sum)

P0_12     <- sum(diag(Table_12))/sum_12
Pe_12     <- sum(sapply(1:length(coef12_12),function(k)coef12_12[k]*coef10_12[k]))/(sum_12)^2

Kappa_12  <- (P0_12-Pe_12)/(1-Pe_12);round(Kappa_12,3)

#####################################  M2  #######################################

activecoef_2_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_2_12_lag0_1[j,which(adjma_3[j,]>0)]))
activecoef_2_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_2_20_lag0[j,which(adjma_3[j,]>0)]))

Table_22  <- table(activecoef_2_12_lag0_1new,activecoef_2_20_lag0new)
sum_22    <- sum(Table_22)
coef10_22 <- apply(Table_22,1,sum)
coef12_22 <- apply(Table_22,2,sum)

P0_22     <- sum(diag(Table_22))/sum_22
Pe_22     <- sum(sapply(1:length(coef12_22),function(k)coef12_22[k]*coef10_22[k]))/(sum_22)^2

Kappa_22  <- (P0_22-Pe_22)/(1-Pe_22);round(Kappa_22,3)

#####################################  M3  #######################################   

activecoef_3_12_lag0_1new <- unlist(sapply(1:12,function(j)activecoef_3_12_lag0_1[j,which(adjma_3[j,]>0)]))
activecoef_3_20_lag0new   <- unlist(sapply(1:12,function(j)activecoef_3_20_lag0[j,which(adjma_3[j,]>0)]))

Table_32  <- table(activecoef_3_12_lag0_1new,activecoef_3_20_lag0new)
sum_32    <- sum(Table_32)
coef10_32 <- apply(Table_32,1,sum)
coef12_32 <- apply(Table_32,2,sum)

P0_32     <- sum(diag(Table_32))/sum_32
Pe_32     <- sum(sapply(1:length(coef12_32),function(k)coef12_32[k]*coef10_32[k]))/(sum_32)^2

Kappa_32  <- (P0_32-Pe_32)/(1-Pe_32);round(Kappa_32,3)
