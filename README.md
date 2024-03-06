# A Modified VAR-deGARCH Model for Asynchronous Multivariate Financial Time Series via Variational Bayesian Inference
This code is the function of the Variational Bayesian Inference for the Modified VAR-deGARCH (M-VAR-deGARCH) Model.
All the simulation results can be reproduced by using this code.

# Simulation:
```
## Read Data

>Beta<-as.matrix(read.table("D:/BetaCoef_m20.csv",sep=",",head=F))
>Data<-as.matrix(read.table("D:/Datam20_rep100.csv",sep=",",head=F))

## Set all the parameters:
```
# References
Lai, W.T., Chen, R.B., Chen, Y., Koch, T. (2022). Variational Bayesian inference for network autoregression models. Computational Statistics & Data Analysis 169, 107406.
