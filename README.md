# Reproducibility for A Modified VAR-deGARCH Model for Asynchronous Multivariate Financial Time Series via Variational Bayesian Inference
This folder contains the R code for the Variational Bayesian Inference function for the Modified VAR-deGARCH (M-VAR-deGARCH) Model.

* This code can reproduce all the simulation results in Section 4.
  + The required R packages include `mvtnorm`, `gtools`, `boot`, `psych`, `Matrix`, `MCMCpack`, `MASS`, `RColorBrewer`, `lattice`, `gridExtra`, `ggplot2`, and `grid`.
  + The `Simulation.R` script reproduces the results from Section 4. The true coefficients are found in `BetaCoef_m20.csv`, and the data from 100 replications is available in `Datam20_rep100.csv`.
  + Since running 100 replications can be time-consuming, the result of the 100th replication is displayed in `Simulation.R`.

