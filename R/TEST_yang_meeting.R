rm(list = ls())
library(invgamma)
library(MASS)
library(nlme)
library(invgamma)
library(MASS)
library(FastGP)
library(emulator)

source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")


data = generate_simulation_data(n_datasets = 100, n_time = 80, n_covariates = 10, seed = 1, seed2 = 1)
print("GOT DATA")

results = gibbs_sampler_r(data_gibbs = data, B = 100, 
                          xi_initial = runif(length(data$w_true[[1]]), -1, 1), 
                          burn_in = 0.7, 
                          NNGP = FALSE, 
                          n_to_store = 100)
save(results, file = "TESTRESULTS.rda")