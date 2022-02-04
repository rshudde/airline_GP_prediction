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


read_data = function(i){
  X_initial = read.csv(file = paste0("t20/X_data_n100_t20_rep",i,".csv"), header = T)
  Y = read.csv(file = paste0("t20/Y_data_n100_t20_rep",i,".csv"), header = T)
  T = 20
  X = list()
  range = nrow(X_initial)/20
  k = 1
  for(j in 1:range){
    X[[j]] = X_initial[(k:(j*T)), ]
    k = (j*T)+1
    
  }
  
  return(list("X" = X, "Y" = Y, "T" = T))
}

data = read_data(i)
print("GOT DATA")
results = gibbs_sampler_r(data_gibbs = data, B = 100000,
                          xi_initial = runif(data$T, -1, 1),
                          burn_in = 0.5,
                          NNGP = FALSE,
                          n_to_store = 20000)
save(results, file = "TESTRESULTS100.rda")


