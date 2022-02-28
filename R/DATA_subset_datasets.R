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

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  stop("NO COMMAND LINE ARGUMENTS PASSED FOR R AND T")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

for (i in 1:n_replicates)
  {
  data = generate_simulation_data(n_datasets = 100, n_time = num_flights, n_covariates = 15, seed = i, seed2 = i, xi_true = 1)
  
  set.seed(i) # seed needs to be set so it's 1-50, otherwise we'll have repeats
  idx = sample(1:num_flights, t_vals, replace = FALSE)
  y_new = data$y[ ,idx]
  x_new = lapply(data$X, function(X) X[c(idx), ])
  data_new_x = do.call(rbind, x_new)
  write.csv(data_new_x, file = paste0("t", t_vals, "/X_data_n100_t",t_vals,"_rep",i,".csv"), row.names = FALSE)
  write.csv(y_new, file = paste0("t", t_vals, "/Y_data_n100_t",t_vals,"_rep",i,".csv"),row.names = FALSE)
  
  print(paste("JUST WROTE", t_vals, "/", i))

}