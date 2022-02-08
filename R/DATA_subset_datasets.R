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

t_vals = c(20,40,60,80)

for(i in 1:50){
  data = generate_simulation_data(n_datasets = 100, n_time = 80, n_covariates = 15, seed = i, seed2 = i, xi_true = 1)
  
  for(j in 1:length(t_vals)){
    set.seed(i) # seed needs to be set so it's 1-50, otherwise we'll have repeats
    idx = sample(1:80,t_vals[j], replace = FALSE)
    y_new = data$y[ ,idx]
    x_new = lapply(data$X, function(X) X[c(idx), ])
    data_new_x= do.call(rbind, x_new)
    write.csv(data_new_x, file = paste0("t", t_vals[j], "/X_data_n100_t",t_vals[j],"_rep",i,".csv"), row.names = FALSE)
    write.csv(y_new, file = paste0("t", t_vals[j], "/Y_data_n100_t",t_vals[j],"_rep",i,".csv"),row.names = FALSE)
    
    # now save all of the data
    data$idx = idx
    file_name = paste0("t", t_vals[j], "/ALL_data_n100_t",t_vals[j],"_rep",i,".rda")
    save(data, file = file_name)
    
    print(paste("JUST WROTE", i, "/", j))
  }
}
