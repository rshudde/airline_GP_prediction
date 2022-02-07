#!/usr/bin/env Rscript
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')
r <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) 

read_data = function(i){
  X_initial = read.csv(file = paste0("t80/X_data_n100_t80_rep",i,".csv"), header = T)
  Y = as.matrix(read.csv(file = paste0("t80/Y_data_n100_t80_rep",i,".csv"), header = T))
  T_val = 80
  X = list()
  range = nrow(X_initial)/80
  k = 1
  for(j in 1:range){
    X[[j]] = as.matrix(X_initial[(k:(j*T_val)), ])
    k = (j*T_val)+1
    
  }
  
  return(list("X" = X, "y" = Y, "T_val" = T_val))
}

data = read_data(r)
print(paste("GOT DATA FOR CASE", r))

# run the simluations
results = gibbs_sampler_r(data_gibbs = data, 
                          B = 50000,
                          xi_initial = runif(data$T_val, -1, 1),
                          burn_in = 0.5,
                          NNGP = FALSE,
                          n_to_store = 15000)

filename = paste("RESULTS/results_n100_t80_rep", r, ".rda", sep = "")
save(results, file = filename)
