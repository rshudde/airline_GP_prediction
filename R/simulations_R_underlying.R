#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

# things passed in 
# r <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) 
# current_t = ""

r = as.numeric(ags[1])
current_t = as.numeric(args[2])

# function to read data
read_data = function(i, t_val){
  X_initial = read.csv(file = paste0("t", t_val, "/X_data_n100_t", t_val, "_rep",i,".csv"), header = T)
  Y = as.matrix(read.csv(file = paste0("t", t_val, "/Y_data_n100_t", t_val, "_rep",i,".csv"), header = T))
  X = list()
  range = nrow(X_initial)/t_val
  k = 1
  for(j in 1:range){
    X[[j]] = as.matrix(X_initial[(k:(j*t_val)), ])
    k = (j*t_val)+1
    
  }
  
  return(list("X" = X, "y" = Y, "t_val" = t_val))
}

data = read_data(r, current_t)
print(paste("GOT DATA FOR CASE", current_t, "/", r))

# run the simluations
results = gibbs_sampler_r(data_gibbs = data, 
                          B = 500,
                          xi_initial = runif(data$t_val, -1, 1),
                          burn_in = 0.5,
                          NNGP = FALSE,
                          n_to_store = 200)

filename = paste("RESULTS/results_n100_t", current_t, "_rep", r, ".rda", sep = "")
save(results, file = filename)
