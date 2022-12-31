#!/usr/bin/env Rscript
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')
r <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) 

## sigma2 simulations
n_datasets_rep = seq(500, 2000, by = 500)
#n_datasets_rep = seq(50, 150, by = 50)

# generate data
idx = 1
data = generate_simulation_data(n_datasets = n_datasets_rep[idx], n_time = 10, n_covariates = 15, seed = 1, seed2 = r)
# # run MCMC
results = gibbs_sampler_r(data_gibbs = data, B = 30000, xi_initial = runif(length(data$xi_true), -1, 1), burn_in = 0.7, NNGP = FALSE, n_to_store = 10000)

# # # save to RDA file
filename = paste("RESULTS_NEW/results", n_datasets_rep[idx], "_", r, ".rda", sep = "")
save(results, file = filename)
print(r)
