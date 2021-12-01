rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')
r <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) 

## sigma2 simulations
n_datasets_rep = seq(500, 1500, by = 500)
# n_datasets_rep = seq(20, 40, by = 10)

idx = 1
r = 1
# generate data
data = generate_simulation_data(n_datasets = n_datasets_rep[idx], n_time = 10, n_covariates = 15, seed = 1, seed2 = r)
# # run MCMC
results = gibbs_sampler_r(data_gibbs = data, B = 30000, xi_initial = runif(length(data$xi_true), -1, 1), burn_in = 0.7, NNGP = FALSE, n_to_store = 10000)
# # # save to RDA file
filename = paste("RESULTS_TEST/TESTresults_", n_datasets_rep[1], "_", r, "rda", sep = "")
save(results, file = filename)
print(r)
  