rm(list = ls())
source('FUNC_woodchan_samples.R')
source('FUNC_paramater_estimates.R')
source('DATA_generate_simulation.R')
source('FUNC_Gibbs_Sampler.R')
source('FUNC_Gibbs_Sampler_r.R')
source('PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

## sigma2 simulations
n_datasets_rep = seq(500, 1500, by = 500)

for (i in n_datasets_rep)
{
  for (r in 1:50)
  {
    # generate data
    data = generate_simulation_data(n_datasets = i, n_time = 10, n_covariates = 15, seed = r)
    # run MCMC
    results = gibbs_sampler_r(data_gibbs = data, B = 15000, xi_initial = runif(length(data$xi_true), -1, 1))
    
    # save to RDA file
    filename = paste("RESULTS/results_", i, "_", r, ".csv", sep = "")
    save(results, file = filename)
  }
}

# TO CHANGE
# set n_covariates = 15
