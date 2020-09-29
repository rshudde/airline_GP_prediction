rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')

# initial estimates for paramaters outside loop
knots = seq(0, 1, length.out = 500) # suggested to be 

# get data
n_covariates = 15
n_datasets = 500
data = generate_simulation_data(n_datasets = n_datasets, n_covariates = n_covariates, knots = knots)

results = gibbs_sampler(data_gibbs = data, knots_gibbs = knots, B = 5000, write = TRUE)