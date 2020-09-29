rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')

# initial estimates for paramaters outside loop
# get data
n_covariates = 10
n_datasets = 500
length_out = floor(n_covariates*n_datasets / 3)

knots = seq(0, 1, length.out = length_out) # suggested to be n_covariates*n_datasets / 2

# generate data
data = generate_simulation_data(n_datasets = n_datasets, n_covariates = n_covariates, knots = knots)

# run gibbs sampler
results = gibbs_sampler(data_gibbs = data, knots_gibbs = knots, B = 1000, write = TRUE, lb_0 = 0.1, lk_0 = 0.1)



