rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/PLOTS_Gibbs_sampler.R')

# initial estimates for paramaters outside loop
# get data
n_covariates = 6
n_datasets = 50
length_out = floor(n_covariates * n_datasets / 3)

knots = seq(0, 1, length.out = length_out) 
l_k_start = 4
sigma_2_start = 0.5
xi_start = rnorm(length(knots), 1, 4)

# generate data
data = generate_simulation_data(n_datasets = n_datasets, n_covariates = n_covariates, knots = knots, l_k = l_k_start, 
                                sigma_2 = sigma_2_start, xi_initial = xi_start)

# run gibbs sampler
xi_start = rnorm(length(knots), 0, 1)
results = gibbs_sampler(data_gibbs = data, knots_gibbs = knots, B = 2000, write = TRUE, lb_0 = 0.1, lk_0 = 0.1, burn_in = 0.3, 
                        xi_initial =  xi_start)

###################
plot_beta(data, results)
plot_sigma(data, results)
plot_lk(data, results)
plot_mu(data, results)
plot_xi(data, results)

# plots of true verses colmeans of the beta
a = data$beta
b = colMeans(results$beta)

par(mfrow = c(1,2))
plot(a, -b, main = "Actual verses -estimated B")
abline(a = 0, b = 1, col = "red")
plot(a, b)


