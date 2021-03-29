rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_Sampler.R')


Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")

# generate data
data = generate_simulation_data(n_datasets = 30, n_time = 10, 
                                n_covariates = 10, seed = 2)

# data_gibbs = data 
# B = 2000 
# mu_initial = data$mu_true
# beta_initial = data$beta_true
# sigma_2_initial = data$sigma_2_true
# xi_initial = runif(length(data$xi_true), -1, 1)
# xi_initial = data$xi_true
# lK_initial = data$lK_true
# lB_initial = data$lB_true

# run gibbs sampler
microbenchmark::microbenchmark(
results = gibbs_sampler(data_gibbs = data, 
                        B = 1000, 
                        mu_initial = data$mu_true,
                        # beta_initial = data$beta_true,
                        sigma_2_initial = data$sigma_2_true,
                        xi_initial = runif(length(data$xi_true), -1, 1),
                        # xi_initial = data$xi_true,
                        lK_initial = data$lK_true,
                        lB_initial = data$lB_true)
, times = 1)

microbenchmark::microbenchmark(
  results = gibbs_sampler_r(data_gibbs = data, 
                          B = 1000, 
                          mu_initial = data$mu_true,
                          # beta_initial = data$beta_true,
                          sigma_2_initial = data$sigma_2_true,
                          xi_initial = runif(length(data$xi_true), -1, 1),
                          # xi_initial = data$xi_true,
                          lK_initial = data$lK_true,
                          lB_initial = data$lB_true)
, times = 1)

# save(results, file = "results.RData")

results = gibbs_sampler_r(data_gibbs = data, 
                        B = 1000, 
                        mu_initial = data$mu_true,
                        # beta_initial = data$beta_true,
                        sigma_2_initial = data$sigma_2_true,
                        xi_initial = runif(length(data$xi_true), -1, 1),
                        # xi_initial = data$xi_true,
                        lK_initial = data$lK_true,
                        lB_initial = data$lB_true)
plot_all(results)

mu_pm = colMeans(results$mu)
beta_pm = colMeans(results$beta)
sigma_2_pm = mean(results$sigma_2)
sigma_2B_pm = mean(results$sigmaB_2)
w_pm = colMeans(results$w)
xi_pm = colMeans(results$xi)
g_pm = colMeans(results$g)
loglhood_pm = mean(results$loglhood)
lk_pm = results$lK
lb_pm = results$lB

