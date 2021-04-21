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
source('R/PLOTS_Gibbs_Sampler.R')
Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")

# generate data
# data = generate_simulation_data(n_datasets = 50, n_time = 20, n_covariates = 10, seed = 71)

# data_gibbs = data
# B = 2000
# mu_initial = data$mu_true
# beta_initial = data$beta_true
# sigma_2_initial = data$sigma_2_true
# xi_initial = runif(length(data$xi_true), -1, 1)
# xi_initial = data$xi_true
# lK_initial = data$lK_true
# lB_initial = data$lB_true
# save(results, file = "results.RData")

data = generate_simulation_data(n_datasets = 300, n_time = 20, n_covariates = 10, seed = sample(1:10000, 1))
start = Sys.time()
results = gibbs_sampler_r(data_gibbs = data, 
                        B = 1000, 
                        # mu_initial = data$mu_true,
                        # beta_initial = data$beta_true,
                        # sigma_2_initial = data$sigma_2_true,
                        xi_initial = runif(length(data$xi_true), -1, 1),
                        # xi_initial = data$xi_true,
                        # lK_initial = data$lK_true,
                        lB_initial = data$lB_true)
saveRDS(results, file = "data")
plot_all(results)
betas = cbind(colMeans(results$beta), data$beta_true)
print(betas)
end = Sys.time()
print(end - start)


# mu_pm = colMeans(results$mu)
# beta_pm = colMeans(results$beta)
# sigma_2_pm = mean(results$sigma_2)
# sigma_2B_pm = mean(results$sigmaB_2)
# w_pm = colMeans(results$w)
# xi_pm = colMeans(results$xi)
# g_pm = colMeans(results$g)
# loglhood_pm = mean(results$loglhood)
# lk_pm = results$lK
# lb_pm = results$lB

