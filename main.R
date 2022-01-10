rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')




n_datasets_sim = 100
data = generate_simulation_data(n_datasets = n_datasets_sim, n_time = 20, n_covariates = 15, seed = 1, seed2 = 1)

results = gibbs_sampler_r(data_gibbs = data, B = 30000, xi_initial = runif(length(data$w_true[[1]]), -1, 1), burn_in = 0.7, NNGP = FALSE, n_to_store = 10000)

# # generate data
# i = 1
# for (i in 1:25)
# {
#   
#   data = generate_simulation_data(n_datasets = 200, n_time = 10,
#                                   n_covariates = 10, seed = i)
#   
#   # run gibbs sampler
#   start = Sys.time()
#   results = gibbs_sampler_r(data_gibbs = data, B = 5000,
#                           # mu_initial = data$mu_true,
#                           # beta_initial = data$beta_true,
#                           # sigma_2_initial = data$sigma_2_true,
#                           xi_initial = runif(length(data$xi_true), -1, 1),
#                           # xi_initial = data$xi_true,
#                           # lK_initial = data$lK_true,
#                           # lB_initial = data$lB_true,
#                           burn_in = 0.7,
#                           NNGP = FALSE
#                           )
#   end = Sys.time()
#   print(end - start)
#   plot_all(results)
#   # par(mfrow = c(1,1));plot(results$lB, type = "l")
# }
# 
# g1 = colMeans(results$g)
# mu = colMeans(results$mu)
# 
# 
# g = vector()
# count = 1
# for (i in seq(1, length(g1), by = 10))
# {
#   g[count] = mean(g1[i:(i+9)])
#   count = count + 1
# }
# g
# 
# data = cbind(g, mu)
