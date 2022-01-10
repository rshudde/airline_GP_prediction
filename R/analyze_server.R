rm(list = ls())
library(matrixStats)
library(MLmetrics)
library(stats)
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

NNGP = FALSE
TEST = FALSE

if (TEST)
{
  n_datasets_rep = seq(50, 450, by = 50)
} else {
  n_datasets_rep = seq(500, 2000, by = 500)
}

beta = list()
sigma = list()
g_mu = list()

beta_truth = list()
sigma_truth = list()
g_mu_truth = list()

beta_high = list()
beta_low = list()

sigma_high = list()
sigma_low = list()

g_mu_high = list()
g_mu_low = list()


beta_MSE = list()
sigma_MSE = list()
g_mu_MSE = list()

beta_bias = list()
sigma_bias = list()
g_mu_bias = list()

for (k in 1:length(n_datasets_rep))
{
  mm = matrix(0, ncol = 15, nrow = 50)
  vv = vector(length = 50)
  ll = matrix(0, nrow = 50, ncol = n_datasets_rep[k])
  gmu = matrix(0, nrow = 50, ncol = n_datasets_rep[k]*10)
  
  
  beta[[k]] = mm
  sigma[[k]] = vv
  g_mu[[k]] = gmu
  
  beta_truth[[k]] = mm
  sigma_truth[[k]] = vv
  g_mu_truth[[k]] = gmu
  
  beta_high[[k]] = mm
  beta_low[[k]] = mm
  
  sigma_high[[k]] = vv
  sigma_low[[k]] = vv
  
  g_mu_high[[k]] = ll
  g_mu_low[[k]] = gmu
  
  
  beta_MSE[[k]] = mm
  sigma_MSE[[k]] = vv
  g_mu_MSE[[k]] = gmu
  
  beta_bias[[k]] = mm
  sigma_bias[[k]] = vv
  g_mu_bias[[k]] = gmu
}

for (idx in 1:length(n_datasets_rep))
{
  for (i in 1:50)
  {
    skip_to_next <- FALSE
    # 
    # if (NNGP) 
    # {
    #   temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS_NNGP/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    # } else {
    #   temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    # }
    
    if (NNGP & !TEST)
    {
      temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS_NNGP/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    } else if (!NNGP & !TEST) {
      temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    } else if (NNGP & TEST)
    {
      temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS_NNGP_TEST/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    } else {
      temp_filename = paste("~/Desktop/Airplane_Paper/RESULTS_TEST/results", n_datasets_rep[idx], "_", i, ".rda", sep = "")
    }
    
    
    data = generate_simulation_data(n_datasets = n_datasets_rep[idx], n_time = 10, n_covariates = 15, seed = 1, seed2 = i)
    
    tryCatch(load(temp_filename), error = function(e) { skip_to_next <<- TRUE})
    
    if (TEST & NNGP) results = results1
    
    if (!skip_to_next)
    {
      beta[[idx]][i,] = colMeans(results$beta)
      beta_truth[[idx]][i,] = data$beta
      CI = t(colQuantiles(results$beta, probs = c(0.05, 0.95)))
      beta_low[[idx]][i,] = CI[1,] # CI - lower
      beta_high[[idx]][i,] = CI[2,] # CI - higher
      beta_bias[[idx]][i, ] = data$beta - colMeans(results$beta)
      
      
      sigma[[idx]][i] = mean(results$sigma_2)
      sigma_truth[[idx]][i] = data$sigma_2_true
      CI = quantile(results$sigma_2, probs = c(0.05, 0.95))
      sigma_low[[idx]][i] = CI[1]
      sigma_high[[idx]][i] = CI[2]
      sigma_bias[[idx]][i] = data$sigma_2_true - mean(results$sigma_2)
      
      # now get g + mu
      mu = rep(colMeans(results$mu), each = 10)
      
      g = colMeans(results$g)
      
      
      g_mu[[idx]][i,] = g + mu
      # CI = quantile(g + mu, probs = c(0.05, 0.95))
      # g_mu_low[[idx]][i,] = CI[1,]
      # g_mu_high[[idx]][i,] = CI[2,]
      # TODO fix this
      g_mu_truth[[idx]][i,] = unlist(data$g_true) + rep(data$mu_true, each = 10) 
      g_mu_bias[[idx]][i,] =  g_mu_truth[[idx]][i,] - g_mu[[idx]][i,]
      
    }
    print(paste(n_datasets_rep[idx], "/", i))
  }
}

idx = 3

# round(colMeans(beta_bias[[idx]]),2)
# round(mean(sigma_bias[[idx]]),5)
# 
# high = colMeans(beta_high[[idx]])
# low = colMeans(beta_low[[idx]])
# total = round(cbind(low, high),2)
# 
# for (i in 1:nrow(total))
# {
#   print(paste("(", total[i,1], ",", total[i, 2], ")", sep = ""))
# }
# 
# print(paste("(", round(mean(sigma_low[[idx]]),2), ",", round(mean(sigma_high[[idx]]),2), ")", sep = ""))
# 
# # high = colMeans(g_mu_high[[idx]])
# # low = colMeans(g_mu_low[[idx]])
# # print(paste("(", round(mean(low),2), ",", round(mean(high),2), ")", sep = ""))
# 
# for (i in 1:15)
# {
#   temp = MSE(beta_truth[[idx]][,i], beta[[idx]][,i])
#   print(round(temp, 2))
# }
# 
# round(MSE(sigma[[idx]], sigma_truth[[idx]]), 7)


# which ones to re-do
for (k in 1:length(n_datasets_rep))
{
  redo = which(rowMeans(beta[[k]]) == 0)
  print(paste("For n =", n_datasets_rep[k]))
  print(paste(redo, collapse = ","))
}

# redo_500 = which(rowMeans(beta[[1]]) == 0)
# redo_1000 = which(rowMeans(beta[[2]]) == 0)
# redo_1500 = which(rowMeans(beta[[3]]) == 0)
# redo_2000 = which(rowMeans(beta[[4]]) == 0)
# 
# print(paste(redo_500, collapse = ","))
# print(paste(redo_1000, collapse = ","))
# print(paste(redo_1500, collapse = ","))
# print(paste(redo_2000, collapse = ","))




if (NNGP & !TEST)
{
  write_filename = "summary_stats/NNGP/"
} else if (!NNGP & !TEST) {
  write_filename = "summary_stats/NORMAL/"
} else if (NNGP & TEST)
{
  write_filename = "summary_stats/NNGP_TEST/"
} else {
  write_filename = "summary_stats/NORMAL_TEST/"
}


for (i in 1:length(n_datasets_rep))
{
  write.csv(beta[[i]], row.names = FALSE, file = paste(write_filename, "beta_", n_datasets_rep[i], ".csv", sep = ""))
  write.csv(beta_truth[[i]], row.names = FALSE, file = paste(write_filename, "beta_truth_500", n_datasets_rep[i], ".csv", sep = ""))
  write.csv(sigma[[i]], row.names = FALSE, file = paste(write_filename, "sigma_", n_datasets_rep[i], ".csv", sep = ""))
  write.csv(sigma_truth[[i]], row.names = FALSE, file = paste(write_filename, "sigma_truth_", n_datasets_rep[i], ".csv", sep = ""))
  write.csv(g_mu[[i]], row.names = FALSE, file = paste(write_filename, "mu_", n_datasets_rep[i], ".csv", sep = ""))
  write.csv(g_mu_truth[[i]], row.names = FALSE, file = paste(write_filename, "gmu_truth_", n_datasets_rep[i], ".csv", sep = ""))
}




