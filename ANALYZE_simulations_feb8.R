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

T_reps = seq(20, 80, by = 20)
n_datasets_sim = 50

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

for (k in 1:length(T_reps))
{
  mm = matrix(0, ncol = 15, nrow = 50)
  vv = vector(length = 50)
  ll = matrix(0, nrow = 50, ncol = n_datasets_sim)
  gmu = matrix(0, nrow = 50, ncol = n_datasets_sim*T_reps[k])
  
  
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

for (idx in 1:length(T_reps))
{
  for (i in 1:50)
  {
    skip_to_next <- FALSE
    
    # getting the filename for the current .rda file
    temp_filename = paste("RESULTS/", T_reps[idx], "_", i, ".rda", sep = "")
    
    # reading in the data that the simulation would have read in - reading in csv
    data = read.csv(paste("t", T_idx[20]))
    
    # tryCatch(load(temp_filename), error = function(e) { skip_to_next <<- TRUE})
    
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
      mu = rep(colMeans(results$mu), each = T_reps[idx])
      
      g = colMeans(results$g)
      
      
      g_mu[[idx]][i,] = g + mu
      # CI = quantile(g + mu, probs = c(0.05, 0.95))
      # g_mu_low[[idx]][i,] = CI[1,]
      # g_mu_high[[idx]][i,] = CI[2,]
      # TODO fix this
      g_mu_truth[[idx]][i,] = unlist(data$g_true) + rep(data$mu_true, each = 10) 
      g_mu_bias[[idx]][i,] =  g_mu_truth[[idx]][i,] - g_mu[[idx]][i,]
      
    }
    print(paste(T_reps[idx], "/", i))
  }
}


# which ones to re-do
for (k in 1:length(T_reps))
{
  redo = which(rowMeans(beta[[k]]) == 0)
  if (length(redo) != 0)
  {
    print(paste("For t =", T_reps[k]))
    print(paste(redo, collapse = ","))
  }
}


write_filename = "output/"

for (i in 1:length(T_reps))
{
  write.csv(beta[[i]], row.names = FALSE, file = paste(write_filename, "beta_", T_reps[i], ".csv", sep = ""))
  write.csv(beta_truth[[i]], row.names = FALSE, file = paste(write_filename, "beta_truth_", T_reps[i], ".csv", sep = ""))
  write.csv(sigma[[i]], row.names = FALSE, file = paste(write_filename, "sigma_", T_reps[i], ".csv", sep = ""))
  write.csv(sigma_truth[[i]], row.names = FALSE, file = paste(write_filename, "sigma_truth_", T_reps[i], ".csv", sep = ""))
  write.csv(g_mu[[i]], row.names = FALSE, file = paste(write_filename, "gpmu_", T_reps[i], ".csv", sep = ""))
  write.csv(g_mu_truth[[i]], row.names = FALSE, file = paste(write_filename, "gmu_truth_", T_reps[i], ".csv", sep = ""))
}




