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
# Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  stop("NO COMMAND LINE ARGUMENTS PASSED FOR R AND T")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
# t_vals = 20

n_datasets_sim = num_flights # how many flights are being observed

# set up size of all initial datasets
mm = matrix(0, ncol = 15, nrow = n_replicates)
vv = vector(length = n_replicates)
ll = matrix(0, nrow = n_replicates, ncol = n_datasets_sim)
gmu = matrix(0, nrow = n_replicates, ncol = n_datasets_sim*t_vals)

beta = mm
sigma = vv
g_mu = gmu

beta_truth = mm
sigma_truth = vv
g_mu_truth = gmu

beta_high = mm
beta_low = mm

sigma_high = vv
sigma_low = vv

g_mu_high = ll
g_mu_low = gmu


beta_MSE = mm
sigma_MSE = vv
g_mu_MSE = gmu

beta_bias = mm
sigma_bias = vv
g_mu_bias = gmu

timing = vector(length = n_replicates)


for (i in 1:n_replicates)
{
  skip_to_next <- FALSE
  
  # getting the filename for the current .rda file
  # temp_filename = paste("RESULTS/", T_reps[idx], "_", i, ".rda", sep = "")
  if (USE_NNGP)
  {
    temp_filename = paste("RESULTSNNGP/results_n100_t", t_vals, "_rep", i, "_NNGP.rda", sep = "")
  } else {
    temp_filename = paste("RESULTS/results_n100_t", t_vals, "_rep", i, ".rda", sep = "")
  }
  
  # reading in the data that the simulation would have read in - reading in csv
  data = generate_simulation_data(n_datasets = 100, n_time = t_vals, n_covariates = 15, seed = i, seed2 = i, xi_true = 1)
  
  tryCatch(load(temp_filename), error = function(e) { skip_to_next <<- TRUE})
  
  if (!skip_to_next)
  {
    beta[i,] = colMeans(results$beta)
    beta_truth[i,] = data$beta_true
    CI = t(colQuantiles(results$beta, probs = c(0.05, 0.95)))
    beta_low[i,] = CI[1,] # CI - lower
    beta_high[i,] = CI[2,] # CI - higher
    beta_bias[i, ] = data$beta - colMeans(results$beta)
    
    
    sigma[i] = mean(results$sigma_2)
    sigma_truth[i] = data$sigma_2_true
    CI = quantile(results$sigma_2, probs = c(0.05, 0.95))
    sigma_low[i] = CI[1]
    sigma_high[i] = CI[2]
    sigma_bias[i] = data$sigma_2_true - mean(results$sigma_2)
    
    # now get g + mu
    mu = rep(colMeans(results$mu), each = t_vals)
    
    g = colMeans(results$g)
    
    
    g_mu[i,] = g + mu
    # CI = quantile(g + mu, probs = c(0.05, 0.95))
    # g_mu_low[i,] = CI[1,]
    # g_mu_high[i,] = CI[2,]
    # TODO fix this
    g_mu_truth[i,] = unlist(data$g_true) + rep(data$mu_true, each = t_vals) 
    g_mu_bias[i,] =  g_mu_truth[i,] - g_mu[i,]
    
    timing[i] = results$time
    
  }
  print(paste(t_vals, "/", i))
}



# which ones to re-do
redo = which(rowMeans(beta) == 0)
if (length(redo) != 0)
{
  print(paste("For t =", t_vals, "there are", length(redo)))
  print(paste(redo, collapse = ","))
}


write_filename = ifelse(USE_NNGP, "outputNNGP/", "output/")
write.csv(beta, row.names = FALSE, file = paste(write_filename, "beta_", t_vals, ".csv", sep = ""))
write.csv(beta_truth, row.names = FALSE, file = paste(write_filename, "beta_truth_", t_vals, ".csv", sep = ""))
write.csv(sigma, row.names = FALSE, file = paste(write_filename, "sigma_", t_vals, ".csv", sep = ""))
write.csv(sigma_truth, row.names = FALSE, file = paste(write_filename, "sigma_truth_", t_vals, ".csv", sep = ""))
write.csv(g_mu, row.names = FALSE, file = paste(write_filename, "gmu_", t_vals, ".csv", sep = ""))
write.csv(g_mu_truth, row.names = FALSE, file = paste(write_filename, "gmu_truth_", t_vals, ".csv", sep = ""))
write.csv(timing, row.names = FALSE, file = paste(write_filename, "timing_", t_vals, ".csv", sep = ""))
