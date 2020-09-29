rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')

# initial estimates for paramaters outside loop
# 
# knots = seq(0, 1, length.out = 10) # suggested to be 
# 
# 
# # get data
# n_covariates = 10
# n_datasets = 100
# data = generate_simulation_data(n_datasets = n_datasets, n_covariates = n_covariates, knots = knots)
# 
# results = gibbs_sampler(data_gibbs = data, knots_gibbs = knots, B = 500)

###################
beta_true = data$beta
beta_post = results$beta
beta_estimates = apply(beta_post, 2, function(x) quantile(x, c(0.025, 0.975)))
beta_means = colMeans(beta_post)

y_grid = c(1:nrow(beta_post))
par(mfrow = c(1,1))
for (i in 1:length(beta_true))
{
  max_value = max(max(beta_post[, i]), beta_true[i], beta_means[i])
  min_value = min(min(beta_post[, i]), beta_true[i], beta_means[i])
  plot(y_grid, beta_post[, i], type = "l", ylim = c(-1.5*abs(min_value), 1.5*max_value), main = paste("Plot of beta[, ", i, "]"), xlab = "MCMC iteration",
       ylab = "beta estimates", xlim = c(10, nrow(beta_post) - 10))
  abline(h = beta_true[i], col = "red")
  abline(h = beta_means[i], col = "blue")
  
  # adding confience intervals
  high = beta_estimates[2, i]
  low = beta_estimates[1, i]
  
  abline(h = high, col = "darkgreen", lty = 3)
  abline(h = low, col = "darkgreen", lty = 3)
  # polygon(c(y_grid, rev(y_grid)), c(rep(high, length(y_grid)), rev(rep(low, length(y_grid)))),  
  #         col = yarrr::transparent("darkgreen", trans.val = .9), border = NA)
  
  legend("topleft", legend=c("Actual Beta Value", "Posterior Mean", "95% CI"), col = c("red", "blue", "darkgreen"), 
         lty=c(1, 1, 3), cex = 0.6)
}



# plot(beta_actual, beta_means, type = "l")


