###################

plot_beta = function(data, results)
{
  beta_true = data$beta
  beta_post = results$beta
  beta_estimates = apply(beta_post, 2, function(x) quantile(x, c(0.025, 0.975)))
  beta_means = colMeans(beta_post)
  
  y_grid = c(1:nrow(beta_post))
  
  
  ## beta plots
  
  par(mfrow = c(2,3))
  for (i in 1:length(beta_true))
  {
    max_value = max(max(beta_post[, i]), beta_true[i], beta_means[i])
    min_value = min(min(beta_post[, i]), beta_true[i], beta_means[i])
    plot(y_grid, beta_post[, i], type = "l", ylim = c(min_value, 1*max_value), main = paste("Plot of beta[, ", i, "]"), xlab = "MCMC iteration",
         ylab = "beta estimates", xlim = c(10, nrow(beta_post) - 10), col = "gray")
    abline(h = beta_true[i], col = "darkgreen", lwd = 2)
    abline(h = beta_means[i], col = "red", lty = 2, lwd = 2)
    
    # adding confience intervals
    high = beta_estimates[2, i]
    low = beta_estimates[1, i]
    
    abline(h = high, col = "darkgreen", lty = 3)
    abline(h = low, col = "darkgreen", lty = 3)
    
    # legend("topleft", legend=c("Actual Beta Value", "Posterior Mean", "95% CI"), col = c("red", "blue", "darkgreen"), 
    #        lty=c(1, 1, 3), cex = 0.4, inset=c(-0.2,0))
  }
}

###################
plot_sigma = function(data, results)
{
  sigma_estimates = results$sigma_2
  sigma_actual = data$sigma_2
  par(mfrow = c(1,1))
  y_grid = c(1:length(sigma_estimates))
  
  plot(y_grid, sigma_estimates, type = "l", main = "Plot of sigma2", col = "gray")
  abline(h = sigma_actual, col = "darkgreen", lwd = 2)
  abline(h = mean(sigma_estimates), col = "red", lty = 2, lwd = 2)
}


###################
plot_lk = function(data, results)
{
  lk_estimates = results$lk
  lb_estimates = results$lb
  
  lk_actual = data$l_k
  max_value_lk = max(max(lk_estimates), lk_actual)
  par(mfrow = c(1,1))
  y_grid = c(1:length(lk_estimates))
  
  par(mfrow = c(1,2))
  
  plot(y_grid, lk_estimates, type = "l", main = "Plot of lk", col = "gray", ylim = c(min(lk_estimates), 1.3*max_value_lk))
  abline(h = lk_actual, col = "darkgreen", lwd = 2)
  abline(h = mean(lk_estimates), col = "red", lty = 2, lwd = 2)
  
  plot(y_grid, lb_estimates, type = "l", main = "Plot of lb", col = "gray", ylim = c(min(lb_estimates), max(lb_estimates)))
}


###################
plot_mu = function(data, results)
{
  mu_true = data$mu
  mu_post = results$mu
  mu_estimates = apply(mu_post, 2, function(x) quantile(x, c(0.025, 0.975)))
  mu_means = colMeans(mu_post)
  y_grid = c(1:nrow(mu_post))
  
  ## mu plots
  samples = sample(1:ncol(mu_post), 6, replace = FALSE)
  
  par(mfrow = c(2,3))
  for (i in samples)
  {
    max_value = max(max(mu_post[, i]), mu_true[i], mu_means[i])
    min_value = min(min(mu_post[, i]), mu_true[i], mu_means[i])
    plot(y_grid, mu_post[, i], type = "l", ylim = c(-1*abs(min_value), 1*max_value), main = paste("Plot of mu[, ", i, "]"), xlab = "MCMC iteration",
         ylab = "mu estimates", xlim = c(10, nrow(mu_post) - 10), col = "gray")
    abline(h = mu_true[i], col = "darkgreen", lwd = 2)
    abline(h = mu_means[i], col = "red", lty = 2, lwd = 2)
    
    # adding confience intervals
    high = mu_estimates[2, i]
    low = mu_estimates[1, i]
    
    abline(h = high, col = "darkgreen", lty = 3)
    abline(h = low, col = "darkgreen", lty = 3)
    # polygon(c(y_grid, rev(y_grid)), c(rep(high, length(y_grid)), rev(rep(low, length(y_grid)))),  
    #         col = yarrr::transparent("darkgreen", trans.val = .9), border = NA)
    # 
    # legend("topleft", legend=c("Actual mu Value", "Posterior Mean", "95% CI"), col = c("red", "blue", "darkgreen"), 
    #        lty=c(1, 1, 3), cex = 0.4, inset=c(-0.2,0))
  }
  
}


###################
plot_xi = function(data, results)
{
  xi_true = data$xi
  xi_post = results$xi
  xi_estimates = apply(xi_post, 2, function(x) quantile(x, c(0.025, 0.975)))
  xi_means = colMeans(xi_post)
  y_grid = c(1:nrow(xi_post))
  
  ## xi plots
  
  samples = sample(1:ncol(xi_post), 6, replace = FALSE)
  
  par(mfrow = c(2,3))
  for (i in samples)
  {
    max_value = max(max(xi_post[, i]), xi_true[i], xi_means[i])
    min_value = min(min(xi_post[, i]), xi_true[i], xi_means[i])
    plot(y_grid, xi_post[, i], type = "l", ylim = c(-1*abs(min_value), 1*max_value), main = paste("Plot of xi[, ", i, "]"), xlab = "MCMC iteration",
         ylab = "xi estimates", xlim = c(10, nrow(xi_post) - 10), col = "gray")
    abline(h = xi_true[i], col = "darkgreen", lwd = 2)
    abline(h = xi_means[i], col = "red", lty = 2, lwd = 2)
    
    # adding confience intervals
    high = xi_estimates[2, i]
    low = xi_estimates[1, i]
    
    abline(h = high, col = "darkgreen", lty = 3)
    abline(h = low, col = "darkgreen", lty = 3)
    # polygon(c(y_grid, rev(y_grid)), c(rep(high, length(y_grid)), rev(rep(low, length(y_grid)))),  
    #         col = yarrr::transparent("darkgreen", trans.val = .9), border = NA)
    # 
    # legend("topleft", legend=c("Actual xi Value", "Posterior Mean", "95% CI"), col = c("red", "blue", "darkgreen"), 
    #        lty=c(1, 1, 3), cex = 0.4, inset=c(-0.2,0))
  }
}



