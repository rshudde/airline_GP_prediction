burn_in = floor(B*.2)

beta_post = beta_gibbs[-c(1:burn_in), ]
beta_actual = data$beta
beta_estimates = apply(beta_post, 2, function(x) quantile(x, c(0.025, 0.975)))
beta_means = colMeans(beta_post)

par(mfrow = c(1,1))
for (i in 1:length(beta_true))
{
  max_value = max(max(beta_post[, i]), beta_true[i], beta_means[i])
  min_value = min(min(beta_post[, i]), beta_true[i], beta_means[i])
  plot(beta_post[, i], type = "l", ylim = c(-1.5*abs(min_value), 1.5*max_value), main = paste("Plot of beta[, ", i, "]"), xlab = "MCMC iteration",
       ylab = "beta estimates" )
  abline(h = beta_true[i], col = "red")
  abline(h = beta_means[i], col = "blue")
  # abline(h = beta_estimates[1, i], col = "lightgreen")
  # abline(h = beta_estimates[2, i], col = "lightgreen")
  legend("topleft", legend=c("Actual Beta Value", "Posterior Mean"), col = c("red", "blue"), lty=(rep(1,2)), cex = 0.6)
}



# plot(beta_actual, beta_means, type = "l")


