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


plot_all = function(results)
{
  
  # posterior mean
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
  
  # plot
  par(mfrow = c(3,4))
  
  # mu
  mu.range = range(c(data$mu_true, mu_pm))
  plot(data$mu_true, mu_pm, pch = 16, main = expression(mu),
       xlab = 'Truth', ylab = 'Posterior mean', ylim = mu.range)
  abline(0, 1, col = 2)
  
  # beta
  beta.range = range(c(data$beta_true, beta_pm))
  plot((data$beta_true), (beta_pm), pch = 16, main = expression(beta),
       xlab = 'Truth', ylab = 'Posterior mean', ylim = beta.range)
  abline(0, 1, col = 2)
  
  # sigma_2
  sigma_2.range = range(c(data$sigma_2_true, results$sigma_2))
  plot(results$sigma_2, type = 'l', col = 'dodgerblue',
       main = expression(sigma^2),
       xlab = 'MCMC iterations', ylab = 'MCMC samples',
       ylim = sigma_2.range)
  abline(h = sigma_2_pm, col = 1, lwd = 2)
  abline(h = data$sigma_2_true, col = 2, lwd = 2)
  
  # sigma_2B
  sigma_2_B.range = range(c(data$sigmaB_2_true, results$sigmaB_2))
  plot(results$sigmaB_2, type = 'l', col = 'dodgerblue',
       main = paste(expression(sigma^2), "B"),
       xlab = 'MCMC iterations', ylab = 'MCMC samples',
       ylim = sigma_2_B.range)
  abline(h = sigma_2B_pm, col = 1, lwd = 2)
  abline(h = data$sigmaB_2_true, col = 2, lwd = 2)
  
  # xi
  if (length(data$xi_true) == length(xi_pm))
  {
    xi.range = range(c(data$xi_true, xi_pm))
    plot(xi_pm, data$xi_true, pch = 16, main = expression(xi),
         xlab = 'Truth', ylab = 'Posterior mean', ylim = xi.range)
    abline(0, 1, col = 2)
  }
  
  # g
  g.range = range(c(data$g_true, g_pm))
  plot(g_pm, unlist(data$g_true), pch = 16, main = expression(g),
       xlab = 'Truth', ylab = 'Posterior mean', ylim = g.range)
  abline(0, 1, col = 2)
  
  # g(w) vs w
  w.range = range(data$w_true, w_pm)
  g.range = range(c(data$g_true, g_pm))
  plot(unlist(data$w_true), unlist(data$g_true), pch = 16,
       main = 'g(w) vs. w',
       col = 2, xlab = 'w', ylab = 'g(w)', xlim = w.range, ylim = g.range)
  points(w_pm, g_pm, pch = 16)
  
  # g(w) vs w
  w.range = range(data$w_true, w_pm)
  g.range = range(c(data$g_true, g_pm))
  plot(unlist(data$w_true), unlist(data$g_true), pch = 16,
       main = 'g(w) vs. w',
       col = 2, xlab = 'w', ylab = 'g(w)', xlim = w.range, ylim = g.range)
  points(w_pm, g_pm, pch = 16)
  
  # log likelihood
  loglhood.range = range(c(data$loglhood_true, results$loglhood))
  plot(results$loglhood, type = 'l', col = 'dodgerblue',
       main = 'Log likelihood',
       xlab = 'MCMC iterations', ylab = 'Log likelihood',
       ylim = loglhood.range)
  abline(h = loglhood_pm, col = 1, lwd = 2)
  abline(h = data$loglhood_true, col = 2, lwd = 2)
  
  #lk plots
  lk.range = range(c(data$lK_true, results$lK))
  plot(results$lK, type = 'l', col = 'dodgerblue',
       main = expression(l_k),
       xlab = 'MCMC iterations', ylab = 'MCMC samples',
       ylim = lk.range)
  abline(h = mean(lk_pm), col = 1, lwd = 2)
  abline(h = data$lK_true, col = 2, lwd = 2)
  
  #lb plot
  lb.range = range(c(data$lB_true, results$lB))
  plot(results$lB, type = 'l', col = 'dodgerblue',
       main = expression(l_b),
       xlab = 'MCMC iterations', ylab = 'MCMC samples',
       ylim = lb.range)
  abline(h = mean(lb_pm), col = 1, lwd = 2)
  abline(h = data$lB_true, col = 2, lwd = 2)
  
  
  ############
  new_g = matrix(g_pm, ncol = length(data$g_true), byrow = T)
  g = rowMeans(new_g)
  
  new_w = matrix(colMeans(results$w), ncol = length(data$g_true), byrow = T)
  w = rowMeans(new_w)
  
  new_add = g + mu_pm
  
  # now get data stuff
  g_new = vector()
  for (i in 1:length(data$g_true))
  {
    g_new[i] = mean(data$g_true[[i]])
  }
  
  old_add = g_new + data$mu_true
  
  plot(old_add, new_add, xlab = "True g + mu", ylab = "Sampled g + mu", main = "g  + mu")
  
  
  # aa = cbind(data$beta_true, beta_pm)
  # colnames(aa) = c("Actual Beta", "Estimated Beta")
  # rownames(aa) = c("Beta 1", "Beta 2", "Beta 3", "Beta 4", "Beta 5", "Beta 6", "Beta 7 ", "Beta 8", "Beta 9", "Beta 10")
  # View(aa)
  
  # truth
  # g(w) vs w
  # w.range = range(data$w_true, w_pm)
  # g.range = range(c(data$g_true, g_pm))
  # plot(unlist(data$w_true), unlist(data$g_true), pch = 16,
  #      main = 'g(w) vs. w',
  #      col = 2, xlab = 'w', ylab = 'g(w)', xlim = w.range, ylim = g.range)
  # points(w_pm, g_pm, pch = 16)
  
  
  # res = vector()
  # for (i in 1:length(data$w_true))
  # {
  #   res[i] = mean(data$w_true[[i]])
  # }
  # plot(unlist(data$w_true), w_pm)
}
