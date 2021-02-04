rm(list = ls())
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/PLOTS_Gibbs_sampler.R')

# generate data
data = generate_simulation_data(n_datasets = 500, n_time = 10, 
                                n_covariates = 10, seed = 6)

# run gibbs sampler
results = gibbs_sampler(data_gibbs = data, B = 10000, 
                        # mu_initial = data$mu_true,
                        # beta_initial = data$beta_true,
                        # sigma_2_initial = data$sigma_2_true,
                        xi_initial = runif(length(data$xi_true), -1, 1),
                        # xi_initial = data$xi_true,
                        # lK_initial = data$lK_true,
                        # lB_initial = data$lB_true
)

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

## plot
par(mfrow = c(3,4))

# mu
mu.range = range(c(data$mu_true, mu_pm))
plot(data$mu_true, mu_pm, pch = 16, main = expression(mu),
     xlab = 'Truth', ylab = 'Posterior mean', ylim = mu.range)
abline(0, 1, col = 2)

# beta
beta.range = range(c(data$beta_true, beta_pm))
plot(data$beta_true, beta_pm, pch = 16, main = expression(beta),
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

# log likelihood
loglhood.range = range(c(data$loglhood_true, results$loglhood))
plot(results$loglhood, type = 'l', col = 'dodgerblue',
     main = 'Log likelihood',
     xlab = 'MCMC iterations', ylab = 'Log likelihood',
     ylim = loglhood.range)
abline(h = loglhood_pm, col = 1, lwd = 2)
abline(h = data$loglhood_true, col = 2, lwd = 2)

#lk plot
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


# g(w) vs w
w.range = range(data$w_true, w_pm)
g.range = range(c(data$g_true, g_pm))
plot(unlist(data$w_true), unlist(data$g_true), pch = 16, 
     main = 'g(w) vs. w',
     col = 2, xlab = 'w', ylab = 'g(w)', xlim = w.range, ylim = g.range)
points(w_pm, g_pm, pch = 16)


############
new_g = matrix(g_pm, ncol = 10, byrow = T)
g = rowMeans(new_g)

new_add = g + mu_pm

# now get old
g_new = vector()
for (i in 1:length(data$g_true))
{
    g_new[i] = mean(data$g_true[[i]])
}
old_add = g_new + data$mu_true

par(mfrow = c(1,1))
plot(old_add, new_add, xlab = "True g + mu", ylab = "Sampled g + mu")










