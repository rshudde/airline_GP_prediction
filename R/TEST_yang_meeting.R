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
Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")


data = generate_simulation_data(n_datasets = 100, n_time = 80, n_covariates = 15, seed = 1, seed2 = 1)
print("GOT DATA")
results = gibbs_sampler_r(data_gibbs = data, B = 100000,
                          xi_initial = runif(length(data$w_true[[1]]), -1, 1),
                          burn_in = 0.5,
                          NNGP = FALSE,
                          n_to_store = 20000)
save(results, file = "TESTRESULTS100.rda")



### 
# load("/Users/rachaelshudde/Desktop/TESTRESULTS.rda")
# plot(results$LOGTEST, main = "Plot of LogLikelihood", xlab = "Iterations", ylab = "Log likelihood", type = "l")
# par(mfrow = c(2,2))
# plot(results$beta[,2], type = "l")
# plot(results$beta[,4], type = "l")
# plot(results$beta[,5], type = "l")
# plot(results$beta[,6], type = "l")
# par(mfrow = c(2,2))
# plot(results$g[,2], type = "l")
# plot(results$g[,4], type = "l")
# plot(results$g[,5], type = "l")
# plot(results$g[,6], type = "l")
# par(mfrow = c(2,2))
# plot(results$mu[,2], type = "l")
# plot(results$mu[,4], type = "l")
# plot(results$mu[,5], type = "l")
# plot(results$mu[,6], type = "l")
# par(mfrow = c(2,2))
# plot(results$w[,2], type = "l")
# plot(results$w[,4], type = "l")
# plot(results$w[,5], type = "l")
# plot(results$w[,6], type = "l")
# par(mfrow = c(2,2))
# plot(results$xi[,2], type = "l")
# plot(results$xi[,4], type = "l")
# plot(results$xi[,5], type = "l")
# plot(results$xi[,6], type = "l")
# par(mfrow = c(1,1))
# plot(results$sigma_2, type = "l")
# plot(results$sigmaB_2, type = "l")
# plot(results$lB, type = "l")
# plot(results$lK, type = "l")


# new = results$g + results$w
# sample1 = sample(1:ncol(new), 4, replace = F)
# sample1 = c(6705, 3891, 2154, 855)
# par(mfrow = c(2,2))
# plot(new[, sample1[1]], type = "l", xlab = "Iterations (burnin removed, thinning applied)", ylab = "", main = "plot of g + w")
# plot(new[, sample1[2]], type = "l", xlab = "Iterations (burnin removed, thinning applied)", ylab = "", main = "plot of g + w")
# plot(new[, sample1[3]], type = "l", xlab = "Iterations (burnin removed, thinning applied)", ylab = "", main = "plot of g + w")
# plot(new[, sample1[4]], type = "l", xlab = "Iterations (burnin removed, thinning applied)", ylab = "", main = "plot of g + w")







