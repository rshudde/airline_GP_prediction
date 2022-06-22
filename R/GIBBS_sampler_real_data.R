rm(list = ls())
library(dplyr)
library(stats)
library(recipes)
library(caret)
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

# load in the data
local = FALSE
if (local)
{
  load("/Users/rachaelshudde/Desktop/X_list.rda")
  load("/Users/rachaelshudde/Desktop/Z_list.rda")
} else {
  load("X_list.rda")
  load("Z_list.rda")
}


# start gibbs sampler
beta_names = colnames(X_list[[1]])
time_idx = apply(Z, 1, function(x) which(!is.na(x)))
start = 1
end = length(X_list)
data_actual = list(y = Z[start:end,], X = X_list[start:end], time_idx = time_idx[start:end])
results = gibbs_sampler_r(data_gibbs = data_actual, 
                          B = 40,
                          xi_initial = runif(ncol(data_actual$X[[1]]), -1, 1),
                          burn_in = 0.5,
                          NNGP = TRUE,
                          nNeighbour = 20,
                          n_to_store = 18000)
# data_gibbs = data_actual; B = 40; n_to_store = 50; runif(ncol(data_actual$X[[1]]), -1, 1)
save(results, file = "REAL_DATA_GIBBS.rda")


##############################
load("/Users/rachaelshudde/Desktop/y.rda")
load("/Users/rachaelshudde/Desktop/m.rda")
load("/Users/rachaelshudde/Desktop/g_gibbs.rda")
load("/Users/rachaelshudde/Desktop/sigma.rda") # ok
load("/Users/rachaelshudde/Desktop/lk.rda")
load("/Users/rachaelshudde/Desktop/time.rda")

lK_post_current = get_lk(y = y, mu = mu_post_current, g = g_gibbs, sigma_2 = sigma_2_post_current, 
                         lk_0 = lK_post_current, time = time_idx)


##############################3 
# get posterior g samples
# g_s = list()
# post_g = colMeans(results$g)
# last = 1
# for (i in 1:end)
# {
#   temp = Z[i,]
#   if (length(which(is.na(temp))) > 0) temp = temp[-which(is.na(temp))]
#   current = length(temp)
#   g_s[[i]] = post_g[last:(last + current - 1)]
#   print(paste("for", i, " - ", last, ":", last + current - 1, "(", length(temp), ")","(", length(last:(last + current - 1)), ")"))
#   last = last + current - 1
# }

# #plots
# par(mfrow = c(3,3))
# plot(colMeans(results$beta), main = "Beta plots")
# plot(results$sigma_2, type = "l", main = "sigma_2 plot")
# plot(results$lB, type  = "l", main = "lb")
# plot(colMeans(results$mu), main = "mu plots")
# plot(results$loglhood, type = "l", main = "loglikelihood")
# plot(colMeans(results$xi), main = "xi results")
# plot(results$beta[,5], type = "l")
# plot(results$beta[,25], type = "l")
# plot(results$beta[,45], type = "l")
# par(mfrow = c(1,1)); plot(colMeans(results$g))
# 
# 
# betas = as.data.frame(beta_names)
# betas$beta_values = colMeans(round(results$beta, 3))
# betas = betas[order(abs(betas$beta_values), decreasing = TRUE),]
# betas
# 
# # look at prediction values of 
# idx = sample(1:end, 1)
# Z2 = Z[idx,]
# Z2 = Z2[-which(is.na(Z2))]
# Z2_pred = colMeans(results$mu)[idx] + g_s[[idx]]
# pred = as.data.frame(cbind(Z2, Z2_pred))
# colnames(pred) = c("actual", "predicted")
# # pred$actual = exp(pred$actual) - adjustment
# # pred$predicted = exp(pred$predicted)
# par(mfrow = c(1,2)); plot(pred$actual, pred$predicted); plot(abs(pred$actual - pred$predicted))
# 

# 
