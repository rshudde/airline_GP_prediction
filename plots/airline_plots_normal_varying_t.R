#setwd("C:\\Users\\anany\\Downloads\\NNGP\\NNGP")
library(ggplot2)

#read in the data for beta
NNGP = FALSE
n_reps = 50
betas = list()
betas_truth = list()
gmu = list()
gmu_truth = list()
sigma = list()
sigma_truth = list()

t_vals = c(20, 40, 60, 80, 100)
# t_vals = c(20, 40)
count = 1
for (i in t_vals)
{
  outputfile = ifelse(NNGP, "outputNNGP", "output")

  betas[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/beta_", i, ".csv", sep = ""))
  betas_truth[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/beta_truth_", i, ".csv", sep = ""))
  
  gmu[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/gmu_", i, ".csv", sep = ""))
  gmu_truth[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/gmu_truth_", i, ".csv", sep = ""))
  
  sigma[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/sigma_", i, ".csv", sep = ""))
  sigma_truth[[count]] = read.csv(paste("/Users/rachaelshudde/Desktop/", outputfile, "/sigma_truth_", i, ".csv", sep = ""))

  # standardize the results
  standardized_beta = apply(as.matrix(abs(betas[[count]]) - abs(betas_truth[[count]])),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(betas_truth[[count]]))

  standardized_gmu = apply(as.matrix(gmu[[count]] - gmu_truth[[count]]),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_truth[[count]]))
  betas[[count]] = standardized_beta
  gmu[[count]] = standardized_gmu
  count = count + 1
}

#make the violin plots for beta
t_vals_plots = rep(t_vals, each = n_reps)
beta_mat = data.frame(unlist(betas),t_vals_plots) # column of betas, column of t_values
gmu_mat = data.frame(unlist(gmu),t_vals_plots) 
sigma_mat = data.frame(unlist(sigma), t_vals_plots)

colnames(beta_mat) = c("beta", "t")
colnames(gmu_mat) = c("gmu", "t")
colnames(sigma_mat) = c("sigma", "t")

beta_mat$t = as.factor(beta_mat$t)
gmu_mat$t = as.factor(gmu_mat$t)
sigma_mat$t = as.factor(sigma_mat$t)


## actually creating plots
beta_plot = ggplot(beta_mat, aes(x = t, y = beta, color = t_vals_plots)) + geom_boxplot(width = 0.5) + ggtitle("Boxplot for beta (normal)")

g_plot = ggplot(gmu_mat, aes(x = t, y = gmu, color = t_vals_plots)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + 
  ggtitle("Violin plot for g + mu (normal)") + ylab(paste("g + ", expression(mu)))

sigma_plot = ggplot(sigma_mat, aes(x = sigma, color = t)) + geom_density(trim = TRUE) + ggtitle("Density plot for sigma (normal)") + 
  geom_vline(xintercept = 0.25) + xlab(expression(sigma^2))

# display plots
beta_plot
g_plot
sigma_plot




