#setwd("C:\\Users\\anany\\Downloads\\NNGP\\NNGP")
library(ggplot2)

#read in the data for beta
t_vals = 10
count = 1
NNGP = TRUE
n_reps = 80
betas = list()
betas_truth = list()
gmu = list()
gmu_truth = list()
sigma = list()
sigma_truth = list()
time = vector(length = length(t_vals))


path = "/Users/rachaelshudde/Desktop/"
# path = "C:/Users/anany/Desktop/Research/Flight_Delay/test_t/"

for (i in t_vals)
{
  outputfile = ifelse(NNGP, "outputNNGP", "output")
 
  betas[[count]] = read.csv(paste(path, outputfile, "/beta_", i, ".csv", sep = ""))
  betas_truth[[count]] = read.csv(paste(path, outputfile, "/beta_truth_", i, ".csv", sep = ""))
  
  gmu[[count]] = read.csv(paste(path, outputfile, "/gmu_", i, ".csv", sep = ""))
  gmu_truth[[count]] = read.csv(paste(path, outputfile, "/gmu_truth_", i, ".csv", sep = ""))
  
  sigma[[count]] = read.csv(paste(path, outputfile, "/sigma_", i, ".csv", sep = ""))
  sigma_truth[[count]] = read.csv(paste(path, outputfile, "/sigma_truth_", i, ".csv", sep = ""))
  
  time[count] = read.csv(paste(path, outputfile, "/timing", i, ".csv", sep = ""))

  # standardize the results
  standardized_beta = apply(as.matrix(abs(betas[[count]]) - abs(betas_truth[[count]])),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(betas_truth[[count]]))

  standardized_gmu = apply(as.matrix(gmu[[count]] - gmu_truth[[count]]),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_truth[[count]]))
  betas[[count]] = standardized_beta
  gmu[[count]] = standardized_gmu
  count = count + 1
}

#make the violin plots for beta
t_vals_plots = as.character(rep(t_vals, each = n_reps))
beta_mat = data.frame(unlist(betas),t_vals_plots) # column of betas, column of t_values
gmu_mat = data.frame(unlist(gmu),t_vals_plots) 
sigma_mat = data.frame(unlist(sigma), t_vals_plots)

colnames(beta_mat) = c("beta", "t")
colnames(gmu_mat) = c("gmu", "t")
colnames(sigma_mat) = c("sigma", "t")

beta_mat$t = factor(beta_mat$t, levels = unique(gmu_mat$t))
gmu_mat$t = factor(gmu_mat$t, levels = unique(gmu_mat$t))
sigma_mat$t = factor(sigma_mat$t, levels = unique(gmu_mat$t))

# #create lable
# lable = ifelse(NNGP,"NNGP", "Full")


## actually creating plots

if(NNGP == TRUE){
  beta_NNGP = ggplot(beta_mat, aes(x = t, y = beta, color = t)) + geom_boxplot(width = 0.5) +
    ggtitle(expression("Consistancy of function estimation for"~beta~"(NNGP Sampler)")) +
    xlab("Time points") + ylab ("Norm Distance") + 
    guides(color=guide_legend(title="Time points"))+ theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))   #ylab(expression("||"~beta~"||"^{2}))
  
  gmu_NNGP = ggplot(gmu_mat, aes(x = t, y = gmu, color = t)) + geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1) + ggtitle(expression("Consistancy of function estimation for g +"~mu~"(NNGP Sampler)")) +
    xlab("Time points") + ylab ("Norm Distance")+
    guides(color=guide_legend(title="Time points"))+ theme_bw() +theme(plot.title = element_text(hjust = 0.5))  #ylab(expression("||"~g + mu~"||"^{2}))
  
  sigma_NNGP = ggplot(sigma_mat, aes(x = sigma, color = t)) + geom_density(trim = TRUE) +
    ggtitle(expression("Density plot for"~sigma^{2}~"(NNGP Sampler)")) +  
    geom_vline(xintercept = 0.25) + xlab("Time points") +  
    ylab(expression(Density~of~sigma^{2}))+ guides(color=guide_legend(title="Time points")) + theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) 
  # display plots
  beta_NNGP
  g_NNGP
  sigma_NNGP
  
  
}else{
  beta_full = ggplot(beta_mat, aes(x = t, y = beta, color = t)) + geom_boxplot(width = 0.5) +
    ggtitle(expression("Consistancy of function estimation for"~beta~"(Full Sampler)")) +
    xlab("Time points") + ylab ("Norm Distance") + 
    guides(color=guide_legend(title="Time points"))+ theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))   #ylab(expression("||"~beta~"||"^{2}))
  
  gmu_full = ggplot(gmu_mat, aes(x = t, y = gmu, color = t)) + geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1) + ggtitle(expression("Consistancy of function estimation for g +"~mu~"(Full Sampler)")) +
    xlab("Time points") + ylab ("Norm Distance")+
    guides(color=guide_legend(title="Time points"))+ theme_bw() +theme(plot.title = element_text(hjust = 0.5))  #ylab(expression("||"~g + mu~"||"^{2}))
  
  sigma_full = ggplot(sigma_mat, aes(x = sigma, color = t)) + geom_density(trim = TRUE) +
    ggtitle(expression("Density plot for"~sigma^{2}~"(Full Sampler)")) +  
    geom_vline(xintercept = 0.25) + xlab("Time points") +  
    ylab(expression(Density~of~sigma^{2}))+ guides(color=guide_legend(title="Time points")) + theme_bw()+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  
  # display plots
  beta_full
  g_full
  sigma_full
  
}





