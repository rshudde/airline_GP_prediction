#setwd("C:\\Users\\anany\\Downloads\\NNGP\\NNGP")
setwd("C:\\Users\\anany\\Desktop\\Research\\Flight_Delay\\output")
library(ggplot2)

#read in the data for beta
beta_20 = read.csv("beta_20.csv")
beta_40 = read.csv("beta_40.csv")
beta_60 = read.csv("beta_60.csv")
beta_80 = read.csv("beta_80.csv")
#beta_100 = read.csv("beta_100.csv")


beta_20_truth = read.csv("beta_truth_20.csv")
beta_40_truth = read.csv("beta_truth_40.csv")
beta_60_truth = read.csv("beta_truth_60.csv")
beta_80_truth = read.csv("beta_truth_80.csv")
#beta_100_truth = read.csv("beta_truth_100.csv")


#calculate the consistancy measure for beta
beta20 = apply(as.matrix(abs(beta_20) - abs(beta_20_truth)),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(beta_20_truth))
beta40 = apply(as.matrix(abs(beta_40) - abs(beta_40_truth)),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(beta_40_truth))
beta60 = apply(as.matrix(abs(beta_60) - abs(beta_60_truth)),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(beta_60_truth))
beta80 = apply(as.matrix(abs(beta_80) - abs(beta_80_truth)),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(beta_80_truth))
#beta100 = apply(as.matrix(abs(beta_100) - abs(beta_100_truth)),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(beta_100_truth))



#make the violin plots for beta
t = rep(c(20, 40, 60, 80), each = 50)
beta_mat = data.frame(c(beta20, beta40, beta60, beta80),t)
colnames(beta_mat) = c("beta", "t")
beta_mat$t = as.factor(beta_mat$t)

#a = c(5,50)
b = c(36, 37, 38, 45, 46)
c = c(1,  4,  5,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
      23, 25, 26, 27, 28, 29, 30, 32, 33, 34, 37, 38, 40, 41, 42, 44, 45, 47, 48, 49, 50)
idx = c(b+50, c+150)

b = ggplot(beta_mat[-idx, ], aes(x = t, y = beta, color = t)) + geom_boxplot(width = 0.5) +ggtitle("Boxplot for beta (normal)");b

#save the plot for beta
png(file="..\\boxplot_beta_varying_t.png", width=600, height=350)
b
dev.off()


#read in the data for gmu
gmu_20 = read.csv("gpmu_20.csv")
gmu_40 = read.csv("gpmu_40.csv")
gmu_60 = read.csv("gpmu_60.csv")
gmu_80 = read.csv("gpmu_80.csv")
#gmu_100 = read.csv("gpmu_100.csv")


gmu_20_truth = read.csv("gmu_truth_20.csv")
gmu_40_truth = read.csv("gmu_truth_40.csv")
gmu_60_truth = read.csv("gmu_truth_60.csv")
gmu_80_truth = read.csv("gmu_truth_80.csv")
#gmu_100_truth = read.csv("gmu_truth_100.csv")

#calculate the consistancy measure for gmu
gmu20 = apply(as.matrix(gmu_20 - gmu_20_truth),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_20_truth))
gmu40 = apply(as.matrix(gmu_40 - gmu_40_truth),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_40_truth))
gmu60 = apply(as.matrix(gmu_60 - gmu_60_truth),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_60_truth))
gmu80 = apply(as.matrix(gmu_80 - gmu_80_truth),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_80_truth))
#gmu100 = apply(as.matrix(gmu_100 - gmu_100_truth),1, function(x){sqrt(sum(x^2))})/sqrt(ncol(gmu_100_truth))

#make the violin plots for gmu
t = rep(c(20, 40, 60, 80), each = 50)
gmu_mat = data.frame(c(gmu20, gmu40, gmu60, gmu80),t)
colnames(gmu_mat) = c("gmu", "t")
gmu_mat$t = as.factor(gmu_mat$t)

g = ggplot(gmu_mat[-idx,], aes(x = t, y = gmu, color = t))+ geom_violin(trim=FALSE) + geom_boxplot(width=0.1)+ ggtitle("Violin plot for g + mu (normal)") + ylab(paste("g + ", expression(mu)));g

#save the plot for gmu
png(file="..\\violin_plot_gmu_normal_varying_t.png", width=600, height=350)
g
dev.off()

#read in the data for sigma
sigma_20 = read.csv("sigma_20.csv")
sigma_40 = read.csv("sigma_40.csv")
sigma_60 = read.csv("sigma_60.csv")
sigma_80 = read.csv("sigma_80.csv")
#sigma_100 = read.csv("sigma_100.csv")


sigma_20_truth = read.csv("sigma_truth_20.csv")
sigma_40_truth = read.csv("sigma_truth_40.csv")
sigma_60_truth = read.csv("sigma_truth_60.csv")
sigma_80_truth = read.csv("sigma_truth_80.csv")
#sigma_100_truth = read.csv("sigma_truth_100.csv")


#sigma = c((sigma_100 - sigma_100_truth)^2, (sigma_1000 - sigma_1000_truth)^2, (sigma_1100 - sigma_1100_truth)^2)
sigma = rbind(sigma_20, sigma_40, sigma_60, sigma_80)
t = rep(c(20, 40, 60, 80), each = 50)
sigma_mat = data.frame(sigma, t)
sigma_mat$t = as.factor(sigma_mat$t)
s = ggplot(sigma_mat[-idx,], aes(x = x, color = t)) + geom_density(trim = TRUE) + ggtitle("Density plot for sigma (normal)") + geom_vline(xintercept = 0.25)+xlab(expression(sigma^2));s

#save the plot for sigma
png(file="..\\density_plot_sigma_normal_varying_t.png", width=600, height=350)
s
dev.off()
