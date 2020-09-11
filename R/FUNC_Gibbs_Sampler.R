# This is a file to start doing paramater estimates for sampling mu
rm(list = ls())
library(invgamma)
library(MASS)
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_woodchan_samples.R')
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_paramater_estimates.R')

# initial estimates for paramaters outside loop


# get data
nrow = 50
ncol = 20
data = list()
y = vector()

n_datasets = 40
for (i in 1:n_datasets)
{
  mean = sample(1:10, 1)
  sd = sample(2:9, 1)
  
  temp_x = normalize_data(matrix(rnorm(nrow*ncol, mean, sd), nrow, ncol))
  
  min_y = sample(1:1000, 1)
  max_y = sample(1:1000, 1)
  temp_y = sample(min_y:max_y, nrow, replace = TRUE)
  
  data[[i]] = temp_x
  y = rbind(y, temp_y)
}


# initialize hyperparamaters
sigma_mu = 2
alpha_mu = 5
knots = seq(0, 1, 0.1)
n = length(knots)
N = nrow(y)
a = 0.1
b = 0.1

# initialize paramaters
alpha_0 = rep(1, ncol(data[[1]]))
xi_0 = rnorm(n) # length of knots
mu_0 = rep(0, nrow(y)) # length is number of flights
sigma_2_0 = 1
l_k_0 = 1
l_b = 1

# empty things
beta = vector() # will be a matrix
sigma_2 = vector()
l_k = vector()
xi = vector() # will be a matrix
mu = vector()

# initialize M and K
M = list()
K = list()
V = list()
g = list()

for (i in 1:nrow(y))
{
  M[[i]] = get_matern(l_k_0, y[i, ])
  K[[i]] = get_K_i(sigma_2_0, M[[i]])
  V[[i]] = get_V_i(sigma_2_0, M[[i]], K[[i]])
}
# initial iteration outside the loop

# getting beta
beta_i = get_beta(alpha_0, y, mu_0, data, xi_0, knots, N, sigma_2_0, l_k, M, K)
beta = rbind(beta, beta_i)

# getting mu
mu_temp = vector()
g = list()
for (i in 1:nrow(y))
{
  g[[i]] = get_g(data[[i]], beta_i, knots, N, xi_0)
  sigma_mu_post_temp = get_sigma_mu_post(sigma_2_0, sigma_mu, V[[i]])
  alpha_mu_post_temp = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post_temp, g[[i]], V[[i]], y[i,])
  mu_temp[i] = get_mu_i(alpha_mu_post_temp, sigma_mu_post_temp)
}
mu = rbind(mu, mu_temp)

# getting sigma 
sigma_2[1] = get_sigma_squared(a, b, y, M, mu, g)

# getting xi 
xi_temp = get_xi(xi_0, y, mu, data, beta_i, knots, N, sigma_2[1], l_k_0, l_b, M, K)
xi = rbind(xi, xi_temp)
  
# get l_k
l_k[1] = get_lk(y, mu, g, sigma_2[1], l_k_0)


# now loop over everything
B = 500
for (b in 2:B)
{
  # updating M and K
  for (i in 1:nrow(y))
  {
    M[[i]] = get_matern(l_k[b-1], y[i, ])
    K[[i]] = get_K_i(sigma_2[b-1], M[[i]])
    V[[i]] = get_V_i(sigma_2[b-1], M[[i]], K[[i]])
  }
  
  # getting beta
  beta_i = get_beta(alpha_0, y, mu[b-1, ], data, xi[b-1, ], knots, N, sigma_2[b-1], l_k[b-1], M, K)
  beta = rbind(beta, beta_i)
  
  # getting mu
  mu_temp = vector()
  g = list()
  for (i in 1:nrow(y))
  {
    g[[i]] = get_g(data[[i]], beta[b, ], knots, N, xi[b-1, ])
    sigma_mu_post_temp = get_sigma_mu_post(sigma_2[b-1], sigma_mu, V[[i]])
    alpha_mu_post_temp = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post_temp, g[[i]], V[[i]], y[i,])
    mu_temp[i] = get_mu_i(alpha_mu_post_temp, abs(sigma_mu_post_temp))
    # print(paste("iteration: ", i, ": ", mu_temp[i]))
  }
  mu = rbind(mu, mu_temp)
  
  # getting sigma 
  sigma_2[b] = get_sigma_squared(a, b, y, M, mu[b, ], g)
  
  # getting xi 
  xi_temp = get_xi(xi[b-1,], y, mu[b, ], data, beta[b, ], knots, N, sigma_2[b], l_k[b-1], l_b, M, K)
  xi = rbind(xi, xi_temp)
  
  # get l_k
  l_k[b] = get_lk(y, mu[b, ], g, sigma_2[b], l_k[b-1])
  
  if (b %% 100 == 0) print(b)
}

burn_in = floor(B*.1)
print(colMeans(beta[-c(1:burn_in), ]))

beta_post = beta[-c(1:burn_in), ]

size = 2
par(mfrow = c(size, size))
for (i in 1:(size*size))
{
  plot(beta_post[, i], type = "l", main = paste("plot of beta[, ", i, "]"))
  abline(h = mean(beta_post[,i]), col = "red")
}



