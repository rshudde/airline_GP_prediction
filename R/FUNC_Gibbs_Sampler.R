# This is a file to start doing paramater estimates for sampling mu
rm(list = ls())
library(invgamma)
library(MASS)
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_woodchan_samples.R')
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_paramater_estimates.R')
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/DATA_generate_simulation.R')
# initial estimates for paramaters outside loop


# get data
data = generate_simulation_data(50, 10)
X = data$X
y = data$y
beta = data$beta

# initialize hyperparamaters
sigma_mu = 2
alpha_mu = 5
knots = seq(0, 1, 0.1)
n = length(knots)
N = nrow(y)
a = 0.1
b = 0.1

# initialize paramaters
alpha_0 = rep(1, ncol(X[[1]]))
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
  M[[i]] = get_matern(l_k_0, rownames(X[[i]]))
  K[[i]] = get_K_i(sigma_2_0, M[[i]])
  V[[i]] = get_V_i(sigma_2_0, M[[i]], K[[i]])
}
# initial iteration outside the loop

# getting beta
beta_i = get_beta(alpha_0, y, mu_0, X, xi_0, knots, N, sigma_2_0, l_k_0, M, K)
beta = rbind(beta, beta_i)

# getting mu
mu_temp = vector()
g = list()
for (i in 1:nrow(y))
{
  g[[i]] = get_g(X[[i]], beta_i, knots, N, xi_0)
  sigma_mu_post_temp = get_sigma_mu_post(sigma_2_0, sigma_mu, V[[i]])
  alpha_mu_post_temp = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post_temp, g[[i]], V[[i]], y[i,])
  mu_temp[i] = get_mu_i(alpha_mu_post_temp, sigma_mu_post_temp)
}
mu = rbind(mu, mu_temp)

# getting sigma 
sigma_2[1] = get_sigma_squared(a, b, y, M, mu, g)

# getting xi 
xi_temp = get_xi(xi_0, y, mu, X, beta_i, knots, N, sigma_2[1], l_k_0, l_b, M, K)
xi = rbind(xi, xi_temp)
  
# get l_k
l_k[1] = get_lk(y, mu, g, sigma_2[1], l_k_0)


# now loop over everything
B = 10000
start = Sys.time()
for (b in 2:B)
{
  # updating M and K
  for (i in 1:nrow(y))
  {
    M[[i]] = get_matern(l_k[b-1], rownames(X[[i]]))
    K[[i]] = get_K_i(sigma_2[b-1], M[[i]])
    V[[i]] = get_V_i(sigma_2[b-1], M[[i]], K[[i]])
  }
  
  # getting beta
  beta_i = get_beta(alpha_0, y, mu[b-1, ], X, xi[b-1, ], knots, N, sigma_2[b-1], l_k[b-1], M, K)
  beta = rbind(beta, beta_i)
  
  # getting mu
  mu_temp = vector()
  g = list()
  for (i in 1:nrow(y))
  {
    g[[i]] = get_g(X[[i]], beta[b, ], knots, N, xi[b-1, ])
    sigma_mu_post_temp = get_sigma_mu_post(sigma_2[b-1], sigma_mu, V[[i]])
    alpha_mu_post_temp = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post_temp, g[[i]], V[[i]], y[i,])
    mu_temp[i] = get_mu_i(alpha_mu_post_temp, abs(sigma_mu_post_temp))
    # print(paste("iteration: ", i, ": ", mu_temp[i]))
  }
  mu = rbind(mu, mu_temp)
  
  # getting sigma 
  sigma_2[b] = get_sigma_squared(a, b, y, M, mu[b, ], g)
  
  # getting xi 
  xi_temp = get_xi(xi[b-1,], y, mu[b, ], X, beta[b, ], knots, N, sigma_2[b], l_k[b-1], l_b, M, K)
  xi = rbind(xi, xi_temp)
  
  # get l_k
  l_k[b] = get_lk(y, mu[b, ], g, sigma_2[b], l_k[b-1])
  
  if (b %% 100 == 0) print(b)
}
print(Sys.time() - start)

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



