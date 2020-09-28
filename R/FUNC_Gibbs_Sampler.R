# This is a file to start doing paramater estimates for sampling mu
rm(list = ls())
library(invgamma)
library(MASS)
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_woodchan_samples.R')
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_paramater_estimates.R')
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/DATA_generate_simulation.R')
# initial estimates for paramaters outside loop

knots_gibbs = seq(0, 1, length.out = 200) # suggested to be 


# get data
n_covariates = 16
n_datasets = 300
data = generate_simulation_data(n_datasets = n_datasets, n_covariates = n_covariates, knots = knots_gibbs)
X = data$X
y = data$y
beta_true = data$beta


################################################################################################################
################################################################################################################

# number of iterations for gibbs sampler
B = 1000

# initialize hyperparamaters
sigma_mu_gibbs = 2
alpha_mu_gibbs = 5
n_gibbs = length(knots_gibbs)
N_gibbs = nrow(y)
a_gibbs = 0.1
b_gibbs = 0.1

# initialize paramaters
alpha_0_gibbs = rep(1, ncol(X[[1]]))
xi_0_gibbs = rnorm(n_gibbs) # length of knots
mu_0_gibbs = rep(0, nrow(y)) # length is number of flights
sigma_2_0_gibbs = 1
l_k_0_gibbs = 1
l_b_0_gibbs = 1

# empty things
# beta = vector() # will be a matrix
beta_gibbs = matrix(rep(0, B * n_covariates), nrow = B, ncol = n_covariates)
alpha_gibbs = matrix(rep(0, B * n_covariates), nrow = B, ncol = n_covariates)

sigma_2_gibbs = vector()
lk_gibbs = vector()
# xi = vector() # will be a matrix
xi_gibbs = matrix(rep(0, B * length(knots_gibbs)), nrow = B, ncol = length(knots_gibbs))
# mu = vector()
mu_gibbs = matrix(rep(0, B * n_datasets), nrow = B, ncol = n_datasets)
lb_gibbs = vector()

# empty things - matrix form


# initialize M and K
M_gibbs = list()
K_gibbs = list()
V_gibbs = list()
g_gibbs = list()

for (i in 1:nrow(y))
{
  M_gibbs[[i]] = get_matern(l_k_0_gibbs, rownames(X[[i]]))
  K_gibbs[[i]] = get_K_i(sigma_2_0_gibbs, M_gibbs[[i]])
  V_gibbs[[i]] = get_V_i(sigma_2_0_gibbs, M_gibbs[[i]], K_gibbs[[i]])
}
# initial iteration outside the loop

# getting beta
alpha_gibbs[1, ] =  get_alpha(alpha_0_gibbs, y, mu_0_gibbs, X, xi_0_gibbs, knots_gibbs, N_gibbs, sigma_2_0_gibbs, l_k_0_gibbs, M_gibbs, K_gibbs)
beta_gibbs[1, ] = alpha_gibbs[1, ] / sqrt(sum(alpha_gibbs[1, ]^2))

# getting mu
mu_temp = vector()
g_gibbs = list()
for (i in 1:nrow(y))
{
  g_gibbs[[i]] = get_g(X[[i]], beta_gibbs[1, ], knots_gibbs, N_gibbs, xi_0_gibbs)
  sigma_mu_post_temp = get_sigma_mu_post(sigma_2_0_gibbs, sigma_mu_gibbs, V_gibbs[[i]])
  alpha_mu_post_temp = get_alpha_mu_post(alpha_mu_gibbs, sigma_mu_gibbs, sigma_mu_post_temp, g_gibbs[[i]], V_gibbs[[i]], y[i,])
  mu_temp[i] = get_mu_i(alpha_mu_post_temp, sigma_mu_post_temp)
}
mu_gibbs[1, ] = mu_temp

# getting sigma 
sigma_2_gibbs[1] = get_sigma_squared(a_gibbs, b_gibbs, y, M_gibbs, mu_gibbs, g_gibbs)

# getting xi 
xi_gibbs[1, ] = get_xi(xi_0_gibbs, y, mu_gibbs, X, beta_gibbs[1, ], knots_gibbs, N_gibbs, sigma_2_gibbs[1], l_k_0_gibbs, l_b_0_gibbs, M_gibbs, K_gibbs)
  
# get l_k
lk_gibbs[1] = get_lk(y, mu_gibbs[1, ], g_gibbs, sigma_2_gibbs[1], l_k_0_gibbs)
blah = get_lb(y, l_b_0_gibbs, xi_gibbs[1, ])
if(is.logical(blah)) stop("something worng with lb_gibbs")
lb_gibbs[1] = blah
print(blah)
################################################################################################################
################################################################################################################


# now loop over everything
start = Sys.time()
for (idx in 2:B)
{
  start_inner = Sys.time()
  # updating M and K
  for (i in 1:nrow(y))
  {
    M_gibbs[[i]] = get_matern(lk_gibbs[idx-1], rownames(X[[i]]))
    K_gibbs[[i]] = get_K_i(sigma_2_gibbs[idx-1], M_gibbs[[i]])
    V_gibbs[[i]] = get_V_i(sigma_2_gibbs[idx-1], M_gibbs[[i]], K_gibbs[[i]])
  }

  # getting beta
  alpha_gibbs[idx, ] = get_alpha(alpha_gibbs[idx - 1, ], y, mu_gibbs[idx-1, ], X, xi_gibbs[idx-1, ], knots_gibbs, N_gibbs,
                               sigma_2_gibbs[idx-1], lk_gibbs[idx-1], M_gibbs, K_gibbs)
  beta_gibbs[idx, ] = alpha_gibbs[idx, ] / sqrt(sum(alpha_gibbs[idx, ]^2))
  
  # getting mu
  mu_temp = vector()
  g_gibbs = list()
  for (i in 1:nrow(y))
  {
    g_gibbs[[i]] = get_g(X[[i]], beta_gibbs[idx, ], knots_gibbs, N_gibbs, xi_gibbs[idx-1, ])
    sigma_mu_post_temp = get_sigma_mu_post(sigma_2_gibbs[idx-1], sigma_mu_gibbs, V_gibbs[[i]])
    alpha_mu_post_temp = get_alpha_mu_post(alpha_mu_gibbs, sigma_mu_gibbs, sigma_mu_post_temp,
                                           g_gibbs[[i]], V_gibbs[[i]], y[i,])
    mu_temp[i] = get_mu_i(alpha_mu_post_temp, abs(sigma_mu_post_temp))
    # print(paste("iteration: ", i, ": ", mu_temp[i]))
  }
  mu_gibbs[idx, ] = mu_temp

  # if (idx > 860) stop("iteration 860")
  # getting sigma
  sigma_2_gibbs[idx] = get_sigma_squared(a_gibbs, b_gibbs, y, M_gibbs, mu_gibbs[idx, ], g_gibbs)

  # getting xi
  # xi_0, y, mu, data, beta, knots, N, sigma_2, l_k, l_b, M, K)
  xi_gibbs[idx, ] = get_xi(xi_gibbs[idx-1,], y, mu_gibbs[idx, ], X, beta_gibbs[idx, ], knots_gibbs, N_gibbs,
                     sigma_2_gibbs[idx], lk_gibbs[idx-1], lb_gibbs[idx-1], M_gibbs, K_gibbs)

  # # get l_k and l_b
  lk_gibbs[idx] = get_lk(y, mu_gibbs, g_gibbs, sigma_2_gibbs[idx], lk_gibbs[idx-1]) # should just be passing it mu
  lb_gibbs[idx] = get_lb(y, lb_gibbs[idx-1], xi_gibbs[idx, ])
  
  # lb_gibbs[idx] = 2


  if (idx %% 50 == 0) print(paste("iteration:", idx, "in", round(Sys.time() - start_inner, 2)))
}
print(round(Sys.time() - start),2)






