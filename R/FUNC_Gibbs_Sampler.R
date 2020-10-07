# This is a file to start doing paramater estimates for sampling mu
library(invgamma)
library(MASS)

# functin to do Gibbs sampling 
gibbs_sampler = function(data_gibbs, knots_gibbs, B = 1000,
                         lk_0 = 0.1, lb_0 = 0.1, a = 10^(-3), b = 10^(-3), sigma_mu = 100, alpha_mu = 0,
                         burn_in = 0.3, write = FALSE)
{
  set.seed(1)

  ## stuff for debugging
  B = 1000
  lk_0 = 0.1
  lb_0 = 0.1
  a = 10^(-3)
  b = 10^(-3)
  sigma_mu = 100 # TODO fix naming here
  alpha_mu = 0
  burn_in = 0.3
  write = FALSE
  data_gibbs = data
  knots_gibbs = knots
  
  # get X and y values from the data
  X = data_gibbs$X
  y = data_gibbs$y
  
  # get the number of datasets and covariates
  n_datasets = length(X)
  n_covariates = ncol(X[[1]])
  
  # initialize hyperparamaters
  sigma_mu_gibbs = sigma_mu
  alpha_mu_gibbs = alpha_mu
  n_gibbs = length(knots_gibbs)
  N_gibbs = nrow(y)
  a_gibbs = a
  b_gibbs = b

  # initialize beta / sigma / lk / xi / mu / lb values
  beta_gibbs = matrix(rep(0, B * n_covariates), nrow = B, ncol = n_covariates)
  alpha_gibbs = matrix(rep(0, B * n_covariates), nrow = B, ncol = n_covariates)
  
  sigma_2_gibbs = as.vector(1)
  lk_gibbs = as.vector(1)
  xi_gibbs = matrix(rep(0, B * length(knots_gibbs)), nrow = B, ncol = length(knots_gibbs))
  mu_gibbs = matrix(rep(0, B * n_datasets), nrow = B, ncol = n_datasets)

  xi_gibbs[1, ] = rnorm(n_gibbs, 1, 4) # length of knots
  alpha_gibbs[1, ] = rep(1, n_covariates)
  alpha_0_gibbs = rep(1, n_covariates)
  beta_gibbs[1, ] = alpha_0_gibbs / sqrt(sum(alpha_0_gibbs^2))
  
  lb_gibbs = as.vector(lb_0)
  lk_gibbs = as.vector(lk_0)
  ################################################################################################################
  ################################################################################################################
  
  # now loop over everything
  start = Sys.time()
  for (idx in 2:B)
  {
    start_inner = Sys.time()
    
    M_gibbs = list()
    K_gibbs = list()
    V_gibbs = list()
    
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
    
    # getting sigma
    sigma_2_gibbs[idx] = get_sigma_squared(a_gibbs, b_gibbs, y, M_gibbs, mu_gibbs[idx, ], g_gibbs)
    
    # getting xi
    xi_gibbs[idx, ] = get_xi(xi_gibbs[idx-1,], y, mu_gibbs[idx, ], X, beta_gibbs[idx, ], knots_gibbs, N_gibbs,
                             sigma_2_gibbs[idx], lk_gibbs[idx-1], lb_gibbs[idx-1], M_gibbs, K_gibbs)
    
    # # get l_k and l_b
    lk_gibbs[idx] = get_lk(y, mu_gibbs, g_gibbs, sigma_2_gibbs[idx], lk_gibbs[idx-1]) # should just be passing it mu
    lb_gibbs[idx] = get_lb(y, lb_gibbs[idx-1], xi_gibbs[idx, ])
  
    # print statement for time  
    if (idx %% 50 == 0) print(paste("iteration:", idx, "in", round(Sys.time() - start_inner, 2)))
  }
  print(round(Sys.time() - start),2)
  
  
  # get burnin and remove 
  burn_in = floor(B*burn_in)
  
  beta_gibbs = beta_gibbs[-c(1:burn_in), ]
  mu_gibbs = mu_gibbs[-c(1:burn_in), ]
  sigma_2_gibbs = sigma_2_gibbs[-c(1:burn_in)]
  lb_gibbs = lb_gibbs[-c(1:burn_in)]
  lk_gibbs = lk_gibbs[-c(1:burn_in)]
  xi_gibbs = xi_gibbs[-c(1:burn_in), ]
  
  
  # write to csv file if running on server
  if (write)
  {
    write.csv(beta_gibbs,"beta_500.csv", row.names = TRUE)
    write.csv(mu_gibbs,"mu_500.csv", row.names = TRUE)
    write.csv(sigma_2_gibbs,"sigma_2_500.csv", row.names = TRUE)
    write.csv(lb_gibbs,"lb_500.csv", row.names = TRUE)
    write.csv(lk_gibbs,"lk_500.csv", row.names = TRUE)
  }
  
  return(list(beta = beta_gibbs, mu = mu_gibbs, sigma_2 = sigma_2_gibbs, xi = xi_gibbs, lb = lb_gibbs,
              lk = lk_gibbs))
}





