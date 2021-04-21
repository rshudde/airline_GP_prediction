# generate data to test on
generate_simulation_data = function(n_datasets, n_time, n_covariates,
                                    mu_true, beta_true, sigma_2_true, lK_true, 
                                    xi_true, sigmaB_2_true, lB_true, seed = 1)
{
  set.seed(seed)
  
  # default truth
  # if(missing(mu_true)) mu_true = c(0, runif(n_datasets-1, -5, 5))
  if (missing(mu_true))
  {
    mu_true = c(0, runif(n_datasets - 1, -5, 5))
    # mu_true = c(5, runif(n_datasets - 1, -5, 5))
  }
    
  if (missing(beta_true))
  {
    # alpha = sample(c(-10:-1, 1:10), n_covariates, replace = T)
    alpha = sample(-10:10, n_covariates, replace = T)
    # alpha = sort(abs(alpha)) # commmented out forcing all alphas to be positive
    
    # set the max to be alpha[1] and positive
    max_value = which.max(abs(alpha))
    alpha = c(abs(alpha[max_value]), alpha[-max_value])
    
    beta_true = alpha / sqrt(sum(alpha^2)) # ||beta|| = 1
    # beta_true = alpha / sum(abs(alpha)) # ||beta|| = 1
    
  }
  if (missing(sigma_2_true)) sigma_2_true = .25
  if (missing(lK_true)) lK_true = .3
  
  if (missing(sigmaB_2_true)) sigmaB_2_true = 5
  if (missing(lB_true)) lB_true = .3
  if (missing(xi_true))
  {
    n_Knots = 20
    knots = seq(0, 1, length.out = n_Knots)
    B_true = sigmaB_2_true*get_matern(lB_true, knots)
    xi_true = as.numeric(mvtnorm::rmvnorm(1, numeric(n_Knots), B_true))
    
  }else{
    
    n_Knots = length(xi_true)
    knots = seq(0, 1, length.out = n_Knots)
  }
  
  # to get x variables first
  X = time_idx = vector("list", n_datasets)
  c.X = rep(NA, n_datasets)
  for (i in 1:n_datasets)
  {
    # n_time_obs_i = sample(2:n_time, 1) # sample how many observations we will have
    # time_idx[[i]] = sort(sample(1:n_time, n_time_obs_i)) # get the indices (basically randomly picking which 'days of the week' we observe)
    n_time_obs_i = n_time
    time_idx[[i]] = 1:n_time
    
    # iid design
    X[[i]] = matrix(rnorm(n_time_obs_i * n_covariates, 0, 1), 
                    n_time_obs_i, n_covariates)
    
    X[[i]][,1] = rep(5, nrow(X[[i]]))
    
    c.X[i] = max(apply(X = X[[i]], MARGIN = 1,
                       FUN = function(r){sqrt(sum(r^2))}))
  }
  
  c.X = max(c.X)
  
  # data generation via equation 5 from writeup 
  y_matrix = matrix(nrow = n_datasets, ncol = n_time)
  Xtilde = w_true = g_true = vector("list", n_datasets)
  loglhood_true = 0
  for (i in 1:n_datasets)
  {
    ## g values from FUNC_paramater_estimates
    Xtilde[[i]] = X[[i]]/c.X
    w_true[[i]] = (as.numeric(Xtilde[[i]] %*% beta_true) + 1)/2
    Hmat_i = get_H_matrix(w_true[[i]], knots, n_Knots)
    g_true[[i]] = get_g(Hmat_i, xi_true)
    
    # third term - eta values 
    M_i = get_matern(lK_true, time_idx[[i]])
    K_i = get_K_i(sigma_2_true, M_i)
    V_i = get_V_i(sigma_2_true, K_i)
    third_term = 
      as.numeric(mvtnorm::rmvnorm(1, numeric(length(time_idx[[i]])),
                                  V_i))
    
    # y values
    y_matrix[i, time_idx[[i]]] = 
      mu_true[i] + g_true[[i]] + third_term
    
    # log likelihood
    loglhood_true = loglhood_true -
      as.numeric(determinant(V_i, log = T)$modulus)/2 -
      emulator::quad.form.inv(V_i,
                              y_matrix[i, time_idx[[i]]] - mu_true[i] - 
                                g_true[[i]])/2
  }
  
  # return our x, y, and beta values
  return(list(X = Xtilde, y = y_matrix, time_idx = time_idx,
              mu_true = mu_true, beta_true = beta_true,
              w_true = w_true, g_true = g_true, xi_true = xi_true,
              lB_true = lB_true, lK_true = lK_true, sigma_2_true = sigma_2_true,
              sigmaB_2_true = sigmaB_2_true, loglhood_true = loglhood_true))
}



