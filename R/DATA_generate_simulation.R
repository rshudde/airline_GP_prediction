# generate data to test on
generate_simulation_data = function(n_datasets, n_covariates, knots = seq(0, 1, length.out = 100), l_k = 2, sigma_2 = 0.5)
{
  # get the number of covariates and the amount of data 
  ncol = n_covariates
  data = list()
  
  # sample betas from -10 to 10
  beta = sample(-10:10, n_covariates, replace = FALSE)
  beta[1] = abs(beta[1]) # force the first beta to be positive 
  beta = beta/sum(beta) # ||beta|| = 1
  
  # generates xi values
  N = length(knots)
  # set.seed(1)
  # xi = rnorm(length(knots), 1, 4) # genreeate xi for testing purposes
  xi = rnorm(length(knots), 0, 1) # genreeate xi for testing purposes
  
  # to get x variables first
  for (i in 1:n_datasets)
  {
    nrow = sample(2:10, 1) # sample how many observations we will have 
    indices = sample(1:10, nrow) # get the indices (basically randomly picking which 'days of the week' we observe)
    # mean = sample(1:10, 1)
    # sd = sample(2:9, 1)
    
    mean = 0
    sd = 1
    
    # sample data from normal 
    temp_x = matrix(rnorm(nrow * ncol, mean, sd), nrow, ncol)
    # new way of normalizing code to get x_tilde 
    temp_x = temp_x / max(sqrt(rowSums(temp_x * temp_x))) 
    rownames(temp_x) = sort(indices)
    
    # add data to list 
    data[[i]] = temp_x
  }

  #  get random starting mu values
  mu = rnorm(length(data), 0, 1)
  
  # data generation via equation 5 from writeup 
  y_matrix = vector()
  for (i in 1:n_datasets)
  {
    # mu terms
    first_term = mu[[i]] %*% rep(1, nrow(data[[i]]))
    
    # g values from FUNC_paramater_estimates
    second_term = get_g(data[[i]], beta, knots, N, xi)
    
    # third term - eta values 
    M_i = get_matern(l_k, rownames(data[[i]]))
    K = get_K_i(sigma_2, M_i)
    third_term = mvtnorm::rmvnorm(1, rep(0, nrow(data[[i]])), K + sigma_2 * diag(1, nrow = nrow(data[[i]])))
    
    temp = rep(NA, 10) # rep the number of covariates
    y_temp = as.vector(first_term + second_term + third_term)
    
    # arrange so our rownames are corresponding to the days of the week we observe
    temp[as.integer(rownames(data[[i]]))] = y_temp
    
    y_matrix = rbind(y_matrix, temp)
  }
  
  # return our x, y, and beta values
  return(list(X = data, y = y_matrix, beta = beta, sigma_2 = sigma_2, l_k = l_k, sigma_2 = sigma_2,
              xi = xi, mu = mu))
}



