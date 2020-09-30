# generate data to test on
generate_simulation_data = function(n_datasets, n_covariates, knots = seq(0, 1, length.out = 100))
{
  ncol = n_covariates
  data = list()
  # xi - values
  beta = sample(1:20, n_covariates)
  beta = beta/sum(beta) # norm beta = 1
  N = length(knots)
  xi = rnorm(length(knots), 0, 1) # genreeate xi for testing purposes
  
  # to get x variables first
  for (i in 1:n_datasets)
  {
    nrow = sample(2:10, 1)
    indices = sample(1:10, nrow)
    mean = sample(1:10, 1)
    sd = sample(2:9, 1)
    
    temp_x = matrix(rnorm(nrow*ncol, mean, sd), nrow, ncol)
    temp_x = temp_x/max(sqrt(rowSums(temp_x*temp_x))) # new way of normalizing code - x_tilda 
    rownames(temp_x) = sort(indices)
    
    data[[i]] = temp_x
  }

  mu = rnorm(length(data), 0, 1)
  
  # assumed l_k and sigma_2 variables 
  l_k = 2
  sigma_2 = 0.5
  
  y_matrix = vector()
  for (i in 1:n_datasets)
  {
    first_term = mu[[i]] %*% rep(1, nrow(data[[i]]))
    
    second_term = get_g(data[[i]], beta, knots, N, xi)
    
    # third term
    M_i = get_matern(l_k, rownames(data[[i]]))
    K = get_K_i(sigma_2, M_i)
    third_term = mvtnorm::rmvnorm(1, rep(0, nrow(data[[i]])), K + sigma_2*diag(1, nrow = nrow(data[[i]])))
    
    temp = rep(NA, 10) # rep the number of covariates
    y_temp = as.vector(first_term + second_term + third_term)
    
    temp[as.integer(rownames(data[[i]]))] = y_temp
    
    y_matrix = rbind(y_matrix, temp)
  }
  
  return(list(X = data, y = y_matrix, beta = beta))
}



