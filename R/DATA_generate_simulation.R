# get data
# nrow = 10

generate_simulation_data = function(n_datasets, n_covariates)
{
  ncol = n_covariates
  data = list()
  
  # to get x variables first
  for (i in 1:n_datasets)
  {
    nrow = sample(2:10, 1)
    indices = sample(1:10, nrow)
    mean = sample(1:10, 1)
    sd = sample(2:9, 1)
    
    temp_x = normalize_data(matrix(rnorm(nrow*ncol, mean, sd), nrow, ncol))
    rownames(temp_x) = sort(indices)
    
    data[[i]] = temp_x
  }
  
  # xi - values
  xi = rnorm(11, 0, 1)
  mu = rnorm(length(data), 0, 1)
  beta = sample(1:20, n_covariates)
  beta = beta/sum(beta)
  knots = seq(0, 1, 0.1)
  N = length(knots)
  
  l_k = 2
  sigma_2 = 4
  
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



