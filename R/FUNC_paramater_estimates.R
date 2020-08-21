# This is a file to start doing paramater estimates for sampling mu

# function to return a matern kernel
get_matern = function(l_k, d)
{
    
}


# function to get our K matrix
get_K = function(sigma_2, l_k)
{
  M = get_matern(d, l_k)
  K = sigma_2*M
  
  return(K)
}


# function to get V
get_V = function(sigma_2, l_k, d)
{
  K = get_K(sigma_2, l_k)
  iden = rep(1, nrow(M))
  iden = diag(iden)
  
  V = K + sigma_2 * iden
  return(V)
}

get_sigma_mu_post = function(simga_2, l_k, sigma_mu, d)
{
  
}

get_h = function(data, beta)
{
  # make sure to scale the data
}

get_g = function(xi, data, beta)
{
  
}

get_alpha_my_post = function(alpha_my, sigma_my, y_response, xi, data, beta, sigma_2, l_k, d)
{
  
}

get_mu = function(alpha_my, sigma_my, y_response, xi, data, beta, sigma_2, l_k, d)
{
  
}