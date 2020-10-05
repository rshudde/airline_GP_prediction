# This is a file to start doing paramater estimates for sampling mu
# rm(list = ls())
library(invgamma)
library(MASS)
library(FastGP)
# source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_woodchan_samples.R')


clean_y = function(y_data)
{
  to_return = y_data[!is.na(y_data)]
  return(to_return)
}

# function to get the individual values for a Matern(5/2) kernel
get_matern_values = function(l_k, r_mj)
{
    term_one = 1 + sqrt(5) * r_mj / l_k + 5 * (r_mj)^2 / (3 * l_k^2)
    exponential = exp(-sqrt(5) * r_mj / l_k)
    return(term_one * exponential)
}

# function to actually create matern matrix
get_matern = function(l_k, time_points)
{
  # get the distance matrix 
  distance_numeric = dist(time_points)
  distance_matrix = as.matrix(distance_numeric)
  
 # get full matrix
  M_i = get_matern_values(l_k, distance_matrix)
  
  # return M matfix
  return(M_i)
}


# function to get our K matrix
get_K_i = function(sigma_2, M_i)
{
  K_i = sigma_2*M_i
  
  return(K_i)
}

# function to get V
get_V_i = function(sigma_2, M_i, K_i)
{
  iden = rep(1, nrow(M_i))
  iden = diag(iden)
  
  V_i = K_i + sigma_2 * iden
  return(V_i)
}

# calculates h vector from equation 6
get_h_j = function(data, beta, knots, N)
{

  # TODO - pass the row as the data 
  # TODO get rid of transpose
  
  # get w_it values 
  w_it = (data %*% beta + 1)/2

  h_return = vector()

  # follow equation 6 
  for (i in 1:length(knots))
  {
    numerator = w_it - knots[i]
    denominator = 1/N
    value = numerator / denominator
    
    inner = ifelse(abs(value) <= 1, value, 0) # indicator function part
    h_return[i] = inner
  }

  return(h_return)
}

# calculatesindividual g values frm euation 6
get_g_i = function(xi, h)
{
  to_return = sum(crossprod(xi, h))
  
  return(to_return)
}

# g function from equation 6
get_g = function(data, beta, knots, N, xi)
{
  # w_it = (data %*% beta + 1) / 2
  w_it = (crossprod(t(data), beta)+ 1) / 2
  
  g = vector()
  for (i in 1:nrow(data))
  {
    h_temp = get_h_j(data[i,], beta, knots, N)
    g_i = get_g_i(xi, h_temp)
    g[i] = g_i
  }
  return(g)
}

# function to get sigma_mu for the calculations of mu_i
get_sigma_mu_post = function(sigma_2, sigma_mu, V_i)
{
  # calculate T_i
  T_i = nrow(V_i)
  ones_vector = rep(1, T_i)
  
  # calculate sigma_mu_post
  inner_part_one = crossprod(ones_vector, chol2inv(V_i))
  inner = inner_part_one %*% ones_vector + sigma_mu^(-2)
  
  # invert to return
  inner_inverted = inner^-1
  
  return(inner_inverted)
}

# function to get alpha_mu for the calculations of mu_i
get_alpha_mu_post = function(alpha_mu, sigma_mu, sigma_mu_post, g_i, V_i, y)
{
  # set up terms needed
  y = y[!is.na(y)] # remove any NA values
  T_i = nrow(V_i)
  
  # term_one = alpha_mu^2 * sigma_mu^(-2) # TODO this may be wrong, fixed possibly below
  term_one = alpha_mu * sigma_mu^(-1) # fixed some things here
  
  term_two_part_one = crossprod(y - g_i, chol2inv(V_i))
  term_two = term_two_part_one %*% rep(1, T_i)
  
  to_return = sigma_mu_post * (term_one + term_two)
  return(to_return)
}

# function to draw the mu values
get_mu_i = function(alpha_mu_post, sigma_mu_post)
{
  # value to return
  to_return = rnorm(1, alpha_mu_post, sigma_mu_post)
  
  return(to_return)
}

# TODO this is where we are no longer sure if we are sane 
# this is the function to calculate Yi = mu_i - gI
vector_differences = function(y, mu_i, g_i)
{
  Ti = length(y)
  ones_vector = rep(1, Ti)
  
  inner = y - mu_i * ones_vector - g_i
  return(inner)
}

# function to sample from an inverse gamma for sigma squared
get_sigma_squared = function(a, b, y, M, mu, g)
{
  # count total numebr of non-NA y values for sum
  Ti = apply(y, 1, function(x) length(which(!is.na(x)))) 

  c = a + sum(Ti) # first value for IG draw 
  d = b
  
  for (i in 1:nrow(y))
  {
    # remove the y values 
    y_noNa = y[i,][!is.na(y[i,])]
    
    # the two terms in the multiplication 
    term_one = vector_differences(y_noNa, mu[i], g[[i]])
    term_one = term_one[!(is.na(term_one))]
    term_two = M[[i]] + diag(rep(1, Ti[i]))

    # do matrix multiplication 
    temp_update_part_one = crossprod(term_one, solve(term_two)) 
    temp_update = temp_update_part_one %*% term_one
    
    d = d + temp_update / 2
  }

  # do inverse gamma draw
  sigma_2 = rinvgamma(1, c, abs(d)) # TODO figure out if this is correct
  
  return(sigma_2)
}

# function to calculate acceptance ratio for l_k
lk_acceptance = function(y, mu, g, sigma_2, lk_prime, l_k)
{
  # indicator function part
  if (lk_prime < 0.1 || lk_prime > 1 || l_k < 0.1 ||  l_k > 1)
  {
    to_return = 0
  } else { # calcualtions assuming indicator = 1
    # calcualte first term outside of the product
    y_noNA = y[1,][!is.na(y[1,])]
    term_one = vector_differences(y_noNA, mu[1], g[[1]])
    term_one = term_one[!(is.na(term_one))]
    
    # get the two v terms necessary
    M_temp = get_matern(l_k, y_noNA)
    M_prime = get_matern(lk_prime, y_noNA)
    
    V_temp = get_V_i(sigma_2, M_temp, get_K_i(sigma_2, M_temp))
    V_prime = get_V_i(sigma_2, M_prime, get_K_i(sigma_2, M_prime))
    
    term_two = chol2inv(V_prime) - chol2inv(V_temp)
    
    # do matrix multiplication
    matrix_part = crossprod(term_one, term_two)
    matrix_part = matrix_part %*% term_one
    
    # calculate ratio
    ratio = exp(-0.5 * as.numeric(matrix_part))
    
    for (i in 2:nrow(y))
    {
      # calcualte proceeding terms in product 
      y_noNA = y[i,][!is.na(y[i,])]
      
      term_one = vector_differences(y_noNA ,  mu[i], g[[i]])
      term_one = term_one[!(is.na(term_one))]
      
      # get two v terms necessary
      M_temp = get_matern(l_k, y_noNA )
      M_prime = get_matern(lk_prime, y_noNA )
      
      V_temp = get_V_i(sigma_2, M_temp, get_K_i(sigma_2, M_temp))
      V_prime = get_V_i(sigma_2, M_prime, get_K_i(sigma_2, M_prime))
      
      term_two = chol2inv(V_prime) - chol2inv(V_temp)
      
      # do matrix multiplication
      matrix_part = crossprod(term_one, term_two)
      matrix_part = matrix_part %*% term_one
      
      # calculate ratio
      ratio = ratio * exp(-0.5 * as.numeric(matrix_part))
    }
  
    # calculate the minimum for the return value
  to_return = min(1, (lk_prime / l_k) * as.numeric(ratio))
  }
  
  return(as.numeric(to_return))
}

get_lk = function(y, mu, g, sigma_2, lk_0)
{
  # small epsilon 
  epsilon = 0.001
  mod_diff = 0.01 # mod_diff starts > epsilon to force into loop
  lk_t = lk_0 

  while (mod_diff > epsilon)
  {
    # step two - draw lk_prime 
    lk_prime = rexp(1, lk_t)
    
    # draw the uniform variable 
    u_t = runif(1, 0, 1)
    
    # calculate the acceptance
    acceptance = lk_acceptance(y, mu, g, sigma_2, lk_prime, lk_t)
    
    # update lk_t1 based on acceptance
    lk_t1 = ifelse(u_t <= acceptance, lk_prime, lk_t)
    
    # reacalculate the absolute value of the differences
    mod_diff = abs(lk_t1 - lk_t)
    
    # update lk_k (either it will have stayed the same or updated)
    lk_t = lk_t1
  }
  
  return(lk_t1)
}

# get the H matrix for calculating the xi
get_H_matrix = function(data, beta, knots, N)
{
  row = nrow(data)
  col = length(knots)
  H = matrix(rep(0, row * col), nrow = row, ncol = col)

  for (i in 1:nrow(data))
  {
    H[i, ] = get_h_j(data[i,], beta, knots, N)
  }
  
  return(H)
}

# function for negative log likelihood for xis
psi_xi = function(y, mu, data, xi, beta, knots, N, sigma_2, l_k, M, K)
{
  sum_term = 0
  for (i in 1:nrow(y))
  {
    # get y values with no NAs
    y_noNA = y[i,][!is.na(y[i,])]
    
    # construct H matrix
    H_term = get_H_matrix(data[[i]], beta, knots, N)
    term_one = y_noNA - mu[i] * rep(1,length(y_noNA)) - H_term
    
    # constrcuting second term
    M_i = get_matern(l_k, y_noNA )
    K_i = get_K_i(sigma_2, M[[i]])
    term_two = get_V_i(sigma_2, M[[i]], K[[i]])

    # matrix multiplication
    matrix_part = crossprod(term_one, chol2inv(term_two))
    matrix_part = matrix_part %*% term_one
    sum_term = sum(matrix_part) + sum_term
  }
  
  to_return = sum_term / 2
  return(to_return)
}

# function to sample the xi values
get_xi = function(xi_0, y, mu, data, beta, knots, N, sigma_2, l_k, l_b, M, K)
{
  # step one
  theta = runif(1, 0, 2*pi)
  gamma = samp.WC(knots, l_b)
  
  xi_list = list()
  xi_proposed = xi_0 * cos(theta) + gamma * sin(theta)

  # step two
  theta_min = theta - 2*pi
  theta_max = theta
  
  # step three
  zeta = runif(1, 0, 1)
  
  # calculate new and old psi values
  psi_old = psi_xi(y, mu, data, xi_0, beta, knots, N, sigma_2, l_k, M, K)
  psi_new = psi_xi(y, mu, data, xi_proposed, beta, knots, N, sigma_2, l_k, M, K)
  acceptance = min(1, exp(psi_old - psi_new))
  
  # continuation of step 3 - don't return until we get something we accept 
  if (acceptance <= zeta)
  {
    while (acceptance <= zeta )
    {
      # step a
      if (theta < 0)
      {
        theta_min = theta
      } else {
        theta_max = theta
      }
      
      # step b 
      theta = runif(1, theta_min, theta_max)
      
      # step c
      xi_proposed = xi_0 * cos(theta) + gamma * sin(theta)
      
      # step d
      psi_old = psi_function(y, mu, data, xi_0, beta, knots, N, sigma_2, l_k, M, K)
      psi_new = psi_function(y, mu, data, xi_proposed, beta, knots, N, sigma_2, l_k, M, K)
      
      # calculate new acceptance
      acceptance = min(1, exp(psi_old - psi_new))
    }
  }
  
  return(xi_proposed)
}

# TODO 
# - ask about alpha term
# - ask about which norm
# function for calculating likelihood of alpha
psi_alpha = function(y, mu, data, xi, alpha, knots, N, sigma_2, l_k, M, K)
{
  # check indicator values
  if (alpha[1] > 0)
  {
    beta = alpha / sum(alpha^2) # changed this to be the l2 norm
    to_return = psi_xi(y, mu, data, xi, beta, knots, N, sigma_2, l_k, M, K)
  } else {
    to_return = 0
  }
  return(to_return)
}

# function to get alpha values
get_alpha = function(alpha_0, y, mu, data, xi, knots, N, sigma_2, l_k, M, K, c_2 = 10^5)
{
  # step one
  theta = runif(1, 0, 2*pi)
  gamma = mvtnorm::rmvnorm(1, mean = rep(0, length(alpha_0)), sigma = c_2 * diag(rep(1, length(alpha_0))))
  gamma = as.vector(gamma)
  alpha_proposed = alpha_0 * cos(theta) + gamma * sin(theta)
  
  # step two
  theta_min = theta - 2*pi
  theta_max = theta
  
  # step 3
  # shortcut to check if we need to go into the loop
  zeta = runif(1, 0, 1)
  
  if (alpha_proposed[1] <= 0)
  {
    acceptance = 0
  } else
  {
    psi_old = psi_alpha(y, mu, data, xi, alpha_0, knots, N, sigma_2, l_k, M, K)
    psi_new = psi_alpha(y, mu, data, xi, alpha_proposed, knots, N, sigma_2, l_k, M, K)
    # calculate new acceptance values
    acceptance = min(1, exp(psi_old - psi_new))
  }
  
  # continuation of step 3 - don't return until we get something we accept 
  if (acceptance <= zeta)
  {
    while (acceptance <= zeta )
    {
      # step a
      if (theta < 0)
      {
        theta_min = theta
      } else {
        theta_max = theta
      }
      
      # step b 
      theta = runif(1, theta_min, theta_max)
      # step c
      alpha_proposed = alpha_0 * cos(theta) + gamma * sin(theta)
      
      # shortcut
      if (alpha_proposed[1] <= 0)
      {
        acceptance = 0
      } else
      {
        psi_old = psi_alpha(y, mu, data, xi, alpha_0, knots, N, sigma_2, l_k, M, K)
        psi_new = psi_alpha(y, mu, data, xi, alpha_proposed, knots, N, sigma_2, l_k, M, K)
        acceptance = min(1, exp(psi_old - psi_new))
      }
    }
  }
  
  return(alpha_proposed)
}

## LB functions below
lb_acceptance = function(y, lb, lb_prime, xi)
{
  # print(knots)
  # indicator function part
  if (lb_prime < 0.1 || lb_prime > 1 || lb < 0.1 ||  lb > 1)
  {
    to_return = 0
  } else { # calcualtions assuming indicator = 1
    # calcualte first term outside of the product
    term_one = xi
    
    # stuff we need to calcualte v_i
    M_lb = get_matern(lb, xi)
    M_lb_prime = get_matern(lb_prime, xi)

    # use tinv to invert M here
    term_two = tinv(M_lb) - tinv(M_lb_prime)
    
    # matrix multiplication
    matrix_part = crossprod(term_one, term_two)
    matrix_part = matrix_part %*% term_one
    ratio = exp(-0.5 * matrix_part)
    
    # minimum value for eturn 
    to_return = min(1, (lb_prime / lb) * as.numeric(ratio))
  }
  
  # return correct type - don't return matrix from matrix multiplication
  to_return = as.numeric(to_return)
  return(to_return)
}

# function to calculate lb updates
get_lb = function(y, lb_0, xi)
{
  epsilon = 0.001
  mod_diff = 0.01 # start mod_diff small to force into loop
  lb_t = lb_0
  
  while (mod_diff > epsilon)
  {
    # step two - draw lb_prime 
    lb_prime = rexp(1, lb_t)

    #  draw uniform valuable
    u_t = runif(1, 0, 1)
    
    # calculate new acceptance 
    acceptance = lb_acceptance(y, lb_t, lb_prime, xi)

    # update lb_t1 based on acceptance
    lb_t1 = ifelse(u_t <= acceptance, lb_prime, lb_t)
    
    # c alculate new mod_difference
    mod_diff = abs(lb_t1 - lb_t)
    
    # update the lb_t1 value (it will either chance or not based on lb_t1)
    lb_t = lb_t1
  }
  
  return(lb_t1)
}

