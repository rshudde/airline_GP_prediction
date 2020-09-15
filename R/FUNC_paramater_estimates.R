# This is a file to start doing paramater estimates for sampling mu
# rm(list = ls())
library(invgamma)
library(MASS)
source('~/Desktop/Summer2020/AirplanePaper/airline_GP_prediction/R/FUNC_woodchan_samples.R')

clean_y = function(y_data)
{
  to_return = y_data[!is.na(y_data)]
  return(to_return)
}

normalize_data = function(data)
{
  # get maximum in rows
  rowmax = apply(data, 1, function(x) max(x))
  
  # calculate the x_i tilda
  data = data / rowmax
  data = data[(complete.cases(data)), ] # remove NA values
  return(data)
}

# function to return a matern kernel
get_matern_values = function(l_k, r_mj)
{
    term_one = 1 + sqrt(5)*r_mj / l_k + 5*(r_mj)^2 / (3*l_k^2)
    exponential = exp(-sqrt(5)*r_mj / l_k)
    return(term_one * exponential)
}

get_matern = function(l_k, time_points)
{
  # call individual values for upper triangular matrix
  # y = y[!is.na(y)]
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

get_h_j = function(data, beta, knots, N)
{

  # TODO - pass the row as the data 
  # make sure to scale the data
  # TODO get rid of transpose
  w_it = (data %*% beta + 1)/2

  h_return = vector()

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

get_g_i = function(xi, h)
{
  # print(xi)
  # print(h)
  to_return = sum(t(xi) %*% h)
  return(to_return)
}

get_g = function(data, beta, knots, N, xi)
{
  w_it = (data %*% beta + 1) / 2
  g = vector()
  for (i in 1:nrow(data))
  {
    h_temp = get_h_j(data[i,], beta, knots, N)
    g_i = get_g_i(xi, h_temp)
    g[i] = g_i
  }
  return(g)
}


get_sigma_mu_post = function(sigma_2, sigma_mu, V_i)
{
  T_i = nrow(V_i)
  ones_vector = rep(1, T_i)
  
  # calculate inner parenthensis
  inner = t(ones_vector) %*% solve(V_i) %*% ones_vector + sigma_mu^(-2)

  # invert to return
  inner_inverted = inner^-1
  
  return(inner_inverted)
}

get_alpha_mu_post = function(alpha_mu, sigma_mu, sigma_mu_post, g_i, V_i, y)
{
  # set up terms needed
  y = y[!is.na(y)]
  T_i = nrow(V_i)
  
  term_one = alpha_mu^2 * sigma_mu^(-2)
  term_two = t(y - g_i) %*% solve(V_i) %*% rep(1, T_i)
  
  to_return = sigma_mu_post*(term_one + term_two)
  return(to_return)
}

get_mu_i = function(alpha_mu_post, sigma_mu_post)
{
  # value to return
  to_return = rnorm(1, alpha_mu_post, sigma_mu_post)
  
  return(to_return)
}

# TODO this is where we are no longer sure if we are sane 
# this is the function to calculate Yi = mu_i - gI-
vector_differences = function(y, mu_i, g_i)
{
  Ti = length(y)
  ones_vector = rep(1, Ti)
  
  inner = y - mu_i*ones_vector - g_i
  return(inner)
}

# function to sample from an inverse gamm for sigma squared
get_sigma_squared = function(a, b, y, M, mu, g)
{
  Ti = apply(y, 1, function(x) length(which(!is.na(x)))) 
  # print(Ti)
  
  c = a + sum(Ti)
  d = b
  
  for (i in 1:nrow(y))
  {
    y_noNa = y[i,][!is.na(y[i,])]
    term_one = vector_differences(y_noNa, mu[i], g[[i]])
    term_one = term_one[!(is.na(term_one))]
    term_two = M[[i]] + diag(rep(1, Ti[i]))

    temp_update = t(term_one) %*% solve(term_two) %*% term_one
    
    # if (temp_update <=0) print(temp_update)
    # print(t(term_one) %*% term_one)
    d = d + temp_update / 2
  }
  # print(paste("c", c, "\td", d))
  
  sigma_2 = rinvgamma(1, c, abs(d)) # TODO figure out if htis is correct
  
  return(sigma_2)
}

# function to calculate acceptance ratio for l_k
lk_acceptance = function(y, mu, g, sigma_2, l_k_prime, l_k)
{
  # indicator function part
  if (l_k_prime < 0.1 || l_k_prime > 1 || l_k < 0.1 ||  l_k > 1)
  {
    to_return = 0
  } else { # calcualtions assuming indicator = 1
    # calcualte first term outside of the product
    
    y_noNA = y[1,][!is.na(y[1,])]
    term_one = vector_differences(y_noNA, mu[1], g[1])
    term_one = term_one[!(is.na(term_one))]
    
    M_temp = get_matern(l_k, y_noNA)
    M_prime = get_matern(l_k_prime, y_noNA)
    
    V_temp = get_V_i(sigma_2, M_temp, get_K_i(sigma_2, M_temp))
    V_prime = get_V_i(sigma_2, M_prime, get_K_i(sigma_2, M_prime))
    
    term_two = solve(V_prime) - solve(V_temp)
    
    ratio = exp(-0.5 * t(term_one) %*% term_two %*% term_one)
    
    for (i in 2:nrow(y))
    {
      # calcualte proceeding terms in product 
      y_noNA = y[i,][!is.na(y[i,])]
      
      term_one = vector_differences(y_noNA ,  mu[i], g[i])
      term_one = term_one[!(is.na(term_one))]
      
      M_temp = get_matern(l_k, y_noNA )
      M_prime = get_matern(l_k_prime, y_noNA )
      
      V_temp = get_V_i(sigma_2, M_temp, get_K_i(sigma_2, M_temp))
      V_prime = get_V_i(sigma_2, M_prime, get_K_i(sigma_2, M_prime))
      
      term_two = solve(V_prime) - solve(V_temp)
      
      ratio = ratio * exp(-0.5 * t(term_one) %*% term_two %*% term_one)
    }
  
  to_return = min(1, (l_k_prime / l_k) * ratio)
  }
  
  return(to_return)
}

get_lk = function(y, mu, g, sigma_2, lk_0)
{
  epsilon = 0.001
  lk_t = lk_0 + 50*epsilon
  mod_diff = abs(lk_0 - lk_t)
    
  while (mod_diff > epsilon)
  {
    # step two - draw lk_prime 
    lk_prime = rexp(1, lk_t)
    u_t = runif(1, 0, 1)
    
    acceptance = lk_acceptance(y, mu, g, sigma_2, lk_prime, lk_t)
    
    lk_t1 = ifelse(u_t <= acceptance, lk_prime, lk_t)
    
    mod_diff = abs(lk_t1 - lk_t)
    lk_t = lk_t1
  }
  
  return(lk_t1)
}

get_H_matrix = function(data, beta, knots, N)
{
  H = vector()
  for (i in 1:nrow(data))
  {
    temp = get_h_j(data[i,], beta, knots, N)
    H = rbind(H, temp)
  }
  
  return(H)
}

# for sampling xis
psi_xi = function(y, mu, data, xi, beta, knots, N, sigma_2, l_k, M, K)
{
  sum_term = 0
  for (i in 1:nrow(y))
  {
    y_noNA = y[i,][!is.na(y[i,])]
    H_term = get_H_matrix(data[[i]], beta, knots, N)
    term_one = y_noNA - mu[i] * rep(1,length(y_noNA)) - H_term
    
    # constrcuting second term
    M_i = get_matern(l_k, y_noNA )
    K_i = get_K_i(sigma_2, M[[i]])
    term_two = get_V_i(sigma_2, M[[i]], K[[i]])

    sum_term = sum(t(term_one) %*% solve(term_two) %*% term_one) + sum_term
  }
  
  to_return = sum_term / 2
  return(to_return)
}

get_xi = function(xi_0, y, mu, data, beta, knots, N, sigma_2, l_k, l_b, M, K)
{
  # count = 1
  # step one
  theta = runif(1, 0, 2*pi)
  gamma = samp.WC(knots, l_b)
  
  xi_list = list()
  xi_proposed = xi_0 * cos(theta) + gamma * sin(theta)
  # xi_list[[count]] = xi_proposed
  
  # step two
  theta_min = theta - 2*pi
  theta_max = theta
  
  # step three
  zeta = runif(1, 0, 1)
  
  psi_old = psi_xi(y, mu, data, xi_0, beta, knots, N, sigma_2, l_k, M, K)
  psi_new = psi_xi(y, mu, data, xi_proposed, beta, knots, N, sigma_2, l_k, M, K)
  acceptance = min(1, exp(psi_old - psi_new))
  
  # continuation of step 3 - don't return unti lwe get something we accept 
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
      acceptance = min(1, exp(psi_old - psi_new))
    }
  }
  
  return(xi_proposed)
}

# TODO 
# - ask about alpha term
# - ask about which norm
psi_alpha = function(y, mu, data, xi, alpha, knots, N, sigma_2, l_k, M, K)
{
  if (alpha[1] > 0)
  {
    beta = alpha / sum(alpha^2)
    to_return = psi_xi(y, mu, data, xi, beta, knots, N, sigma_2, l_k, M, K)
  } else {
    to_return = 0
  }
  return(to_return)
}

get_beta = function(alpha_0, y, mu, data, xi, knots, N, sigma_2, l_k, M, K, c_2 = 10^5)
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
  # shortcut 
  zeta = runif(1, 0, 1)
  
  if (alpha_proposed[1] <= 0)
  {
    acceptance = 0
  } else
  {
    psi_old = psi_alpha(y, mu, data, xi, alpha_0, knots, N, sigma_2, l_k, M, K)
    psi_new = psi_alpha(y, mu, data, xi, alpha_proposed, knots, N, sigma_2, l_k, M, K)
    acceptance = min(1, exp(psi_old - psi_new))
  }
  
  # continuation of step 3 - don't return unti lwe get something we accept 
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

g_lb_value = function(lb_value, beta, knots, N, xi)
{
  
  # TODO - pass the row as the data 
  # make sure to scale the data
  # TODO get rid of transpose

  h_return = vector()
  
  for (i in 1:length(knots))
  {
    numerator = lb_value - knots[i]
    denominator = 1/N
    value = numerator / denominator
    
    inner = ifelse(abs(value) <= 1, value, 0) # indicator function part
    h_return[i] = inner
  }

    to_return = sum(t(xi) %*% h_return)
  
  return(h_return)
}

lb_acceptance = function(lb_prime, lb, y)
{
  indicator = TRUE
  product_value = 0
  
  for (i in 1:nrow(y) && indicator)
  {
    if (indicator)
    {
      if (((l_b < 1) || (l_b > 0.01)) || ((l_k < 1) || (l_k > 0.01)))
      {
        indicator = FALSE
      }
      
      # if indicator is still good, calculate term in product
      g_lb = g_lb_value(lb, beta, knots, N, xi)
      g_lb_prime = g_lb_value(lb_prime, beta, knots, N, xi)
      
      term_one = (g_lb - g_lb_prime)
      
      # stuff we need to calcualte v_i
      M_temp = get_matern(l_k, y_noNA )
      M_prime = get_matern(l_k_prime, y_noNA )
      
      V_temp = get_V_i(sigma_2, M_temp, get_K_i(sigma_2, M_temp))
      V_prime = get_V_i(sigma_2, M_prime, get_K_i(sigma_2, M_prime))
      
      
    }
  }
  
  # return 0 if the indicator condition is not met
  to_return = indicator * product
  
  return(product_value)
}




# 
# beta = c(0.3, 0.2, 0.4, 0.1)
# N = 10
# xi_0 = rnorm(N+1)
# 
# set.seed(1)
# l_k = 0.9
# sigma_2 = 2
# sigma_mu = 5
# alpha_mu = 3
# x1 = matrix(rnorm(12, 5, 5), 3, 4)
# x1 = rbind(x1, rep(NA, 4))
# y1 = c(5,6,7, NA)
# 
# # remove NA
# y1_noNA = y1[!is.na(y1)]
# 
# x2 = matrix(rnorm(16, 6, 6), 4, 4)
# y2 = c(8,9,20,11)
# beta = c(0.3, 0.2, 0.4, 0.1)
# 
# x1 = normalize_data(x1)
# x2 = normalize_data(x2)
# knots = seq(0, 1, 0.1)
# 
# M1_i = get_matern(l_k, y1_noNA)
# K1_i = get_K_i(sigma_2, M1_i)
# V1_i = get_V_i(sigma_2, M1_i, K1_i)
# 
# M2_i = get_matern(l_k, y2)
# K2_i = get_K_i(sigma_2, M2_i)
# V2_i = get_V_i(sigma_2, M2_i, K2_i)
# 
# g1 = get_g(x1, beta, knots, N, xi_0)
# g2 = get_g(x2, beta, knots, N, xi_0)
# 
# 
# # now for the main calculations
# sigma_mu_post1 = get_sigma_mu_post(sigma_2, sigma_mu, V1_i)
# alpha_mu_post1 = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post1, g1, V1_i, y1_noNA)
# 
# sigma_mu_post2 = get_sigma_mu_post(sigma_2, sigma_mu, V2_i)
# alpha_mu_post2 = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post2, g2, V2_i, y2)
# 
# mu_1i = get_mu_i(alpha_mu_post1, sigma_mu_post1)
# mu_2i = get_mu_i(alpha_mu_post2, sigma_mu_post2)
# 
# 
# y = matrix(c(y1, y2), ncol = 4, byrow = T)
# M = list(M1_i, M2_i)
# mu = c(mu_1i, mu_2i)
# g = c(g1, g2)
# data = list(x1, x2)
# sigma_2 = get_sigma_squared(0.1, 0.1, y, M, mu, g)
# 
# lk = get_lk(y, mu, g, sigma_2, 0.5)
# 
# 
# ###
# psi = psi_xi(y, mu, data, xi, beta, knots, N, sigma_2, lk)
# xi = get_xi(xi_0, y, mu, data, beta, knots, N, sigma_2, l_k, l_b = 1)
# 
# # psi = psi_alpha(y, mu, data, xi, alpha = 1, knots, N, sigma_2, l_k)
# get_beta(alpha_0 = rnorm(4, 0, 1), y, mu, data, xi, knots, N, sigma_2, l_k, c_2 = 10^5)
# 
