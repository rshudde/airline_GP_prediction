# This is a file to start doing paramater estimates for sampling mu
rm(list = ls())
library(invgamma)
library(MASS)

# function to return a matern kernel
get_matern_values = function(l_k, r_mj)
{
    term_one = 1 + sqrt(5)*r_mj / l_k + 5*r_mj^2 / l_k^2
    exponential = exp(-sqrt(5*r_mj / l_k))
    
    return(term_one * exponential)
}

get_matern = function(l_k, y)
{
  # call individual values for upper triangular matrix
  distance_numeric = dist(y)
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


get_h = function(data, beta)
{
  # make sure to scale the data
  w_it = (t(data) %*% beta + 1)/2
  inner = 1 - abs(w_it)
  inner = ifelse(abs(w_it) <= 1, inner, 0) # indicator function part
  
  return(inner)
}

get_g_i = function(xi, h)
{
  to_return = sum(xi * h)
  return(to_return)
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
  
  c = a + sum(Ti)
  d = b
  
  for (i in 1:nrow(y))
  {
    term_one = vector_differences(y[i,], mu[i], g[i])
    term_one = term_one[!(is.na(term_one))]
    term_two = M[[i]] + diag(rep(1, Ti[i]))
    
    temp_update = t(term_one) %*% solve(term_two) %*% term_one
    d = d + temp_update / 2
  }
  
  sigma_2 = rinvgamma(1, c, d) # TODO figure out if htis is correct
  
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


# function to calculate acceptance ratio for l_k
lb_acceptance = function(y, mu, g, sigma_2, l_k_prime, l_k)
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
    
    to_return = min(1, ratio)
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

set.seed(1)
l_k = 0.9
sigma_2 = 2
sigma_mu = 5
alpha_mu = 3
x1 = matrix(rnorm(12, 5, 5), 3, 4)
x1 = x1 / max(x1)
x1 = rbind(x1, rep(NA, 4))
y1 = c(5,6,7, NA)
beta1 = c(0.3, 0.2, 0.5, NA)
x1i = rnorm(4)

# remove NA
x1_noNA = x1[(complete.cases(x1)), ]
y1_noNA = y1[!is.na(y1)]
beta1_noNA = beta1[!is.na(beta1)]


x2 = matrix(rnorm(16, 6, 6), 4, 4)
x2 = x2 / max(x2)
y2 = c(8,9,20,11)
beta2 = c(0.3, 0.2, 0.4, 0.1)
x2i = rnorm(4)

M1_i = get_matern(l_k, y1_noNA)
K1_i = get_K_i(sigma_2, M1_i)
V1_i = get_V_i(sigma_2, M1_i, K1_i)

M2_i = get_matern(l_k, y2)
K2_i = get_K_i(sigma_2, M2_i)
V2_i = get_V_i(sigma_2, M2_i, K2_i)

h1 = get_h(x1_noNA, beta1_noNA)
h2 = get_h(x2, beta2)
g1_i = get_g_i(x1i, h1)
g2_i = get_g_i(x2i, h2)

# now for the main calculations
sigma_mu_post1 = get_sigma_mu_post(sigma_2, sigma_mu, V1_i)
alpha_mu_post1 = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post1, g1_i, V1_i, y1_noNA)

sigma_mu_post2 = get_sigma_mu_post(sigma_2, sigma_mu, V2_i)
alpha_mu_post2 = get_alpha_mu_post(alpha_mu, sigma_mu, sigma_mu_post2, g2_i, V2_i, y2)

mu_1i = get_mu_i(alpha_mu_post1, sigma_mu_post1)
mu_2i = get_mu_i(alpha_mu_post2, sigma_mu_post2)


y = matrix(c(y1, y2), ncol = 4, byrow = T)
M = list(M1_i, M2_i)
mu = c(mu_1i, mu_2i)
g = c(g1_i, g2_i)
sigma_2 = get_sigma_squared(0.1, 0.1, y, M, mu, g)

lk = get_lk(y, mu, g, sigma_2, 0.5)

