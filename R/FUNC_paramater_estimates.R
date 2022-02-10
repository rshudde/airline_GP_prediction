# This is a file to start doing paramater estimates for sampling




# function to get the individual values for a Matern(5/2) kernel
get_matern_values = function(l_k, r_mj)
{
  term_one = 1 + (sqrt(5) * r_mj) / l_k + (5 * (r_mj^2)) / (3 * (l_k^2))
  exponential = exp(-(sqrt(5) * r_mj) / l_k)
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
get_V_i = function(sigma_2, K_i)
{
  if (!is.matrix(K_i)) 
  {
    stop("K_i needs to be a matrix for get_V_i()")
  }
  V_i = K_i + sigma_2 * diag(nrow(K_i))
  return(V_i)
}


# function to draw the mu values
get_mu = function(y, n_datasets, g, V_mat, time_idx, sigma_2_mu, alpha_mu)
{
  mu = rep(NA, n_datasets)
  
  for (i in 1:n_datasets)
  {
    sigma_2_mu_post_i = 1/(quad.form.inv(V_mat[[i]], rep(1, length(time_idx[[i]]))) + 1/sigma_2_mu)
    alpha_mu_post_i = sigma_2_mu_post_i*(alpha_mu/sigma_2_mu +
         quad.3form.inv(V_mat[[i]], y[i,time_idx[[i]]] - g[[i]],
                        rep(1, length(time_idx[[i]]))))

    mu[i] = rnorm(1, alpha_mu_post_i, sqrt(sigma_2_mu_post_i))
  }
  
  return(mu)
}


# get the H matrix
get_H_matrix = function(w, knots, n_Knots)
{
  t(mapply(i = 1:length(w),
           FUN = function(i){
             
             pmax(1 - abs((w[i] - knots)*(n_Knots - 1)), 0)
           }))
}


# g function from equation 6
get_g = function(H_mat, xi)
{
  as.numeric(Matrix::tcrossprod(H_mat, t(xi)))
}


# function for calculating likelihood of alpha
psi_alpha = function(alpha, y, n_datasets, time_idx, dat,
                     mu, xi, V_mat, knots, n_Knots)
{
  beta = alpha/sqrt(sum(alpha^2))
  # beta = alpha / sum(abs(alpha))
  w = H_mat = g = list()
  negloglhood = 0
  for (i in 1:n_datasets)
  {
    w[[i]] = (as.numeric(Matrix::tcrossprod(dat[[i]], t(beta))) + 1)/2

    H_mat[[i]] = get_H_matrix(w[[i]], knots, n_Knots)
    g[[i]] = as.numeric(Matrix::tcrossprod(H_mat[[i]], t(xi)))
    negloglhood = negloglhood +
      emulator::quad.form.inv(V_mat[[i]],
                              y[i,time_idx[[i]]] - mu[i] - g[[i]])
  }
  
  return(list("negloglhood" = negloglhood/2, "beta" = beta, 
              "w" = w, "H_mat" = H_mat, "g" = g))
}


# function for MH sampling to get alpha values
get_alpha = function(alpha_0, y, n_datasets, time_idx, dat, n_covariates,
                     mu, xi, V_mat, knots, n_Knots, c_prior)
{
  # step one
  theta = runif(1, 0, 2*pi)
  alpha_prior = rnorm(n_covariates, 0, c_prior)
  alpha_proposed = cos(theta) * alpha_0 + sin(theta) * alpha_prior
  
  # step two
  theta_min = theta - 2 * pi
  theta_max = theta
  
  # old negative loglikelihood
  psi_out_old = psi_alpha(alpha_0, y, n_datasets, time_idx, dat,
                          mu, xi, V_mat, knots, n_Knots)
  
  # new negative loglikelihood
  if (alpha_proposed[1] <= 0)
  {
    acceptance = 0
  } else
  {
    psi_out_new = psi_alpha(alpha_proposed, y, n_datasets, time_idx, dat,
                            mu, xi, V_mat, knots, n_Knots)
    
    # calculate new acceptance value
    acceptance = min(1, exp(psi_out_old$negloglhood - psi_out_new$negloglhood))
  }
  
  # step 3
  zeta = runif(1, 0, 1)
  
  # continuation of step 3 - don't return until we get something we accept 
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
    alpha_proposed = cos(theta) * alpha_0 + sin(theta) * alpha_prior
    
    # new negative loglikelihood
    # reject if alpha[1] is negative because we want beta[1] to be positive
    if (alpha_proposed[1] <= 0)
    {
      acceptance = 0
    } else
    {
      psi_out_new = psi_alpha(alpha_proposed, y, n_datasets, time_idx, dat,
                              mu, xi, V_mat, knots, n_Knots)
      
      # calculate new acceptance value
      acceptance = min(1, exp(psi_out_old$negloglhood - psi_out_new$negloglhood))
    }
  }
  
  return(c(psi_out_new, list("alpha" = alpha_proposed)))
}


# function to sample from an inverse gamma for sigma squared
get_sigma_2 = function(a, b, y, n_datasets, n_nonNA_y, time_idx,
                       mu, M_mat, g)
{
  # calculating shape and rate of inv gamma
  rate_term = 0
  for (i in 1:n_datasets)
  {
    rate_term = rate_term +
      emulator::quad.form.inv(M_mat[[i]] + diag(length(time_idx[[i]])),
                              y[i,time_idx[[i]]] - mu[i] - g[[i]])
  }
  
  # do inverse gamma draw
  sigma_2 = invgamma::rinvgamma(n = 1, shape = a + n_nonNA_y/2,
                                rate = b + rate_term/2)
  
  return(list("sigma_2" = sigma_2, "rate" = rate_term))
}


# function for calculating negative loglikelihood for xi
psi_xi = function(xi, y, n_datasets, time_idx,
                  mu, H_mat, V_mat)
{
  g = list()
  negloglhood = 0
  for (i in 1:n_datasets)
  {
    g[[i]] = as.numeric(Matrix::tcrossprod(H_mat[[i]], t(xi)))
    negloglhood = negloglhood +
      (emulator::quad.form.inv(V_mat[[i]],
                               y[i,time_idx[[i]]] - mu[i] - g[[i]]))/2
  }
  
  return(list("negloglhood" = negloglhood, "g" = g))
}


# function for MH sampling to get xi values
get_xi = function(xi_0, sigmaB_2, y, n_datasets, time_idx,
                  mu, H_mat, V_mat, lb, knots, nNeighbour, NNGP = FALSE)
{
  # step one
  theta = runif(1, 0, 2*pi)
  
  if (!NNGP)
  {
    # WC sampler
    xi_prior = samp.WC(knots, lb, 5/2, sigmaB_2) 
  } else { # non WC sampler
    nknots = length(knots)
    xi_prior = rep(NA, nknots)
    xi_prior[1] = rnorm(1, 0, sigmaB_2)
    for (j in 2:nknots)
    {
      b_j = sigmaB_2*get_matern_values(lb, abs(knots[max(1,j - nNeighbour):(j - 1)] - knots[j]))
      A_j = forwardsolve(sigmaB_2*get_matern_values(lb, as.matrix(dist(knots[max(1,j - nNeighbour):(j - 1)]))),b_j)
      d_j = sigmaB_2 - sum(A_j*b_j)
      xi_prior[j] = sum(xi_prior[max(1,j - nNeighbour):(j - 1)]*A_j) +
        rnorm(1, 0, sqrt(d_j))
    }
  }
  
  xi_proposed = cos(theta) * xi_0 + sin(theta) * xi_prior
  
  # step two
  theta_min = theta - 2 * pi
  theta_max = theta
  
  # old negative loglikelihood
  psi_out_old = psi_xi(xi_0, y, n_datasets, time_idx, mu,
                       H_mat, V_mat)
  
  # new negative loglikelihood
  psi_out_new = psi_xi(xi_proposed, y, n_datasets, time_idx, mu,
                       H_mat, V_mat)
  
  # calculate new acceptance value
  acceptance = min(1, exp(psi_out_old$negloglhood - psi_out_new$negloglhood))
  
  # step 3
  zeta = runif(1, 0, 1)
  
  # continuation of step 3 - don't return until we get something we accept 
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
    xi_proposed = cos(theta) * xi_0 + sin(theta) * xi_prior
    
    # new negative loglikelihood
    psi_out_new = psi_xi(xi_proposed, y, n_datasets, time_idx, mu,
                         H_mat, V_mat)
    
    # calculate new acceptance value
    acceptance = min(1, exp(psi_out_old$negloglhood - psi_out_new$negloglhood))
  }
  
  return(c(psi_out_new, list("xi" = xi_proposed)))
}


# function to sample from an inverse gamma for sigmaB squared
get_sigmaB_2 = function(a, b, xi, lb, knots, n_Knots)
{
  # calculating shape and rate of inv gamma
  rate_term = (emulator::quad.form.inv(get_matern(lb, knots), xi))/2
  
  # do inverse gamma draw
  sigmaB_2 = invgamma::rinvgamma(n = 1, 
                                 shape = a + n_Knots/2,
                                 rate = b + rate_term)
  
  return(sigmaB_2)
}

# function for MH sampling to get xi values
get_xi_c = function(xi_0, sigmaB_2, y, n_datasets, time_idx,
                  mu, H_mat, V_mat, lb, knots, nNeighbour, NNGP = FALSE)
{
  # step one
  if (!NNGP)
  {
    xi_prior = samp.WC(knots, lb, 5/2, sigmaB_2)
  } else {
    nknots = length(knots)
    xi_prior = rep(NA, nknots)
    xi_prior[1] = rnorm(1, 0, sigmaB_2)
    for (j in 2:nknots)
    {
      b_j = sigmaB_2*get_matern_values(lb, abs(knots[max(1,j - nNeighbour):(j - 1)] - knots[j]))
      A_j = forwardsolve(sigmaB_2*get_matern_values(lb, as.matrix(dist(knots[max(1,j - nNeighbour):(j - 1)]))),b_j)
      d_j = sigmaB_2 - sum(A_j*b_j)
      xi_prior[j] = sum(xi_prior[max(1,j - nNeighbour):(j - 1)]*A_j) +
        rnorm(1, 0, sqrt(d_j))
    }
  }

  # step two
  temp = get_xi_backend_c(xi_0, sigmaB_2, y, n_datasets, time_idx,
                          mu, H_mat, V_mat, lb, knots, xi_prior)

  return(temp)
  # return(c(psi_out_new, list("xi" = xi_proposed)))
}


# # function to remove NA values, not actually used
# clean_y = function(y_data)
# {
#   to_return = y_data[!is.na(y_data)]
#   return(to_return)
# }


# # calculates h vector from equation 6
# get_h_j = function(data, beta, knots)
# {
#   
#   N = length(knots)
#   
#   # get w_it values 
#   w_it = (data %*% beta + 1)/2
# 
#   h_return = vector()
# 
#   # follow equation 6 
#   for (i in 1:N)
#   {
#     numerator = w_it - knots[i]
#     denominator = 1/N
#     value = numerator / denominator
#     
#     inner = ifelse(abs(value) <= 1, 1 - abs(value), 0) # indicator function part
#     h_return[i] = inner
#   }
# 
#   return(h_return)
# }
# 
# 
# # calculates individual g values frm euation 6
# get_g_i = function(xi, h)
# {
#   to_return = sum(crossprod(xi, h))
#   
#   return(to_return)
# }


# # function to get sigma_mu for the calculations of mu_i
# get_sigma_mu_post = function(sigma_2, sigma_2_mu, V_i)
# {
#   # calculate T_i
#   T_i = nrow(V_i)
#   ones_vector = rep(1, T_i)
#   
#   # calculate sigma_mu_post
#   inner_part_one = crossprod(ones_vector, chol2inv(V_i)) ### !!!!!!!
#   inner = inner_part_one %*% ones_vector + sigma_mu^(-1) # this should not be ^(-2)
#   
#   # invert to return
#   inner_inverted = inner^-1
#   
#   return(inner_inverted)
# }


# # function to get alpha_mu for the calculations of mu_i
# get_alpha_mu_post = function(alpha_mu, sigma_mu, sigma_mu_post, g_i, V_i, y)
# {
#   # set up terms needed
#   y = y[!is.na(y)] # remove any NA values
#   T_i = nrow(V_i)
# 
#   # term_one = alpha_mu^2 * sigma_mu^(-2) # TODO this may be wrong, fixed possibly below
#   term_one = alpha_mu * sigma_mu^(-1) # fixed some things here
# 
#   term_two_part_one = crossprod(y - g_i, chol2inv(V_i))
#   term_two = term_two_part_one %*% rep(1, T_i)
# 
#   to_return = sigma_mu_post * (term_one + term_two)
#   return(to_return)
# }


# # TODO this is where we are no longer sure if we are sane 
# # this is the function to calculate Yi = mu_i - gI
vector_differences = function(y, mu_i, g_i)
{
  Ti = length(y)

  inner = y - rep(mu_i, Ti) - g_i
  return(inner)
}


# # function to calculate acceptance ratio for l_k
lk_acceptance = function(y, mu, g, sigma_2, lk_prime, l_k, time)
{
  # indicator function part
  if (lk_prime < 0.1 || lk_prime > 1 || l_k < 0.1 ||  l_k > 1)
  {
    to_return = 0
  } else { # calcualtions assuming indicator = 1
    # calcualte first term outside of the product
    y_noNA = y[1,][!is.na(y[1,])]
    term_one = y_noNA - rep(mu[1], length(y_noNA)) - g[[1]]

    # get the two v terms necessary
    M_temp = get_matern(l_k, time[[1]])
    M_prime = get_matern(lk_prime, time[[1]])
    
    K_temp = get_K_i(sigma_2, M_temp)
    K_prime = get_K_i(sigma_2, M_prime)
    
    # FIXED - removed sigma_2 superfluous arguments
    V_temp = get_V_i(sigma_2, K_temp) 
    V_prime = get_V_i(sigma_2, K_prime)

    term_two = chol2inv(V_prime) - chol2inv(V_temp)

    # do matrix multiplication
    matrix_part = crossprod(term_one, term_two)
    matrix_part = matrix_part %*% term_one

    # calculate ratio
    ratio = exp(-0.5 * as.numeric(matrix_part))

    for (i in 2:nrow(y))
    {
      y_noNA = y[i,][!is.na(y[i,])]
      
      # calcualte proceeding terms in product
      y_noNA = y[i,][!is.na(y[i,])]
      
      term_one = y_noNA - rep(mu[i], length(y_noNA)) - g[[i]]

      # get two v terms necessary
      M_temp = get_matern(l_k, time[[i]] )
      M_prime = get_matern(lk_prime, time[[i]] )
      
      K_temp = get_K_i(sigma_2, M_temp)
      K_prime = get_K_i(sigma_2, M_prime)
      
      # FIXED - removed sigma_2 superfluous arguments
      V_temp = get_V_i(sigma_2, K_temp) 
      V_prime = get_V_i(sigma_2, K_prime)

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


# # function for MH sampling of lk
get_lk = function(y, mu, g, sigma_2, lk_0, time)
{
  # small epsilon
  lk_t = lk_0

  # step two - draw lk_prime
  lk_prime = rexp(1, lk_t)

  # draw the uniform variable
  u_t = runif(1, 0, 1)

  # calculate the acceptance
  acceptance = lk_acceptance(y, mu, g, sigma_2, lk_prime, lk_t, time)

  # update lk_t1 based on acceptance
  lk_t1 = ifelse(u_t <= acceptance, lk_prime, lk_t)

  # update lk_k (either it will have stayed the same or updated)
  lk_t = lk_t1

  return(lk_t1)
}

# # function for MH sampling of lk
get_lk_c = function(y, mu, g, sigma_2, lk_0, time)
{
  # small epsilon
  lk_t = lk_0
  
  # step two - draw lk_prime
  lk_prime = rexp(1, lk_t)
  
  # draw the uniform variable
  u_t = runif(1, 0, 1)
  
  # calculate the acceptance
  acceptance = lk_acceptance_c(y, mu, g, sigma_2, lk_prime, lk_t, time)
  
  # update lk_t1 based on acceptance
  lk_t1 = ifelse(u_t <= acceptance, lk_prime, lk_t)
  
  # update lk_k (either it will have stayed the same or updated)
  lk_t = lk_t1
  
  return(lk_t1)
}



## LB functions below
lb_acceptance = function(y, lb, lb_prime, xi, knots) # depends on lb, lb', g values (this comes from H matrix and xi), v matrix
{
  # print(knots)
  # indicator function part
  if (lb_prime < 0.1 || lb_prime > 1 || lb < 0.1 ||  lb > 1)
  {
    to_return = 0
  } else { # calcualtions assuming indicator = 1
    # calcualte first term outside of the product
    term_one = xi
    # xi_lb = get_g()
    # xi_lb_prime = get_g()

    # stuff we need to calcualte v_i
    M_lb = get_matern(lb, knots)
    M_lb_prime = get_matern(lb_prime, knots) # this should be dependent on the knot points, not hte lb 

    # use tinv to invert M here
    # term_two = tinv(M_lb) - tinv(M_lb_prime)
    term_two = chol2inv(M_lb) - chol2inv(M_lb_prime)
    

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
get_lb = function(y, lb_0, xi, knots)
{
  lb_t = lb_0

  # step two - draw lb_prime
  lb_prime = rexp(1, lb_t)
  
  #  draw uniform valuable
  u_t = runif(1, 0, 1)
  
  # calculate new acceptance
  acceptance = lb_acceptance(y, lb_t, lb_prime, xi, knots)
  
  # update lb_t1 based on acceptance
  lb_t1 = ifelse(u_t <= acceptance, lb_prime, lb_t)

  return(lb_t1)
}

# function to calculate lb updates
get_lb_c = function(y, lb_0, xi, knots)
{
  lb_t = lb_0
  
  # step two - draw lb_prime
  lb_prime = rexp(1, lb_t)
  
  #  draw uniform valuable
  u_t = runif(1, 0, 1)
  
  # calculate new acceptance
  acceptance = lb_acceptance_c(y, lb_t, lb_prime, xi, knots)
  
  # update lb_t1 based on acceptance
  lb_t1 = ifelse(u_t <= acceptance, lb_prime, lb_t)
  
  return(lb_t1)
}
