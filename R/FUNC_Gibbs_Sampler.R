# This is a file to start doing paramater estimates for sampling mu
library(invgamma)
library(MASS)
library(nlme)

# functin to do Gibbs sampling 
gibbs_sampler = function(data_gibbs, B = 1000, 
                         mu_initial, beta_initial, sigma_2_initial, lK_initial,
                         xi_initial, sigmaB_2_initial, lB_initial,
                         a = 10^(-3), b = 10^(-3), alpha_normal_prior = 0,
                         sigma_normal_prior = 1000, burn_in = 0.5)
{
  #### data ####
  # get X and y values from the data
  X = data_gibbs$X
  y = data_gibbs$y
  time_idx = data_gibbs$time_idx
  
  # get the number of datasets, covariates and knots
  n_datasets = length(X)
  n_covariates = ncol(X[[1]])
  n_nonNA_y = 0
  for(i in 1:n_datasets)
  {
    n_nonNA_y = n_nonNA_y + length(time_idx[[i]])
  }
  
  if(missing(xi_initial)){
    
    n_Knots_gibbs = ceiling(n_nonNA_y/2) + 1
    
  }else{n_Knots_gibbs = length(xi_initial)}
  # n_Knots_gibbs = ifelse(missing(xi_initial), 
  #                        ceiling(n_nonNA_y/2) + 1,
  #                        length(xi_initial))
  knots_gibbs = seq(0, 1, length.out = n_Knots_gibbs)
  
  #### fixed hyperparamaters ####
  sigma_2_mu_gibbs = sigma_normal_prior^2
  alpha_mu_gibbs = alpha_normal_prior
  a_gibbs = a
  b_gibbs = b
  
  #### storage ####
  mu_post = matrix(nrow = B+1, ncol = n_datasets)
  alpha_post = beta_post = matrix(nrow = B+1, ncol = n_covariates)
  xi_post = matrix(nrow = B+1, ncol = n_Knots_gibbs)
  sigmaB_2_post = sigma_2_post = lK_post = lB_post = rep(NA, B+1)
  w_post = g_post = matrix(nrow = B+1, ncol = n_nonNA_y)
  loglhood_gibbs = numeric(B+1)
  
  #### initialize chain ####
  if(missing(mu_initial)) mu_initial = numeric(n_datasets)
  if(missing(beta_initial)){
    
    beta_initial = rep(1, n_covariates)
    beta_initial = beta_initial/sqrt(sum(beta_initial^2))
  }
  if(missing(sigma_2_initial)) sigma_2_initial = 1
  if(missing(sigmaB_2_initial)) sigmaB_2_initial = 1
  if(missing(xi_initial)) xi_initial = runif(n_Knots_gibbs, -5, 5) #rep(1, n_Knots_gibbs)
  if(missing(lK_initial)) lK_initial = 0.1
  if(missing(lB_initial)) lB_initial = 0.1
  
  # initializing
  mu_post[1, ] = mu_initial
  alpha_post[1, ] = beta_initial
  beta_post[1, ] = beta_initial
  sigma_2_post[1] = sigma_2_initial
  xi_post[1, ] = xi_initial # length of knots
  sigmaB_2_post[1] = sigmaB_2_initial
  lK_post[1] = lK_initial
  lB_post[1] = lB_initial
  
  
  ## updating parameter related quantities
  # M, K, V, H matrix
  M_gibbs = K_gibbs = V_gibbs = w_gibbs = H_gibbs = g_gibbs = list()
  for (i in 1:n_datasets)
  {
    M_gibbs[[i]] = get_matern(lK_post[1], time_idx[[i]])
    K_gibbs[[i]] = get_K_i(sigma_2_post[1], M_gibbs[[i]])
    V_gibbs[[i]] = get_V_i(sigma_2_post[1], K_gibbs[[i]])
    w_gibbs[[i]] = (as.numeric(X[[i]] %*% beta_post[1, ]) + 1)/2
    H_gibbs[[i]] = get_H_matrix(w_gibbs[[i]], knots_gibbs, 
                                n_Knots_gibbs)
    g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post[1, ])
  }
  
  
  ################################################################################################################
  ################################################################################################################
  
  #### starting Gibbs ####
  start = Sys.time()
  for (idx in 2:(B+1))
  {
    start_inner = Sys.time()
    set.seed(idx)
    
    
    #### getting mu ####
    mu_post[idx, ] = get_mu(y, n_datasets, g_gibbs, V_gibbs, time_idx,
                            sigma_2_mu_gibbs, alpha_mu_gibbs)
    mu_post[idx, 1] = 0
    
    # mu_post[idx, ] = data_gibbs$mu_true
    
    
    #### getting beta ####
    alpha_gibbs_out = get_alpha(alpha_post[idx-1, ], y, n_datasets, time_idx, X,
                                n_covariates, mu_post[idx, ], xi_post[idx-1, ],
                                V_gibbs, knots_gibbs, n_Knots_gibbs,
                                sigma_normal_prior)
    alpha_post[idx, ] = alpha_gibbs_out$alpha
    beta_post[idx, ] = alpha_gibbs_out$beta

    # updating beta related term
    w_gibbs = alpha_gibbs_out$w
    H_gibbs = alpha_gibbs_out$H_mat
    g_gibbs = alpha_gibbs_out$g
    
    # beta_post[idx, ] = data_gibbs$beta_true
    # for (i in 1:n_datasets)
    # {
    #   w_gibbs[[i]] = (as.numeric(X[[i]] %*% beta_post[idx, ]) + 1)/2
    #   H_gibbs[[i]] = get_H_matrix(w_gibbs[[i]], knots_gibbs,
    #                               n_Knots_gibbs)
    #   g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post[idx - 1, ])
    # }

    
    #### getting sigma_2 ####
    sigma_2_gibbs_out = get_sigma_2(a_gibbs, b_gibbs, y, n_datasets, n_nonNA_y,
                                    time_idx, mu_post[idx, ], M_gibbs, g_gibbs)
    sigma_2_post[idx] = sigma_2_gibbs_out$sigma_2
    
    # sigma_2_post[idx] = data_gibbs$sigma_2_true
    
    # updating sigma_2 related term
    for (i in 1:n_datasets)
    {
      K_gibbs[[i]] = get_K_i(sigma_2_post[idx], M_gibbs[[i]])
      V_gibbs[[i]] = get_V_i(sigma_2_post[idx], K_gibbs[[i]])
    }
    
    
    #### get l_k ####
    lK_post[idx] = get_lk(y, mu_post[idx, ], g_gibbs, sigma_2_post[idx - 1], lK_post[idx-1], time_idx)
    
    # lK_post[idx] = data_gibbs$lK_true
    
    # updating l_k related term
    for (i in 1:n_datasets)
    {
      M_gibbs[[i]] = get_matern(lK_post[idx], time_idx[[i]])
      K_gibbs[[i]] = get_K_i(sigma_2_post[idx], M_gibbs[[i]])
      V_gibbs[[i]] = get_V_i(sigma_2_post[idx], K_gibbs[[i]])
    }
    
    
    #### get l_b ####
    lB_post[idx] = get_lb(y, lB_post[idx-1], xi_post[idx-1, ])
    # lB_post[idx] = data_gibbs$lB_true
    
    
    #### getting sigmaB_2 ####
    sigmaB_2_post[idx] = get_sigmaB_2(a_gibbs, b_gibbs, xi_post[idx - 1, ],
                                      lB_post[idx], knots_gibbs, n_Knots_gibbs)
    
    
    #### getting xi ####
    xi_gibbs_out = get_xi(xi_post[idx - 1, ], sigmaB_2_post[idx], y, n_datasets,
                          time_idx, mu_post[idx, ], H_gibbs, V_gibbs,
                          lB_post[idx], knots_gibbs)
    xi_post[idx, ] = xi_gibbs_out$xi
    
    # updating xi related term
    g_gibbs = xi_gibbs_out$g
    
    # xi_post[idx, ] = data_gibbs$xi_true
    # for (i in 1:n_datasets)
    # {
    #   g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post[idx, ])
    # }
    
    
    #### get log likelihood ####
    for(i in 1:n_datasets)
    {
      loglhood_gibbs[idx] = loglhood_gibbs[idx] -
        as.numeric(determinant(V_gibbs[[i]], log = T)$modulus)/2 -
        (emulator::quad.form.inv(V_gibbs[[i]],
                                 y[i,time_idx[[i]]] - mu_post[idx, i] -
                                   g_gibbs[[i]]))/2
    }
    # loglhood_gibbs[idx] = loglhood_gibbs[idx] - xi_gibbs_out$negloglhood
    
    
    #### g ####
    w_post[idx, ] = unlist(w_gibbs)
    g_post[idx, ] = unlist(g_gibbs)
    
    
    # print statement for time  
    if (idx %% 1000 == 0) print(paste("iteration:", idx, "in", round(Sys.time() - start_inner, 2)))
  }
  # print(round(Sys.time() - start),2)
  
  
  #### get burnin to remove #### 
  burn_in = floor(B * burn_in)
  
  return(list(beta = beta_post[(burn_in+2):(B+1), ], 
              mu = mu_post[(burn_in+2):(B+1), ],
              sigma_2 = sigma_2_post, #sigma_2_post[(burn_in+2):(B+1)], 
              xi = xi_post[(burn_in+2):(B+1), ],
              w = w_post[(burn_in+2):(B+1), ],
              g = g_post[(burn_in+2):(B+1), ],
              sigmaB_2 = sigmaB_2_post[(burn_in+2):(B+1)],
              lB = lB_post[(burn_in+2):(B+1)], 
              lK = lK_post[(burn_in+2):(B+1)],
              loglhood = loglhood_gibbs[(burn_in+2):(B+1)]))
}





