# This is a file to start doing paramater estimates for sampling mu
library(lubridate)

# functin to do Gibbs sampling 
gibbs_sampler_r = function(data_gibbs, B = 1000, 
                         mu_initial, beta_initial, sigma_2_initial, lK_initial,
                         xi_initial, sigmaB_2_initial, lB_initial, nNeighbour,
                         a = 10^(-3), b = 10^(-3), alpha_normal_prior = 0,
                         sigma_normal_prior = 1000, burn_in = 0.5, cpp = TRUE, 
                         NNGP = FALSE, n_to_store = 100)
{
  # # for testing
  # B = 1000
  # a = 10^(-3)
  # b = 10^(-3)
  # alpha_normal_prior = 0
  # sigma_normal_prior = 1000
  # burn_in = 0.5
  # sigma_2_initial = 1
  # sigmaB_2_initial = 1
  # mu_initial = data$mu_true
  # beta_initial = data$beta_true
  # sigma_2_initial = data$sigma_2_true
  # xi_initial = runif(length(data$xi_true), -1, 1)
  # lK_initial = data$lK_true
  # lB_initial = data$lB_true
  # NNGP = FALSE
  # n_to_store = 100

  #### data ####
  # get X and y values from the data
  X = data_gibbs$X
  y = data_gibbs$y
  if (is.null(data_gibbs$time_idx))
  {
    time_idx = list()
    for(i in 1:length(X)) time_idx[[i]] = c(1:nrow(X[[i]]))
    
  } else {
    time_idx = data_gibbs$time_idx
  }
  
  # get the number of datasets, covariates and knots
  n_datasets = length(X)
  n_covariates = ncol(X[[1]])
  n_nonNA_y = sum(unlist(lapply(time_idx, length)))
  
  if (missing(xi_initial))
  {
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
  mu_post = matrix(nrow = n_to_store, ncol = n_datasets)
  alpha_post = beta_post = matrix(nrow = n_to_store, ncol = n_covariates)
  xi_post = matrix(nrow = n_to_store, ncol = n_Knots_gibbs)
  sigmaB_2_post = sigma_2_post = lK_post = lB_post = rep(NA, n_to_store)
  w_post = g_post = matrix(nrow = n_to_store, ncol = n_nonNA_y)
  loglhood_gibbs = numeric(n_to_store)
  loglhood_gibbs_current = 0
  
  #### initialize chain ####
  if (missing(mu_initial)) mu_initial = numeric(n_datasets)
  if (missing(beta_initial))
  {
    beta_initial = rep(1, n_covariates)
    beta_initial = beta_initial/sqrt(sum(beta_initial^2))
  }
  if (missing(sigma_2_initial)) sigma_2_initial = 1
  if (missing(sigmaB_2_initial)) sigmaB_2_initial = 1
  if (missing(xi_initial)) xi_initial = runif(n_Knots_gibbs, -5, 5) #rep(1, n_Knots_gibbs)
  if (missing(lK_initial)) lK_initial = 0.5
  if (missing(lB_initial)) lB_initial = 1/length(X)
  
  ##  initializing for testing
  # mu_initial = numeric(n_datasets)
  # beta_initial = rep(1, n_covariates)
  # beta_initial = beta_initial/sqrt(sum(beta_initial^2))
  # sigma_2_initial = 1
  # sigmaB_2_initial = 1
  # xi_initial = runif(n_Knots_gibbs, -5, 5) #rep(1, n_Knots_gibbs)
  # lK_initial = 0.5
  # lB_initial = 0.5

  
  mu_post_current = mu_initial
  alpha_post_current = beta_initial
  beta_post_current = beta_initial
  sigma_2_post_current = sigma_2_initial
  xi_post_current = xi_initial # length of knots
  sigmaB_2_post_current = sigmaB_2_initial
  lK_post_current = lK_initial
  lB_post_current = lB_initial
  
  ## updating parameter related quantities
  # M, K, V, H matrix
  M_gibbs = K_gibbs = V_gibbs = w_gibbs = H_gibbs = g_gibbs = list()
  for (i in 1:n_datasets)
  {
    M_gibbs[[i]] = get_matern(lK_post_current, time_idx[[i]])
    K_gibbs[[i]] = get_K_i(sigma_2_post_current, M_gibbs[[i]])
    V_gibbs[[i]] = get_V_i(sigma_2_post_current, K_gibbs[[i]])
    w_gibbs[[i]] = (as.numeric(X[[i]] %*% beta_post_current) + 1)/2
    H_gibbs[[i]] = get_H_matrix(w_gibbs[[i]], knots_gibbs, 
                                n_Knots_gibbs)
    g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post_current)
  }
  
  
  ################################################################################################################
  ################################################################################################################
  idx = 2
  count = 1
  cutoff = B - n_to_store + 1
  loglikelihoodtest = vector(length = B + 1)
  start_outer = proc.time()
  start_inner = proc.time()
  #### starting Gibbs ####
  for (idx in 2:(B + 1))
  {
    set.seed(idx)
  
    #### getting mu ####
    mu_post_current = get_mu_c(y, n_datasets, g_gibbs, V_gibbs, time_idx,
                              sigma_2_mu_gibbs, alpha_mu_gibbs)
    mu_post_current[1] = 0
    mu_post_current = as.vector(mu_post_current)
    
    if (idx > cutoff)
    {
      mu_post[count, ] = mu_post_current
    }
    
    # mu_post[idx, ] = data_gibbs$mu_true
    
    # #### getting beta ####
    alpha_gibbs_out = get_alpha_c(alpha_post_current, y, n_datasets, time_idx, X,
                                  n_covariates, mu_post_current, xi_post_current,
                                  V_gibbs, knots_gibbs, n_Knots_gibbs,
                                  sigma_normal_prior)
    alpha_post_current = alpha_gibbs_out$alpha
    beta_post_current = alpha_gibbs_out$beta
    
    if (idx > cutoff)
    {
      alpha_post[count, ] = alpha_gibbs_out$alpha
      beta_post[count, ] = alpha_gibbs_out$beta
    }
    
    # updating beta related term
    w_gibbs = alpha_gibbs_out$w # list of vectors of length n_covairates
    H_gibbs = alpha_gibbs_out$H_mat # list of matrices of length n_datasets, dim n_covariates x n_Knots_gibbs
    g_gibbs = alpha_gibbs_out$g # list of vectors of length n_covairates
    
    # beta_post[idx, ] = data_gibbs$beta_true
    for (i in 1:n_datasets)
    {
      w_gibbs[[i]] = (as.numeric(X[[i]] %*% beta_post_current) + 1)/2
      H_gibbs[[i]] = get_H_matrix(w_gibbs[[i]], knots_gibbs,
                                  n_Knots_gibbs)
      g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post_current)
    }
    
    
    #### getting sigma_2 ####
    print(paste("THIS IS: ", sigma_2_post_current))
    sigma_2_gibbs_out = get_sigma_2_c(a_gibbs, b_gibbs, y, n_datasets, n_nonNA_y,
                                      time_idx, mu_post_current, M_gibbs, g_gibbs)
    sigma_2_post_current = sigma_2_gibbs_out$sigma_2
    
    if (idx > cutoff)
    {
      sigma_2_post[count] = sigma_2_gibbs_out$sigma_2
    }
    
    # sigma_2_post[idx] = data_gibbs$sigma_2_true
    
    # updating sigma_2 related term
    for (i in 1:n_datasets)
    {
      K_gibbs[[i]] = get_K_i(sigma_2_post_current, M_gibbs[[i]])
      V_gibbs[[i]] = get_V_i(sigma_2_post_current, K_gibbs[[i]])
    }
    
    
    #### get l_k ####
    save(y, file="y.rda")
    save(mu_post_current, file="m.rda")
    save(g_gibbs, file="g_gibbs.rda")
    save(sigma_2_post_current, file="sigma.rda")
    save(lK_post_current, file = "lk.rda")
    save(time, file = "time.rda")
    lK_post_current = get_lk(y = y, mu = mu_post_current, g = g_gibbs, sigma_2 = sigma_2_post_current, 
                               lk_0 = lK_post_current, time = time_idx)
    
    # y = y; mu = mu_post_current; g = g_gibbs; sigma_2 = sigma_2_post_current; lk_0 = lK_post_current; time = time_idx
    
    if (idx > cutoff)
    {
      lK_post[count] = lK_post_current
    }
    
    
    # lK_post[idx] = data_gibbs$lK_true
    
    # updating l_k related term
    for (i in 1:n_datasets)
    {
      M_gibbs[[i]] = get_matern(lK_post_current, time_idx[[i]]) # only c function is fastest here
      K_gibbs[[i]] = get_K_i(sigma_2_post_current, M_gibbs[[i]])
      V_gibbs[[i]] = get_V_i(sigma_2_post_current, K_gibbs[[i]])
    }
    
    #### get l_b ####
    lB_post_current = get_lb(y, lB_post_current, xi_post_current, knots_gibbs)
    # lB_post[idx] = data_gibbs$lB_true
    
    if (idx > cutoff)
    {
      lB_post[count] = lB_post_current
    }
    
    
    #### getting sigmaB_2 ####
    sigmaB_2_post_temp = get_sigmaB_2(a_gibbs, b_gibbs, xi_post_current, # CHANGED FOMR _c
                                        lB_post_current, knots_gibbs, n_Knots_gibbs)
    if (!is.na(sigmaB_2_post_temp)) sigmaB_2_post_current = sigmaB_2_post_temp
    
    if (idx > cutoff)
    {
      sigmaB_2_post[count] = sigmaB_2_post_current
    }
    
    
    #### getting xi ####
    xi_gibbs_out = get_xi_c(xi_post_current, sigmaB_2_post_current, y, n_datasets,
                            time_idx, mu_post_current, H_gibbs, V_gibbs,
                            lB_post_current, knots_gibbs, nNeighbour, NNGP)
    xi_post_current = xi_gibbs_out$xi

    
    if (idx > cutoff)
    {
      xi_post[count,] = xi_post_current
    }
    
    # updating xi related term
    g_gibbs = xi_gibbs_out$g
    
    # for debugging
    # xi_post[idx, ] = data_gibbs$xi_true
    # for (i in 1:n_datasets)
    # {
    #   g_gibbs[[i]] = get_g(H_gibbs[[i]], xi_post[idx, ])
    # }
    
    
    #### get log likelihood ####
    loglhood_gibbs_current = 0
    for (i in 1:n_datasets)
    {
      loglhood_gibbs_current = loglhood_gibbs_current - as.numeric(determinant(V_gibbs[[i]], log = T)$modulus)/2 -
        (emulator::quad.form.inv(V_gibbs[[i]], y[i,time_idx[[i]]] - mu_post_current[i] - g_gibbs[[i]]))/2
    }
    # loglhood_gibbs[idx] = loglhood_gibbs[idx] - xi_gibbs_out$negloglhood
    
    #### g ####
    w_post_current = unlist(w_gibbs)
    g_post_current = unlist(g_gibbs)
    
    if (idx > cutoff)
    {
      loglhood_gibbs[count] = loglhood_gibbs_current
      w_post[count, ] = w_post_current
      g_post[count, ] = g_post_current
      count = count + 1
    }
    
    loglikelihoodtest[idx] = loglhood_gibbs_current
    # print(paste("GOT TO:", idx))
    # print statement for time  

    if (idx %% 100 == 0) 
    {
      end_inner = proc.time()
      time_inner = round(as.numeric(end_inner[3] - start_inner[3]),2)
      print(paste("iteration:", idx, "in", seconds_to_period(time_inner)))
      start_inner = proc.time()
    }
      
      
    print(paste(idx, "|", lK_post_current))
  }
  # print(round(Sys.time() - start),2)
  end_outer = proc.time()
  
  #### get burnin to remove #### 
  burn_in = floor(B * burn_in)
  thinning = seq(1, n_to_store, by = 2)
  time = as.numeric(end_outer[3] - start_outer[3])
  print(paste("TOTAL TIME", seconds_to_period(time)))
  
  return(list(beta = beta_post[thinning, ],
              mu = mu_post[thinning, ],
              sigma_2 = sigma_2_post[thinning],
              xi = xi_post[thinning, ],
              w = w_post[thinning, ],
              g = g_post[thinning, ],
              sigmaB_2 = sigmaB_2_post[thinning],
              # lB = lB_post[thinning],
              lB = lB_post,
              lK = lK_post[thinning],
              loglhood = loglhood_gibbs[thinning],
              knots = knots_gibbs,
              LOGTEST = loglikelihoodtest[-1], 
              time = time))
}





