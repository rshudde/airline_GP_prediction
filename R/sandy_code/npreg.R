rm(list = ls())
# source codes and libraries:
require(Rcpp)
library(matrixcalc)
source("R/sandy_code/all_functions.R")
# sourceCpp("inv_chol.cpp")
Rcpp::sourceCpp('src/inv_chol.cpp')

exactGP.sampler = function(y, x, N, mcmc = 1000, brn = 4000, thin = 1, 
                           nu.in, l.in, tau2.in, sig2.in, xi.in,
                           nu.fix, l.fix, tau2.fix, sig2.fix, xi.fix,
                           verbose = T, print.at = 500){
  
  if(length(y)!=length(x)) stop("y and x should be of same length!")
  n = length(y)
  if(missing(N)) N = floor(n/2)
  my_knots = (0:N)/N
  X = t(mapply(i = 1:n,
               FUN = function(i){
                 
                 pmax(1 - abs((x[i] - my_knots)*N), 0)
                 
               }))
  em = thin*mcmc + brn
  ef = mcmc
  
  if(!missing(nu.fix)) nu.in = nu.fix
  if(!missing(l.fix)) l.in = l.fix
  if(!missing(tau2.fix)) tau2.in = tau2.fix
  if(!missing(sig2.fix)) sig2.in = sig2.fix
  if(!missing(xi.fix)) xi.in = xi.fix
  
  if(missing(nu.fix) && missing(nu.in)) nu.in = 1
  if(missing(l.fix) && missing(l.in)) l.in = l_est(nu.in, c(0,1), 0.05)
  if(missing(tau2.fix) && missing(tau2.in)) tau2.in = 1
  if(missing(sig2.fix) && missing(sig2.in)) sig2.in = 1
  
  Kmat = round(covmat(my_knots, nu.in, l.in), 5)
  Kmat = Kmat + diag(10e-3, nrow = nrow(Kmat))
  Linv = round(inv_chol(Kmat), 5)
         
  if(missing(xi.fix) && missing(xi.in)){
    
    xi.in = as.numeric(mvtnorm::rmvnorm(n = 1, mean = numeric(N+1),
                                        sigma = tau2.in*Kmat))
  }
  
  range.nu = c(0.5, 1)
  range.l = c(0.1, 1)
  sd.proposal.nu = 0.05
  sd.proposal.l = 0.1
  
  nu_sam = l_sam = tau2_sam = sig2_sam = time.per.iter = rep(NA, ef)
  xi_sam = matrix(nrow = N+1, ncol = ef)
  fhat_sam = matrix(nrow = n, ncol = ef)
  
  frac.seq = rep(1, em)
  # frac.seq[1:brn] = exp(seq(log(0.5), log(1), length.out = brn))
  
  for(i in 1:em){
    set.seed(i)
    
    ptm.start_i = proc.time()
    
    # sampling from nu or l
    if(!missing(nu.fix) & missing(l.fix)){
      print('1')
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      if(dproposed.l>0){
        
        Kmat.proposed = round(covmat(my_knots, nu.in, l.proposed), 5)
        Linv.proposed = round(inv_chol(Kmat.proposed), 5)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (l.proposed/l.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
      
    }else if(missing(nu.fix) & !missing(l.fix)){
      print('2')
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if(dproposed.nu>0){
        
        Kmat.proposed = round(covmat(my_knots, nu.proposed, l.in), 5)
        Linv.proposed = round(inv_chol(Kmat.proposed), 5)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          nu.in = nu.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
      
    }else if(missing(nu.fix) & missing(l.fix)){
      print('3')
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if((dproposed.l>0)||(dproposed.nu>0)){
        
        if(dproposed.l==0) l.proposed = l.in
        if(dproposed.nu==0) nu.proposed = nu.in
        
        Kmat.proposed = round(covmat(my_knots, nu.proposed, l.proposed), 5)
        Linv.proposed = round(inv_chol(Kmat.proposed), 5)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (l.proposed/l.in)*(nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          nu.in = nu.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
    }
    
    # sampling \tau^2:
    if(missing(tau2.fix)) tau2.in = invgamma::rinvgamma(n = 1, shape = (N+1)/2,
                                                        rate = sum((Matrix::crossprod(Linv, xi.in))^2)/2)
    # sampling \sigma^2:
    if(missing(sig2.fix)) sig2.in = invgamma::rinvgamma(n = 1, shape = n/2,
                                                        rate = sum((y - as.numeric(X%*%xi.in))^2)/2)
    
    # sampling Xi:
    if(missing(xi.fix)){
    
      # Sigma.xipost = round(solve(Matrix::crossprod(X)/sig2.in +
      #                              Matrix::tcrossprod(Linv)/tau2.in), 5)
      a = round(Matrix::crossprod(X)/sig2.in, 5)
      b = round(Matrix::tcrossprod(Linv)/tau2.in, 5)
      c = round(a + b, 5)
      Sigma.xipost = (solve(c) + t(solve(c)))/2
      mu.xipost = as.numeric(Sigma.xipost%*%(Matrix::crossprod(X, y)/sig2.in))
      xi.in = as.numeric(mvtnorm::rmvnorm(n = 1, mean = mu.xipost, sigma = Sigma.xipost))
    }
    
    ptm.end_i = proc.time()
    
    # storing MCMC samples:
    if(i > brn && i%%thin == 0){
      
      nu_sam[(i-brn)/thin] = nu.in
      l_sam[(i-brn)/thin] = l.in
      tau2_sam[(i-brn)/thin] = tau2.in
      sig2_sam[(i-brn)/thin] = sig2.in
      xi_sam[,(i-brn)/thin] = xi.in
      fhat_sam[,(i-brn)/thin] = as.numeric(X%*%xi.in)
      time.per.iter[(i-brn)/thin] = as.numeric((ptm.end_i-ptm.start_i)[3])
    }
    
    if(i%%print.at==0 && verbose) print(i)
  }
  
  f.postsumm = apply(fhat_sam, 1,
                     function(x){c(mean(x),
                                   quantile(x, c(0.025, 0.975)))})
  
  return(list('nu' = nu_sam, 'l' = l_sam, 'tau2' = tau2_sam,
              'sig2' = sig2_sam, 'knots' = my_knots, 'nknots' = N+1, 'xi' = xi_sam, 
              'fhat' = fhat_sam, 'fmean' = f.postsumm[1,],
              'flow' = f.postsumm[2,], 'fup' = f.postsumm[3,],
              'time.per.iter' = time.per.iter))
}

WCGP.sampler = function(y, x, N, mcmc = 1000, brn = 4000, thin = 1, 
                        nu.in, l.in, tau2.in, sig2.in, xi.in,
                        nu.fix, l.fix, tau2.fix, sig2.fix, xi.fix,
                        verbose = T, print.at = 500){
  
  if(length(y)!=length(x)) stop("y and x should be of same length!")
  n = length(y)
  if(missing(N)) N = floor(n/2)
  my_knots = (0:N)/N
  X = t(mapply(i = 1:n,
               FUN = function(i){
                 
                 pmax(1 - abs((x[i] - my_knots)*N), 0)
                 
               }))
  em = thin*mcmc + brn
  ef = mcmc
  
  if(!missing(nu.fix)) nu.in = nu.fix
  if(!missing(l.fix)) l.in = l.fix
  if(!missing(tau2.fix)) tau2.in = tau2.fix
  if(!missing(sig2.fix)) sig2.in = sig2.fix
  if(!missing(xi.fix)) xi.in = xi.fix
  
  if(missing(nu.fix) && missing(nu.in)) nu.in = 1
  if(missing(l.fix) && missing(l.in)) l.in = l_est(nu.in, c(0,1), 0.05)
  if(missing(tau2.fix) && missing(tau2.in)) tau2.in = 1
  if(missing(sig2.fix) && missing(sig2.in)) sig2.in = 1
  
  Kmat = covmat(my_knots, nu.in, l.in)
  Linv = inv_chol(Kmat)
  
  if(missing(xi.fix) && missing(xi.in)){
    
    xi.in = as.numeric(mvtnorm::rmvnorm(n = 1, mean = numeric(N+1),
                                        sigma = tau2.in*Kmat))
  }
  
  range.nu = c(0.5, 1)
  range.l = c(0.1, 1)
  sd.proposal.nu = 0.05
  sd.proposal.l = 0.1
  
  nu_sam = l_sam = tau2_sam = sig2_sam = time.per.iter = rep(NA, ef)
  xi_sam = matrix(nrow = N+1, ncol = ef)
  fhat_sam = matrix(nrow = n, ncol = ef)
  
  frac.seq = rep(1, em)
  # frac.seq[1:brn] = exp(seq(log(0.5), log(1), length.out = brn))
  
  for(i in 1:em){
    
    set.seed(i)
    
    ptm.start_i = proc.time()
    
    # sampling from nu or l
    if(!missing(nu.fix) & missing(l.fix)){
      
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      if(dproposed.l>0){
        
        Kmat.proposed = covmat(my_knots, nu.in, l.proposed)
        Linv.proposed = inv_chol(Kmat.proposed)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (l.proposed/l.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
      
    }else if(missing(nu.fix) & !missing(l.fix)){
      
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if(dproposed.nu>0){
        
        Kmat.proposed = covmat(my_knots, nu.proposed, l.in)
        Linv.proposed = inv_chol(Kmat.proposed)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          nu.in = nu.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
      
    }else if(missing(nu.fix) & missing(l.fix)){
      
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if((dproposed.l>0)||(dproposed.nu>0)){
        
        if(dproposed.l==0) l.proposed = l.in
        if(dproposed.nu==0) nu.proposed = nu.in
        
        Kmat.proposed = covmat(my_knots, nu.proposed, l.proposed)
        Linv.proposed = inv_chol(Kmat.proposed)
        acceptance.ratio =
          min(1, exp(sum(log(diag(Linv.proposed))) - 
                       sum((Matrix::crossprod(Linv.proposed, xi.in))^2)/(2*tau2.in) -
                       sum(log(diag(Linv))) +
                       sum((Matrix::crossprod(Linv, xi.in))^2)/(2*tau2.in))*
                (l.proposed/l.in)*(nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          nu.in = nu.proposed
          Kmat = Kmat.proposed
          Linv = Linv.proposed
        }
      }
    }
    
    # sampling \tau^2:
    if(missing(tau2.fix)) tau2.in = invgamma::rinvgamma(n = 1, shape = (N+1)/2,
                                                        rate = sum((Matrix::crossprod(Linv, xi.in))^2)/2)
    
    # sampling \sigma^2:
    if(missing(sig2.fix)) sig2.in = invgamma::rinvgamma(n = 1, shape = n/2,
                                                        rate = sum((y - as.numeric(X%*%xi.in))^2)/2)
    
    # sampling Xi:
    if(missing(xi.fix)){
      
      # step one
      theta = runif(1, 0, 2*pi)
      xi.prior = as.numeric(samp.WC(my_knots, nu.in, l.in, tau2.in))
      xi.proposed = cos(theta)*xi.in + sin(theta)*xi.prior
      
      # step two
      theta_min = theta - 2 * pi
      theta_max = theta
      
      # old negative loglikelihood
      negloglhood.old = (frac.seq[i]*sum((y-as.numeric(X%*%xi.in))^2))/(2*sig2.in)
      
      # new negative loglikelihood
      negloglhood.new = (frac.seq[i]*sum((y-as.numeric(X%*%xi.proposed))^2))/(2*sig2.in)
      
      # calculate new acceptance value
      acceptance = min(1, exp(negloglhood.old - negloglhood.new))
      
      # step 3
      zeta = runif(1, 0, 1)
      
      # continuation of step 3 - don't return until we get something we accept 
      while (acceptance <= zeta ){
        
        # step a
        if(theta < 0){
          
          theta_min = theta
          
        }else{theta_max = theta}
        
        # step b 
        theta = runif(1, theta_min, theta_max)
        
        # step c
        xi.proposed = cos(theta)*xi.in + sin(theta)*xi.prior
        
        # new negative loglikelihood
        negloglhood.new = (frac.seq[i]*sum((y-as.numeric(X%*%xi.proposed))^2))/(2*sig2.in)
        
        # calculate new acceptance value
        acceptance = min(1, exp(negloglhood.old - negloglhood.new))
      }
      
      xi.in = xi.proposed
    }
    
    ptm.end_i = proc.time()
    
    # storing MCMC samples:
    if(i > brn && i%%thin == 0){
      
      nu_sam[(i-brn)/thin] = nu.in
      l_sam[(i-brn)/thin] = l.in
      tau2_sam[(i-brn)/thin] = tau2.in
      sig2_sam[(i-brn)/thin] = sig2.in
      xi_sam[,(i-brn)/thin] = xi.in
      fhat_sam[,(i-brn)/thin] = as.numeric(X%*%xi.in)
      time.per.iter[(i-brn)/thin] = as.numeric((ptm.end_i-ptm.start_i)[3])
    }
    
    if(i%%print.at==0 && verbose) print(i)
  }
  
  f.postsumm = apply(fhat_sam, 1,
                     function(x){c(mean(x),
                                   quantile(x, c(0.025, 0.975)))})
  
  return(list('nu' = nu_sam, 'l' = l_sam, 'tau2' = tau2_sam,
              'sig2' = sig2_sam, 'knots' = my_knots, 'nknots' = N+1, 'xi' = xi_sam, 
              'fhat' = fhat_sam, 'fmean' = f.postsumm[1,],
              'flow' = f.postsumm[2,], 'fup' = f.postsumm[3,],
              'time.per.iter' = time.per.iter))
}

NNGP.sampler = function(y, x, N, nNeighbour = 30, mcmc = 1000, brn = 4000, thin = 1, 
                        nu.in, l.in, tau2.in, sig2.in, xi.in,
                        nu.fix, l.fix, tau2.fix, sig2.fix, xi.fix,
                        verbose = T, print.at = 500){
  
  if(length(y)!=length(x)) stop("y and x should be of same length!")
  n = length(y)
  if(missing(N)) N = floor(n/2)
  my_knots = (0:N)/N
  X = t(mapply(i = 1:n,
               FUN = function(i){
                 
                 pmax(1 - abs((x[i] - my_knots)*N), 0)
                 
               }))
  em = thin*mcmc + brn
  ef = mcmc
  
  if(!missing(nu.fix)) nu.in = nu.fix
  if(!missing(l.fix)) l.in = l.fix
  if(!missing(tau2.fix)) tau2.in = tau2.fix
  if(!missing(sig2.fix)) sig2.in = sig2.fix
  if(!missing(xi.fix)) xi.in = xi.fix
  
  if(missing(nu.fix) && missing(nu.in)) nu.in = 1
  if(missing(l.fix) && missing(l.in)) l.in = l_est(nu.in, c(0,1), 0.05)
  if(missing(tau2.fix) && missing(tau2.in)) tau2.in = 1
  if(missing(sig2.fix) && missing(sig2.in)) sig2.in = 1
  
  # K matrix and some related terms
  Kout.list = lapply(X = 1:(N+1),
                     FUN = function(X){
                       
                       if(X==1){
                         
                         Linv_X = NA
                         v_X = NA
                         a_X = 0
                         d_X = 1
                         
                       }else if(X==2){
                         
                         Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                         nu.in, l.in)
                         Linv_X = 1
                         v_X = Kmat_X[1,2]
                         a_X = Kmat_X[1,2]
                         d_X = 1 - a_X*Kmat_X[1,2]
                         
                       }else{
                         
                         Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                         nu.in, l.in)
                         Kmat_X = Kmat_X + diag(10e-5, nrow = nrow(Kmat_X))
                      
                         Linv_X = inv_chol(Kmat_X[1:min(nNeighbour,X-1),
                                                  1:min(nNeighbour,X-1)])
                         v_X = Matrix::crossprod(Linv_X, Kmat_X[1:min(nNeighbour,X-1),
                                                                min(nNeighbour,X-1)+1])
                         a_X = as.numeric(Matrix::crossprod(t(Linv_X), v_X))
                         d_X = 1 - sum(v_X^2)
                       }
                       
                       list('Linv' = Linv_X, 'v' = v_X, 'a' = a_X, 'd' = d_X)
                     })
  
  if(missing(xi.fix) && missing(xi.in)){
    
    # initializing xi from the NNGP
    xi.in = rep(NA, N+1)
    xi.in[1] = rnorm(1, 0, sqrt(tau2.in))
    for(j in 2:(N+1)){
      
      xi.in[j] = sum(Kout.list[[j]]$a*xi.in[max(1,j-nNeighbour):(j-1)]) +
        rnorm(1, 0, sqrt(tau2.in*Kout.list[[j]]$d))
    }
  }
  
  range.nu = c(0.5, 1)
  range.l = c(0.1, 1)
  sd.proposal.nu = 0.05
  sd.proposal.l = 0.1
  
  nu_sam = l_sam = tau2_sam = sig2_sam = time.per.iter = rep(NA, ef)
  xi_sam = matrix(nrow = N+1, ncol = ef)
  fhat_sam = matrix(nrow = n, ncol = ef)
  
  frac.seq = rep(1, em)
  # frac.seq[1:brn] = exp(seq(log(0.5), log(1), length.out = brn))
  
  for(i in 1:em){
    
    set.seed(i)
    
    ptm.start_i = proc.time()
    
    # sampling from nu or l using MH
    if(!missing(nu.fix) & missing(l.fix)){
      
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      if(dproposed.l>0){
        
        # new K matrix and some related terms
        Kout.list.proposed = lapply(X = 1:(N+1),
                                    FUN = function(X){
                                      
                                      if(X==1){
                                        
                                        Linv_X = NA
                                        v_X = NA
                                        a_X = 0
                                        d_X = 1
                                        
                                      }else if(X==2){
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.in, l.proposed)
                                        Linv_X = 1
                                        v_X = Kmat_X[1,2]
                                        a_X = Kmat_X[1,2]
                                        d_X = 1 - a_X*Kmat_X[1,2]
                                        
                                      }else{
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.in, l.proposed)
                                        Linv_X = inv_chol(Kmat_X[1:min(nNeighbour,X-1),
                                                                 1:min(nNeighbour,X-1)])
                                        v_X = Matrix::crossprod(Linv_X, Kmat_X[1:min(nNeighbour,X-1),
                                                                               min(nNeighbour,X-1)+1])
                                        a_X = as.numeric(Matrix::crossprod(t(Linv_X), v_X))
                                        d_X = 1 - sum(v_X^2)
                                      }
                                      
                                      list('Linv' = Linv_X, 'v' = v_X, 'a' = a_X, 'd' = d_X)
                                    })
        
        acceptance.ratio =
          min(1, exp(Reduce('+',
                            lapply(X = 2:(N+1),
                                   FUN = function(X){
                                     
                                     dnorm(xi.in[X], 
                                           sum(Kout.list.proposed[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                           sqrt(tau2.in*Kout.list.proposed[[X]]$d), log = T) -
                                       dnorm(xi.in[X], 
                                             sum(Kout.list[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                             sqrt(tau2.in*Kout.list[[X]]$d), log = T)
                                     
                                   })))*(l.proposed/l.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          Kout.list = Kout.list.proposed
        }
      }
      
    }else if(missing(nu.fix) & !missing(l.fix)){
      
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if(dproposed.nu>0){
        
        # new K matrix and some related terms
        Kout.list.proposed = lapply(X = 1:(N+1),
                                    FUN = function(X){
                                      
                                      if(X==1){
                                        
                                        Linv_X = NA
                                        v_X = NA
                                        a_X = 0
                                        d_X = 1
                                        
                                      }else if(X==2){
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.proposed, l.in)
                                        Linv_X = 1
                                        v_X = Kmat_X[1,2]
                                        a_X = Kmat_X[1,2]
                                        d_X = 1 - a_X*Kmat_X[1,2]
                                        
                                      }else{
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.proposed, l.in)
                                        Linv_X = inv_chol(Kmat_X[1:min(nNeighbour,X-1),
                                                                 1:min(nNeighbour,X-1)])
                                        v_X = Matrix::crossprod(Linv_X, Kmat_X[1:min(nNeighbour,X-1),
                                                                               min(nNeighbour,X-1)+1])
                                        a_X = as.numeric(Matrix::crossprod(t(Linv_X), v_X))
                                        d_X = 1 - sum(v_X^2)
                                      }
                                      
                                      list('Linv' = Linv_X, 'v' = v_X, 'a' = a_X, 'd' = d_X)
                                    })
        
        acceptance.ratio =
          min(1, exp(Reduce('+',
                            lapply(X = 2:(N+1),
                                   FUN = function(X){
                                     
                                     dnorm(xi.in[X], 
                                           sum(Kout.list.proposed[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                           sqrt(tau2.in*Kout.list.proposed[[X]]$d), log = T) -
                                       dnorm(xi.in[X], 
                                             sum(Kout.list[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                             sqrt(tau2.in*Kout.list[[X]]$d), log = T)
                                     
                                   })))*(nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          nu.in = nu.proposed
          Kout.list = Kout.list.proposed
        }
      }
      
    }else if(missing(nu.fix) & missing(l.fix)){
      
      l.proposed = exp(log(l.in) + rnorm(1, 0, sd.proposal.l))
      nu.proposed = exp(log(nu.in) + rnorm(1, 0, sd.proposal.nu))
      dproposed.l = dunif(l.proposed, range.l[1], range.l[2])
      dproposed.nu = dunif(nu.proposed, range.nu[1], range.nu[2])
      if((dproposed.l>0)||(dproposed.nu>0)){
        
        if(dproposed.l==0) l.proposed = l.in
        if(dproposed.nu==0) nu.proposed = nu.in
        
        # new K matrix and some related terms
        Kout.list.proposed = lapply(X = 1:(N+1),
                                    FUN = function(X){
                                      
                                      if(X==1){
                                        
                                        Linv_X = NA
                                        v_X = NA
                                        a_X = 0
                                        d_X = 1
                                        
                                      }else if(X==2){
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.proposed, l.proposed)
                                        Linv_X = 1
                                        v_X = Kmat_X[1,2]
                                        a_X = Kmat_X[1,2]
                                        d_X = 1 - a_X*Kmat_X[1,2]
                                        
                                      }else{
                                        
                                        Kmat_X = covmat(my_knots[max(1,X-nNeighbour):X],
                                                        nu.proposed, l.proposed)
                                        Linv_X = inv_chol(Kmat_X[1:min(nNeighbour,X-1),
                                                                 1:min(nNeighbour,X-1)])
                                        v_X = Matrix::crossprod(Linv_X, Kmat_X[1:min(nNeighbour,X-1),
                                                                               min(nNeighbour,X-1)+1])
                                        a_X = as.numeric(Matrix::crossprod(t(Linv_X), v_X))
                                        d_X = 1 - sum(v_X^2)
                                      }
                                      
                                      list('Linv' = Linv_X, 'v' = v_X, 'a' = a_X, 'd' = d_X)
                                    })
        
        acceptance.ratio =
          min(1, exp(Reduce('+',
                            lapply(X = 2:(N+1),
                                   FUN = function(X){
                                     
                                     dnorm(xi.in[X], 
                                           sum(Kout.list.proposed[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                           sqrt(tau2.in*Kout.list.proposed[[X]]$d), log = T) -
                                       dnorm(xi.in[X], 
                                             sum(Kout.list[[X]]$a*xi.in[max(1,X-nNeighbour):(X-1)]),
                                             sqrt(tau2.in*Kout.list[[X]]$d), log = T)
                                     
                                   })))*(l.proposed/l.in)*(nu.proposed/nu.in))
        
        if(runif(1)<acceptance.ratio){
          
          l.in = l.proposed
          nu.in = nu.proposed
          Kout.list = Kout.list.proposed
        }
      }
    }
    
    # sampling \tau^2:
    if(missing(tau2.fix)) tau2.in = invgamma::rinvgamma(n = 1, shape = (N+1)/2,
                                                        rate = (xi.in[1]^2 +
                                                                  Reduce('+',
                                                                         lapply(X = 2:(N+1),
                                                                                FUN = function(X){
                                                                                  
                                                                                  ((xi.in[X] -
                                                                                      sum(Kout.list[[X]]$a*
                                                                                            xi.in[max(1,X-nNeighbour):(X-1)]))^2)/
                                                                                    Kout.list[[X]]$d
                                                                                })))/2)
    
    # sampling \sigma^2:
    if(missing(sig2.fix)) sig2.in = invgamma::rinvgamma(n = 1, shape = n/2,
                                                        rate = (frac.seq[i]*sum((y - as.numeric(X%*%xi.in))^2))/2)
    
    # sampling Xi:
    if(missing(xi.fix)){
      
      # draw from the NNGP prior
      theta = runif(1, 0, 2*pi)
      if(i==1) xi.prior = rep(NA, N+1)
      xi.prior[1] = rnorm(1, 0, sqrt(tau2.in))
      for(j in 2:(N+1)){
        
        xi.prior[j] = sum(Kout.list[[j]]$a*xi.prior[max(1,j-nNeighbour):(j-1)]) +
          rnorm(1, 0, sqrt(tau2.in*Kout.list[[j]]$d))
      }
      xi.proposed = cos(theta)*xi.in + sin(theta)*xi.prior
      
      # step two
      theta_min = theta - 2 * pi
      theta_max = theta
      
      # old negative loglikelihood
      negloglhood.old = (frac.seq[i]*sum((y-as.numeric(X%*%xi.in))^2))/(2*sig2.in)
      
      # new negative loglikelihood
      negloglhood.new = (frac.seq[i]*sum((y-as.numeric(X%*%xi.proposed))^2))/(2*sig2.in)
      
      # calculate new acceptance value
      acceptance = min(1, exp(negloglhood.old - negloglhood.new))
      
      # step 3
      zeta = runif(1, 0, 1)
      
      # continuation of step 3 - don't return until we get something we accept 
      while (acceptance <= zeta){
        
        # step a
        if(theta < 0){
          
          theta_min = theta
          
        }else{theta_max = theta}
        
        # step b 
        theta = runif(1, theta_min, theta_max)
        
        # step c
        xi.proposed = cos(theta)*xi.in + sin(theta)*xi.prior
        
        # new negative loglikelihood
        negloglhood.new = (frac.seq[i]*sum((y-as.numeric(X%*%xi.proposed))^2))/(2*sig2.in)
        
        # calculate new acceptance value
        acceptance = min(1, exp(negloglhood.old - negloglhood.new))
      }
      
      xi.in = xi.proposed
    }
    
    ptm.end_i = proc.time()
    
    # storing MCMC samples:
    if(i > brn && i%%thin == 0){
      
      nu_sam[(i-brn)/thin] = nu.in
      l_sam[(i-brn)/thin] = l.in
      tau2_sam[(i-brn)/thin] = tau2.in
      sig2_sam[(i-brn)/thin] = sig2.in
      xi_sam[,(i-brn)/thin] = xi.in
      fhat_sam[,(i-brn)/thin] = as.numeric(X%*%xi.in)
      time.per.iter[(i-brn)/thin] = as.numeric((ptm.end_i-ptm.start_i)[3])
    }
    
    if(i%%print.at==0 && verbose) print(i)
  }
  
  f.postsumm = apply(fhat_sam, 1,
                     function(x){c(mean(x),
                                   quantile(x, c(0.025, 0.975)))})
  
  return(list('nu' = nu_sam, 'l' = l_sam, 'tau2' = tau2_sam,
              'sig2' = sig2_sam, 'knots' = my_knots, 'nknots' = N+1, 'xi' = xi_sam, 
              'fhat' = fhat_sam, 'fmean' = f.postsumm[1,],
              'flow' = f.postsumm[2,], 'fup' = f.postsumm[3,],
              'time.per.iter' = time.per.iter))
}

