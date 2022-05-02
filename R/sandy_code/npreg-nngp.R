rm(list = ls())
source('R/sandy_code/npreg.R')
library(doParallel)

args=(commandArgs(TRUE))
if(length(args)==0){
  stop("NO COMMAND LINE ARGUMENTS PASSED")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

R = num_replications
nObs = num_observations
mcmc = num_mcmc
burnin = num_burnin

## A monotone function (as in Maatouk & Bay)
f1 = function(x){log(20*x + 1)}
x.seq = seq(0, 1, by = .01)
f.true = f1(x.seq)

set.seed(1)
xObs = sort(runif(nObs, 0, 1))
f1_x = f1(xObs)
sig.true = .2

# f2 = function(x){5*((x-0.5)^2)}

# generating the data
doParallel::registerDoParallel(cores = 18)
NNGP100 = foreach::foreach(r = 1:R, .combine = 'rbind', .multicombine = T) %dopar% {
  start = proc.time()
  
  set.seed(r)
  yObs = f1_x + rnorm(nObs, 0, sig.true)
  
  MCMCout = NNGP.sampler(y = yObs, x = xObs, nu.fix = 5/2, l.fix = 1,
                         print.at = 1000, mcmc = mcmc, brn = burnin)
  
  # summaries
  MBmat = t(mapply(i = 1:length(x.seq),
                   FUN = function(i){
                     
                     pmax(1 - abs((x.seq[i] - MCMCout$knots)*MCMCout$nknots), 0)
                     
                   }))
  fpost = MBmat%*%MCMCout$xi
  fpost.summ = rbind(rowMeans(fpost),
                     apply(fpost, 1,
                           function(v){quantile(x = v, probs = c(0.025, 0.975))}))
  fpost.coverage = as.numeric((fpost.summ[2,]<=f.true)&(fpost.summ[3,]>=f.true))
  
  print(r)
  end = proc.time()
  time = as.numeric(end[3] - start[3])
  c(sqrt(mean((fpost.summ[1,]-f.true)^2)), mean(fpost.coverage), mean(MCMCout$time.per.iter), time = time)
}

print(head(NNGP100))
save(NNGP100, file = 'TIMING/NNGP100.RData')

