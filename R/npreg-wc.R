
source('npreg.R')
library(doParallel)

## A monotone function (as in Maatouk & Bay)
f1 = function(x){log(20*x + 1)}
x.seq = seq(0, 1, by = .01)
f.true = f1(x.seq)

set.seed(1)
nObs = 100
xObs = sort(runif(nObs, 0, 1))
f1_x = f1(xObs)
sig.true = .2

# f2 = function(x){5*((x-0.5)^2)}

# generating the data
R = 10
doParallel::registerDoParallel(cores = 18)
WCGP100 = foreach::foreach(r = 1:R, .combine = 'rbind', .multicombine = T) %dopar% {
  
  set.seed(r)
  yObs = f1_x + rnorm(nObs, 0, sig.true)
  
  MCMCout = WCGP.sampler(y = yObs, x = xObs, nu.fix = 5/2, l.fix = 1,
                            print.at = 100, mcmc = 100, brn = 400
  )
  
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
  
  c(sqrt(mean((fpost.summ[1,]-f.true)^2)), mean(fpost.coverage), mean(MCMCout$time.per.iter))
}

save(WCGP100, file = 'WCGP100.RData')

