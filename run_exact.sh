#!/bin/bash

num_observations=400
num_replications=20
num_burnin=$num_observations-$num_observations/4
num_mcmc=5000

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args num_observations=$num_observations num_replications=$num_replications num_burnin=$num_burnin num_mcmc=$num_mcmc" R/sandy_code/npreg-exact.R exact.out &

