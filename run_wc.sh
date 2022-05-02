#!/bin/bash

num_observations=100
num_replications=20
num_burnin=100
num_mcmc=100

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args num_observations=$num_observations num_replications=$num_replications num_burnin=$num_burnin num_mcmc=$num_mcmc" R/sandy_code/npreg-wc.R wc.out &

