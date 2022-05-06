#!/bin/bash

num_observations=800
num_replication=20
num_burnin=500
num_mcmc=1000

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args num_observations=$num_observations num_replications=$num_replications num_burnin=$num_burnin num_mcmc=$num_mcmc" R/sandy_code/npreg-nngp.R nngp.out &

