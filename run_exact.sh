#!/bin/bash

num_observations=10
num_replications=10

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args num_observations=$num_observations num_replications=$num_replications" R/sandy_code/npreg-exact.R exact.out &

