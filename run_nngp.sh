#!/bin/bash

num_observations=100
num_replications=20

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args num_observations=$num_observations num_replications=$num_replications" R/sandy_code/npreg-nngp.R nngp.out &

