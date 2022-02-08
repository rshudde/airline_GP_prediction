#!/bin/bash
mkdir -p t20
mkdir -p t40
mkdir -p t60
mkdir -p t80
mkdir -p RESULTS
mkdir -p output
echo "made directories"

# now create the data
nohup R CMD BATCH R/DATA_subset_datasets.R R/DATA_subset_datasets.out &
