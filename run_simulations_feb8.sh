#!/bin/bash
mkdir -p t20
mkdir -p t40
mkdir -p t60
mkdir -p t80
mkdir -p RESULTS
mkdir -p output

# now create the data
nohup R CMD BATCH R/DATA_subset_datasets.R R/DATA_subset_datasets.out &

# now run all of the datasets for x = 20 to see how this goes yay
for i in {1..3}; do \
first="--args r="
second=" t=20"
c="${first}$i ${second}"
echo "${c}"
nohup R CMD BATCH --no-save --no-restore "${first}$i ${second}" BLAH.R BLAH.out &
done

