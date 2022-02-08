#!/bin/bash

# now run all of the datasets for x = 20 to see how this goes yay
for i in {1..50}; do \
first="--args r="
second=" t=20 B_VAL=300 STORE_VAL=100"
c="${first}$i ${second}"
echo "${c}"
nohup R CMD BATCH --no-save --no-restore "${first}$i ${second}" BLAH.R BLAH.out &
done