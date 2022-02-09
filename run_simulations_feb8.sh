#!/bin/bash

# now run all of the datasets for x = 20 to see how this goes yay
for i in {1..25}; do \
first="--args r="
second=" t=60 B_VAL=10000 STORE_VAL=1000"
c="${first}$i ${second}"
echo "${c}"
nohup R CMD BATCH --no-save --no-restore "${first}$i ${second}" BLAH.R BLAH60.out &
done