#!/bin/bash

# now run all of the datasets for x = 20 to see how this goes yay
for i in {1..30}; do \
first="--args r="
second=" t=40 B_VAL=50000 STORE_VAL=20000"
c="${first}$i ${second}"
echo "${c}"
nohup R CMD BATCH --no-save --no-restore "${first}$i ${second}" BLAH.R BLAH40.out &
done