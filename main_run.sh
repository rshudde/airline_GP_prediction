#!/bin/bash

Tnum=20
NNGP="FALSE"
n_replicates=10
MCMCiterations=4000
nsave=2000

# create the single folder to store results as well as the output and results files
mkdir -p t$Tnum
mkdir -p RESULTS
mkdir -p output
mkdir -p outputNNGP
mkdir -p RESULTSNNGP

# populate the folder if it's not empty
if ! [ "$(ls -A t$Tnum)" ]; then
	nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum n_replicates=$n_replicates" R/DATA_subset_datasets.R OUTPUT_creating_data.out &
fi
wait $BACK_PID 

# run the Gibbs samplers (loop over 50)
for i in $(seq 1 $n_replicates); do \
nohup R CMD BATCH --no-save --no-restore "--args r=$i t=$Tnum B_VAL=$MCMCiterations STORE_VAL=$nsave USE_NNGP=$NNGP" main_R_file.R BLAHOUTPUT_$Tnum$NNGP.out &
echo "--args r=$i t=$Tnum B_VAL=$MCMCiterations STORE_VAL=$nsave USE_NNGP=$NNGP"
done



