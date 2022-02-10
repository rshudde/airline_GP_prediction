#!/bin/bash

Tnum=20
MCMCiterations=200
nsave=50

# create the single folder to store results as well as the output and results files
mkdir -p t$Tnum
mkdir -p RESULTS
mkdir -p output

# populate the folder if it's not empty
if ! [ "$(ls -A t$Tnum)" ]; then
	nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum" R/DATA_subset_datasets.R OUTPUT_creating_data.out &
fi
wait $BACK_PID 

# run the Gibbs samplers (loop over 50)
for i in {1..50}; do \
nohup R CMD BATCH --no-save --no-restore "--args r=$i t=$Tnum B_VAL=$MCMCiterations STORE_VAL=$nsave" BLAH.R OUTPUT_$Tnum.out &
echo "--args r=$i t=$Tnum B_VAL=$MCMCiterations STORE_VAL=$nsave"
done
wait $BACK_PID 

# delete the single folder 
rm -r t$Tnum

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum" ANALYZE_simulations_single.R OUTPUT_analyzing.out &
wait $BACK_PID 

# delete the .rda files 
rm -r RESULTS

