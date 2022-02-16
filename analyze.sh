#!/bin/bash
# delete the single folder
Tnum=10
n_replicates=50
rm -r t$Tnum

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum n_replicates=$n_replicates" ANALYZE_simulations_single.R OUTPUT_analyzing.out &
wait $BACK_PID 

# delete the .rda files 
# rm -r RESULTS
rm -rf RESULTS/results_n100_t$Tnum*