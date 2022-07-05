#!/bin/bash
# delete the single folder

Tnum=20
NNGP="FALSE"
n_replicates=50
num_flights=100
max_T=100
MCMCiterations=50000
# rm -rf t$Tnum

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum n_replicates=$n_replicates USE_NNGP=$NNGP num_flights=$num_flights" ANALYZE_simulations_single_.R OUTPUT_analyzing.out &
wait $BACK_PID 

# delete the .rda files 
# rm -r RESULTS
# rm -rf RESULTS/results_n100_t$Tnum*

if [[ $NNGP == "TRUE" ]]; then
echo "DELETING NNGP TRUE FOR $Tnum"
rm -rf RESULTSNNGP/results_n100_t$Tnum*
else
echo "DELETING NNGO FALSE FOR $Tnum"
rm -rf RESULTS/results_n100_t$Tnum*
fi