#!/bin/bash
# delete the single folder 
rm -r t$Tnum

# do analysis 
nohup R CMD BATCH --no-save --no-restore "--args t_vals=$Tnum" ANALYZE_simulations_single.R OUTPUT_analyzing.out &
wait $BACK_PID 

# delete the .rda files 
rm -r RESULTS