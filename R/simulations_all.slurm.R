#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=simulations_output     #Set the job name to "JobExample2"
#SBATCH --time=4:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=28                   #Request 1 node
#SBATCH --ntasks=1         #Request 8 tasks/cores per node
#SBATCH --ntasks-per-node=1         #Request 8 tasks/cores per node
#SBATCH --mem=32G                     #Request 8GB per node 
#SBATCH --output=output/simulations_output.out      #Send stdout/err to "Example2Out.[jiobID]" 
#SBATCH --array=1-50

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --account=132744323415             #Set billing account to 123456
##SBATCH --mail-type=ALL              #Send email on all job events
##SBATCH --mail-user=ananya.rc94@tamu.edu    #Send all emails to email_address 
#SBATCH --mail-type=END,FAIL,BEGIN

module load GCC/10.3.0
module load R/4.1.0-foss-2021a
R CMD BATCH --no-save --no-restore simulations_20.R simulations_20_$SLURM_ARRAY_TASK_ID

