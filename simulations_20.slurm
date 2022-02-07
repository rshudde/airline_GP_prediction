#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=simulations_output     #Set the job name to "JobExample2"
#SBATCH --time=1:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=28                   #Request 1 node
#SBATCH --ntasks=1         #Request 8 tasks/cores per node
#SBATCH --ntasks-per-node=1         #Request 8 tasks/cores per node
#SBATCH --mem=32G                     #Request 8GB per node 
#SBATCH --output=simulations_output.out      #Send stdout/err to "Example2Out.[jobID]" 
#SBATCH --array=1-4

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --account=132752207705             #Set billing account to 123456
##SBATCH --mail-type=ALL              #Send email on all job events
##SBATCH --mail-user=rachaelshudde@tamu.edu    #Send all emails to email_address 
#SBATCH --mail-type=END,FAIL,BEGIN

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2
R CMD BATCH --no-save --no-restore simulations_20.R simulations_20_$SLURM_ARRAY_TASK_ID

