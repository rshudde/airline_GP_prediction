rm(list = ls())
library(dplyr)
library(stats)
library(recipes)
library(caret)
source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_sampler.R')
Rcpp::sourceCpp('src/FUNC_paramater_estimates_c.cpp')

# file
load("/Users/rachaelshudde/Desktop/REAL_DATA_GIBBS.rda")
print(results$time)

time_in_seconds = results$time
time_in_minutes = time_in_seconds / 60
time_in_hours = time_in_minutes / 60

print(paste("Time in Minutes:", round(time_in_minutes,2), "| Time in hours:", round(time_in_hours,2)))

num_iterations = 30000 / 10
total_time_in_minutes = num_iterations * time_in_minutes
total_time_in_hours = num_iterations * time_in_hours
total_time_in_days = total_time_in_hours / 24

print(paste("Total expected time in days: ", round(total_time_in_days, 1)))
