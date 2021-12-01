rm(list = ls())
library(dplyr)
library(stats)
library(recipes)
library(caret)
folder = "/Users/rachaelshudde/Desktop/Data/AirlineData/"
jan18 = read.csv(paste(folder, "January2018.csv", sep = ""))
feb18 = read.csv(paste(folder, "February2018.csv", sep = ""))
march18 = read.csv(paste(folder, "March2018.csv", sep = ""))
april18 = read.csv(paste(folder, "April2018.csv", sep = ""))
may18 = read.csv(paste(folder, "May2018.csv", sep = ""))
june18 = read.csv(paste(folder, "June2018.csv", sep = ""))
july18 = read.csv(paste(folder, "July2018.csv", sep = ""))

aug = read.csv(paste(folder, "August2019.csv", sep = ""))
sep = read.csv(paste(folder, "September2019.csv", sep = ""))
oct = read.csv(paste(folder, "October2019.csv", sep = ""))
nov = read.csv(paste(folder, "November2019.csv", sep = ""))
dec = read.csv(paste(folder, "December2019.csv", sep = ""))

selected_columns = c("Month", "DayOfWeek", "Marketing_Airline_Network", "FlightDate",
                     "Tail_Number", "Flight_Number_Operating_Airline", "Origin", "Dest", "CRSDepTime",
                     "CRSArrTime", "AirTime", "Distance", "DestState", "OriginState", "DepDelay")

# combined = rbind(jan18, feb18, march18, april18, may18, june18, july18, aug, sep, oct, nov, dec)
combined = rbind(jan18, june18, sep, oct, dec)
rm(jan18)
rm(feb18)
rm(march18)
rm(april18)
rm(may18)
rm(june18)
rm(july18)
rm(aug)
rm(sep)
rm(oct)
rm(nov)
rm(dec)
X_dec = combined[, selected_columns]
idx = which(abs(X_dec$DepDelay) > 120)
X_dec = X_dec[-idx,]

data = X_dec
data = data[complete.cases(data), ] # remove data with missing covariates

# REMOVE THE LESS FREQUENT AIRLINES / DEPARTURES / ORIGINS
origin_sort = as.data.frame(sort(table(data$Origin), decreasing = TRUE))[-c(1:15),]
new_origin_list = origin_sort$Var1
data = data[-(which(data$Origin %in% new_origin_list)), ]

deset_sort = as.data.frame(sort(table(data$Dest), decreasing = TRUE))[-c(1:15),]
new_dest_list = deset_sort$Var1
data = data[-(which(data$Dest %in% new_dest_list)), ]

data = droplevels(data)

# now update arrival / departure time to be factors
data$Month = as.factor(data$Month)
data$DayOfWeek = as.factor(data$DayOfWeek)
data$Tail_Number = NULL # removing this for now due to data size
data$CRSArrTime = data$CRSArrTime/2400
data$CRSDepTime = data$CRSDepTime/2400
max_dist = max(data$Distance)
data$Distance = data$Distance/max(data$Distance) # this can stay continuous right now


# set up Y variables
columns = c("Flight_Number_Operating_Airline", "DepDelay", "FlightDate", "Marketing_Airline_Network", "Origin", "Dest")
Y = data[, columns]
colnames(Y) = c("num", "delay", "date", "airline", "origin", "dest")
Y = Y[order(as.Date(Y$date, format="%Y-%m-%d")),]
Y$date = as.Date(Y$date, format="%Y-%m-%d")

Y = Y %>% group_by(num, date, airline, origin, dest) %>% summarize(delay = mean(delay, na.rm = TRUE))
Y = data.frame(Y)
head(Y)


# take log and remove any delays over 2 hours 
adjustment = abs(min(Y$delay)) + 1
temp = Y$delay + adjustment
log_mean = mean(log(temp))
Y$delay = log(temp) - log_mean

# rm(data)

# get minimum and miximum for range
min = min(Y$date)
max = max(Y$date)
range = length(numeric(max-min)) + 1

# set up range of dates
temp_year = as.numeric(format(min, format = "%Y"))
temp_month = as.numeric(format(min, format = "%m"))
# temp_day = format(dates, format = "%Y")
dates_list = format( seq(c(ISOdate(temp_year,temp_month,1)), by = "DSTday", length.out = range),"%Y-%m-%d")
flights = unique(Y$num)
date_range = data.frame(dates_list)
colnames(date_range) = "date"
date_range$date = as.Date(date_range$date, format="%Y-%m-%d")

# get length - add origin / destination 
Y$unique_flights = paste(Y$airline, Y$num, Y$origin, Y$dest, sep = "")

# initialize the overall dataframe
Z = matrix(NA, nrow = length(unique(Y$unique_flights)), ncol = length(dates_list))
colnames(Z) = dates_list

# get the Z
count = 1
for (i in unique(Y$unique_flights))
{
  #print(paste("Attempt?ing to figure out for: ", i))
  temp = Y[which(Y$unique_flights == i), ]
  a = merge(temp, date_range, by = "date", all.y = TRUE)
  Z[count, ] = a$delay
  count = count + 1
  
  if (count %% 200 == 0) print(i)
}

rownames(Z) = unique(Y$unique_flights)
# save(Z,file = "Z_matrix_redone.Rda")

print(dim(Z))
percent_missing = apply(Z, 1, function(row) sum(is.na(row)) / length(row)*100)
hist(percent_missing, main = "Percent of days a flight is not flown", xlab = "Percent Missing")
summary(percent_missing)

# remove flights which didn't fly
over_missing_30 = which(percent_missing > 30)
over_missing_50 = which(percent_missing > 50) # reasonable 
over_missing_75 = which(percent_missing > 75)
over_missing_80 = which(percent_missing > 85)
over_missing_90 = which(percent_missing > 90)

print(length(over_missing_30))
print(length(over_missing_50))
print(length(over_missing_75))
print(length(over_missing_80))
print(length(over_missing_90))

old_Z = Z
Z = Z[-over_missing_75, ]
unique_flights = unique(rownames(Z))

# if a flight takes over 2 hours, assume it is "not observed"
# Z[Z > 120] = NA

X_dec = data # may need to update
X_dec$ID = paste(X_dec$Marketing_Airline_Network, X_dec$Flight_Number_Operating_Airline, X_dec$Origin, X_dec$Dest, sep = "")
rm(data)
# ## try to get rid of srme of this
# na_count <- rowSums(is.na(Z))
# table(na_count)

X_list = list()
count = 1
c_X = vector(length = length(unique_flights))
for (i in unique_flights)
{
  temp = X_dec[X_dec$ID == i,]
  temp$date = temp$FlightDate
  temp$FlightDate = NULL
  temp$date = as.Date(temp$date, format="%Y-%m-%d")
  temp = merge(temp, date_range, by = "date", all.y = TRUE)
  temp$DepDelay = NULL
  temp$AirTime = NULL
  temp$ID = NULL
  rownames(temp) = temp$date
  temp$date = NULL
  temp$Flight_Number_Operating_Airline = NULL
  # temp[1:10, 1:ncol(temp)]
  
  temp = temp[complete.cases(temp), ]
  dmy <- dummyVars(" ~ .", data = temp)
  done <- data.frame(predict(dmy, newdata = temp))
  #continuous 
  continuous = c("CRSArrTime", "CRSDepTime", "Distance")
  temp_cont = temp[, continuous]
  # do index matrices
  temp2 = temp[, -c(6:8)]
  rec_obj = recipe(~., data = temp2) %>% step_dummy(all_nominal_predictors())
  prep_step = prep(rec_obj, temp2)
  done = bake(prep_step, temp2)
  
  done = as.matrix(done)
  rownames(done) = rownames(temp)
  done = cbind(done, temp_cont)
  done = as.matrix(done)

  # get c_max 
  c_X[i] = max(apply(X = done, MARGIN = 1,
                     FUN = function(r){sqrt(sum(r^2))}))
  
  X_list[[count]] = done
  count = count + 1
  
  if (count %% 100 == 0) print(count)
}

c.X = max(c_X)
for (i in 1:length(X_list))
{
  X_list[[i]] = X_list[[i]]/c.X
}
  
# keep distance continuous 
# Z transform


source('R/FUNC_woodchan_samples.R')
source('R/FUNC_paramater_estimates.R')
source('R/DATA_generate_simulation.R')
source('R/FUNC_Gibbs_Sampler.R')
source('R/FUNC_Gibbs_Sampler_r.R')
source('R/PLOTS_Gibbs_Sampler.R')
Rcpp::sourceCpp("src/FUNC_paramater_estimates_c.cpp")

beta_names = colnames(X_list[[1]])
time_idx = apply(Z, 1, function(x) which(!is.na(x)))
start = 1
end = 20
data_actual = list(y = Z[start:end,], X = X_list[start:end], time_idx = time_idx[start:end])
# data_gibbs = data_actual; B = 1000; n_to_store = 500
results = gibbs_sampler_r(data_gibbs = data_actual, B = 10000, n_to_store = 5000)

save(results)
# get posterior g samples
# g_s = list()
# post_g = colMeans(results$g)
# last = 1
# for (i in 1:end)
# {
#   temp = Z[i,]
#   if (length(which(is.na(temp))) > 0) temp = temp[-which(is.na(temp))]
#   current = length(temp)
#   g_s[[i]] = post_g[last:(last + current - 1)]
#   print(paste("for", i, " - ", last, ":", last + current - 1, "(", length(temp), ")","(", length(last:(last + current - 1)), ")"))
#   last = last + current - 1
# }

# #plots
par(mfrow = c(3,3))
plot(colMeans(results$beta), main = "Beta plots")
plot(results$sigma_2, type = "l", main = "sigma_2 plot")
plot(results$lB, type  = "l", main = "lb")
plot(colMeans(results$mu), main = "mu plots")
plot(results$loglhood, type = "l", main = "loglikelihood")
plot(colMeans(results$xi), main = "xi results")
plot(results$beta[,5], type = "l")
plot(results$beta[,25], type = "l")
plot(results$beta[,45], type = "l")
par(mfrow = c(1,1)); plot(colMeans(results$g))


betas = as.data.frame(beta_names)
betas$beta_values = colMeans(round(results$beta, 3))
betas = betas[order(abs(betas$beta_values), decreasing = TRUE),]
betas

# look at prediction values of 
idx = sample(1:end, 1)
Z2 = Z[idx,]
Z2 = Z2[-which(is.na(Z2))]
Z2_pred = colMeans(results$mu)[idx] + g_s[[idx]]
pred = as.data.frame(cbind(Z2, Z2_pred))
colnames(pred) = c("actual", "predicted")
# pred$actual = exp(pred$actual) - adjustment
# pred$predicted = exp(pred$predicted)
par(mfrow = c(1,2)); plot(pred$actual, pred$predicted); plot(abs(pred$actual - pred$predicted))


# 
