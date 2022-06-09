rm(list = ls())
library(dplyr)
library(stats)
library(recipes)
library(caret)

should_write = FALSE
if (should_write) {
  folder = "/Users/rachaelshudde/Desktop/Data/AirlineData/"
  
  months = c("January", "February", "March", "April", "May", "June", "July", "August", 
             "September", "October", "November", "December")
  
  combined = c()
  for (k in months) {
    title = paste(folder, k, "2018.csv", sep = "")
    temp = read.csv(title)
    print(k)
    combined = rbind(combined, temp)
    rm(temp)
  }
  
  save(combined, file = "combineddata.rda")
} else 
{
  print("LOADING FILE")
  load("combineddata.rda")
}


selected_columns = c("Month", "DayOfWeek", "Marketing_Airline_Network", "FlightDate",
                     "Tail_Number", "Flight_Number_Operating_Airline", "Origin", "Dest", "CRSDepTime",
                     "CRSArrTime", "AirTime", "Distance", "DestState", "OriginState", "DepDelay")

################# clean data
X_dec = combined[, selected_columns]
idx = which(abs(X_dec$DepDelay) > 120)
print(paste("Initial data is of dimension ", nrow(X_dec), "rows and", ncol(X_dec), "columns")) 
X_dec = X_dec[-idx,]
print(paste("Data after removing large delays is of dimension ", nrow(X_dec), "rows and", ncol(X_dec), 
            "columns, removing", length(idx), "rows")) 

data = X_dec
data = data[complete.cases(data), ] # remove data with missing covariates
print(paste("Data after removing incomplete cases is of dimension ", nrow(data), "rows and", ncol(data), 
            "columns, removing", nrow(X_dec) - nrow(data), "rows")) 
temp = nrow(data)
# REMOVE THE LESS FREQUENT AIRLINES / DEPARTURES / ORIGINS
origin_sort = as.data.frame(sort(table(data$Origin), decreasing = TRUE))[-c(1:15),]
new_origin_list = origin_sort$Var1
data = data[-(which(data$Origin %in% new_origin_list)), ]

deset_sort = as.data.frame(sort(table(data$Dest), decreasing = TRUE))[-c(1:15),]
new_dest_list = deset_sort$Var1
data = data[-(which(data$Dest %in% new_dest_list)), ]

data = droplevels(data)
print(paste("Data after removing infrequent airlines cases is of dimension ", nrow(data), "rows and", ncol(data), 
            "columns, removing", temp - nrow(data), "rows")) 

################# finish clean data

# now update arrival / departure time to be factors
data$Month = as.factor(data$Month)
data$DayOfWeek = as.factor(data$DayOfWeek)
data$Tail_Number = NULL # removing this for now due to data size
data$CRSArrTime = data$CRSArrTime/2400
data$CRSDepTime = data$CRSDepTime/2400
max_dist = max(data$Distance)
data$Distance = data$Distance/max(data$Distance) # this can stay continuous right now


################# set up Y variables
columns = c("Flight_Number_Operating_Airline", "DepDelay", "FlightDate", "Marketing_Airline_Network", "Origin", "Dest")
Y = data[, columns]
colnames(Y) = c("num", "delay", "date", "airline", "origin", "dest")
Y = Y[order(as.Date(Y$date, format="%Y-%m-%d")),]
Y$date = as.Date(Y$date, format="%Y-%m-%d")

Y = Y %>% group_by(num, date, airline, origin, dest) %>% summarize(delay = mean(delay, na.rm = TRUE))
Y = data.frame(Y)
head(Y)


#################  take log of delay and remove any delays over 2 hours 
adjustment = abs(min(Y$delay)) + 1
temp = Y$delay + adjustment
log_mean = mean(log(temp))
Y$delay = log(temp) - log_mean

# rm(data)

# get minimum and miximum for range
min = min(Y$date)
max = max(Y$date)
range = length(numeric(max-min)) + 1

################# set up range of dates
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

################# get the Z matrix
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

print(paste("Number of flights missing 30% of observations", length(over_missing_30)))
print(paste("Number of flights missing 50% of observations", length(over_missing_50)))
print(paste("Number of flights missing 75% of observations", length(over_missing_75)))
print(paste("Number of flights missing 80% of observations", length(over_missing_80)))
print(paste("Number of flights missing 90% of observations", length(over_missing_90)))

old_Z = Z
Z = Z[-over_missing_75, ]
unique_flights = unique(rownames(Z))

# if a flight takes over 2 hours, assume it is "not observed"
# Z[Z > 120] = NA

X_dec = data # may need to update
X_dec$ID = paste(X_dec$Marketing_Airline_Network, X_dec$Flight_Number_Operating_Airline, X_dec$Origin, X_dec$Dest, sep = "")
# rm(data)
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
  
save(X_list, file = "X_list.rda")
save(Z, file = "Z_list.rda")
# keep distance continuous 
# Z transform

