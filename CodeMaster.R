################################################################################
# Nonparametric Take Home Final Code
################################################################################

################################################################################
# Library Importation
################################################################################
library(tidyverse)
library(data.table)
library(stats)
library(pspline)
library(mapproj)
library(lubridate)
library(ASSIST)
library(TSA)
library(nlme)
library(bigsplines)
library(splines)
library(mgcv)
library(ggplot2)
library(forecast)

################################################################################
# Data Importation
################################################################################

data_dir <- paste0(getwd(), "/archive")

stations_file <- paste0(data_dir, "/dot_traffic_stations_2015.txt")

stations <- read.csv(stations_file, sep=',')
stations$fips_county_code <- sapply(stations$fips_county_code, function(x) {str_pad(x, 3, side='left', pad="0")})
stations$fips_state_code <- sapply(stations$fips_state_code, function(x) {str_pad(x, 2, side='left', pad="0")})

df <- read.csv(paste0(data_dir, "/trimmed.csv"), sep=',')
df$station_id <- str_pad(df$station_id, 6, side='left', pad='0')

################################################################################
# Data Cleaning and Preparation
################################################################################

pdx <- '003011' #Get a station in the Portland metro area

#Filtering for chosen monitoring station ids
# df <- df[df$station_id %in% c(beach, pdx), ]

# #Filter dataframe to desired columns, create hour columns
# hours <- paste0("traffic_volume_counted_after_", 
# 	str_pad(str_pad(0:23, 2, side='left', pad='0'), 4, side='right', pad='0'), 
# 	"_to_", 
# 	str_pad(str_pad(1:24, 2, side='left', pad='0'), 4, side='right', pad='0'))

# #Get list of columns to filter the original dataframe by
# filter_cols <- c(c("X...date", "direction_of_travel", 
# 	'direction_of_travel_name', 'fips_state_code', 'functional_classification',
# 	'functional_classification_name', 'lane_of_travel', 'record_type',
# 	'restrictions', 'station_id'), hours)
# #Filter the dataframe
# df <- df[, filter_cols]

# #Create a vector of new column names
# new_cols <- c('date', 'travel_direction', 'travel_direction_name',
# 	'FIPS_state', 'classif', 'classif_name', 'travel_lane',
# 	'record_type', 'restrictions', 'station_id')

# #Rename the columns in the filtered dataframe
# names(df) <- c(new_cols, 1:24)

# #Pivot the data into long format
# df <- df %>% pivot_longer(cols=str_pad(1:24, 1), 
# 	names_to='HourEnding', values_to='TrafficCount')

#Create a datetime column from date and time columns
# df$datetime <- with(df, ymd(date) + hms(paste0(HourEnding, ".00.00")))

#Aggregate the data by datetime and station_id
data <- aggregate(df$TrafficCount, by=list(df$datetime, df$station_id), FUN=sum)
names(data) <- c("datetime", 'station_id', 'TrafficCount')

data <- data[data$station_id == pdx, ] #Filter to chosen station id

data$Daily <- rep(1:24, 365) #Create a list of repeating hours of the day
data$Weekly <- wday(data$datetime, week_start=1) #Get the day of the week
data$Month <- month(data$datetime) #Get the month number
#Convert the datetime to a posix variable
data$datetime <- as.POSIXct(strptime(data$datetime, "%Y-%m-%d %H:%M:%S"), tz='GMT')

data$Time <- 1:dim(data)[1] #Column for the time index
data$station_name <- "Portland Metro" #Name the station
data$Date <- as.Date(data$datetime) #Get just the date portion of the datetime

#Get a set of lags of the Traffic Count variable to handle autocorrelation
data$Lag1 <- shift(data$TrafficCount, 1)
data$Lag2 <- shift(data$TrafficCount, 2)
data$Lag3 <- shift(data$TrafficCount, 3)
data$Lag24 <- shift(data$TrafficCount, 24)
data$Lag48 <- shift(data$TrafficCount, 48)
data$Lag168 <- shift(data$TrafficCount, 168)

################################################################################
# Exploratory Analysis of Data
################################################################################

#Generate a random starting index for a 2 week length for initial plot
len <- 24*7*2
start <- runif(1, 169, dim(data)[1])

ymin <- min(data$TrafficCount) #Get minimum of the data
ymax <- max(data$TrafficCount) #Used for defining the y-axis maximum

png("InitYearPlot.png") #Call the png function in base R
layout(matrix(1:2, nrow=2)) #Make a plotting grid
#Plot the full dataset of the traffic count per hour
plot(data$datetime, data$TrafficCount, type='l', col='black',
	xlab='Datetime', ylab='TrafficCount', ylim=c(0, ymax),
	main="Hourly Traffic Count, 2015",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
#Plot just the randomly chosen two week subset
plot(data[start:(len + start),]$datetime, 
	data[start:(len + start),]$TrafficCount, type='l', col='black',
	xlab='Datetime', ylab='TrafficCount', ylim=c(0, ymax),
	main="Hourly Traffic Count, 2 Weeks",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off() #Store the plot as a png file

#Get the periodogram dataframe, showing the spectral density at frequencies
ft <- periodogram(data$TrafficCount, plot=FALSE)
ft <- data.frame(ft$freq, ft$spec) #Put into a dataframe
names(ft) <- c("freq", 'spec') #Rename the columns
ft$Hours <- 1/ft$freq #The inverse of the frequency is the number of hours

png("InitSpec.png") #Develop a plot of the spectral density
layout(matrix(1:2, nrow=2, ncol=1)) #Two column plotting grid
#Plot the full frequency series of the spectral density
plot(ft$Hours, ft$spec, type='l', col='black', xlab='Hours', 
	ylab='Spectral Density', main="Spectral Density, All Frequencies",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
#Plot just the first week, 168 hours of the spectral density
plot(ft[ft$Hours <= 168,]$Hours, ft[ft$Hours <= 168,]$spec, type='l', 
	col='black', xlab='Hours', ylab='Spectral Density', 
	main="Spectral Density, 1 Week", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off() #Store as a png

################################################################################
# Generalized Additive Model Estimation
################################################################################

	############################################################################
	# First GAM model, testing chosen variables groupings
	############################################################################

#Specifying cyclic cubic splines with the correct periods, also constraining 
# them to have an additive relationship.
mod1 <- gam(TrafficCount ~ s(Daily, bs='cc', k=24) + s(Weekly, bs='cc', k=7) +
	s(Month, bs='cc', k=12),
	data = data, family = gaussian)

png("Mod1_2D.png", units='mm', res=300)
layout(matrix(1:3, nrow = 1)) #2 column plotting grid
plot(mod1, shade = TRUE, ylab='TrafficCount', cex.lab=1.5, cex.axis=1.5, 
	cex.main=4, main='Cyclic Splines') #Plot the model
dev.off()

png("Mod1_3D.png", units='mm', res=300)
vis.gam(mod1, view=c('Daily', 'Weekly'), n.grid = 50, theta = 135, phi = 32, zlab = "", ticktype = "detailed", color = "topo", 
	main="Model1 Daily, Weekly Additive",
	cex.lab=2.5, cex.axis=2.5, cex.main=4)
dev.off()

#Get the residuals, fitted values 
res_1 <- mod1$residuals
mod1_fit <- mod1$fitted.values
min1 <- min(ymin, min(mod1_fit)) #Get min and max for plotting
max1 <- max(ymax, max(mod1_fit))

png("Model1.png", units='mm', res=300) #Store the model
layout(matrix(1:4, nrow = 2, ncol=2)) #Get 2x2 grid for plotting
#Plot a random subset of the fitted values and the true values
plot(data[start:(start+len),]$datetime, 
	mod1_fit[start:(start+len)], ylab='TrafficCount', xlab='Datetime', 
	ylim=c(min1, max1), col='blue', type='l', main="Model 1 Fitted/True Values",
	cex.lab=1.5, cex.axis=1.5, cex.main=3)
lines(data[start:(start+len),]$datetime, data[start:(start+len),]$TrafficCount, ylab='TrafficCount', xlab='Datetime', ylim=c(min1, max1), col='red')
#Get the legend
legend('topright', c("Fitted", 'True'), fill=c("blue", 'red'))
#Plot the ACF up to 200 lags
plot(acf(res_1, lag=200, plot=FALSE), ylim=c(-0.4, 0.8), 
	main="ACF of Model 2 Residuals", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

#Plot the random subset for the true values
plot(data$datetime, res_1, 
	ylab='GAM Residuals', xlab='Datetime', ylim=c(min(res_1), max(res_1)), 
	col='green', type='l', main="Model 1 Residuals",
	cex.lab=1.5, cex.axis=1.5, cex.main=3)
legend('bottomright', 'Residuals', fill='green')

#Plot the PACF values for lags up to 200
plot(pacf(res_1, lag=200, plot=FALSE), ylim=c(-0.4, 0.8), 
	main="PACF of Model 2 Residuals",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

#See what AR and MA terms should be included based on the residuals
arima1 <- auto.arima(mod1$residuals, stationary=TRUE, seasonal=TRUE)

	############################################################################
	# Second GAM model relaxing additive model assumption, including AR terms
	############################################################################

#Add AR terms to model the autocorrelation present in the time series
mod2 <- gam(TrafficCount ~ t2(Daily, Weekly, bs=c('cc', 'cc'), k=c(24, 7)) +
			t2(Month, bs='cs') +
			t2(Lag1, bs='cs') + t2(Lag2, bs='cs') + t2(Lag3, bs='cs') +
			t2(Lag24, bs='cs') + t2(Lag168, bs='cs'),
			data = data, family=gaussian)

png("Mod2_3D.png", units='mm', res=300)
vis.gam(mod2, view=c('Daily', 'Weekly'), n.grid = 50, theta = 135, 
	phi = 32, zlab = "", ticktype = "detailed", color = "topo",
	main='Model2 Daily, Weekly with AR Terms',
	cex.lab=2.5, cex.axis=2.5, cex.main=4)
dev.off()

#Get the fitted values and residuals from the second model
mod2_fit <- mod2$fitted.values
res_2 <- mod2$residuals
min2 <- min(ymin, min(mod2_fit))
max2 <- max(ymax, max(mod2_fit))

png("Model2.png", units='mm', res=300) #Store the plot in a png file
layout(matrix(1:4, nrow = 2)) #Layout a 2x2 grid for plotting
#Plot the fitted and true values for the chosen subset of data
plot(data[start:(start+len),]$datetime, mod2_fit[start:(start+len)], 
	ylab='TrafficCount', xlab='Datetime', ylim=c(min2, max2), col='blue', 
	type='l', main="Model 2 Fitted/True Values",
	cex.lab=1.5, cex.axis=1.5, cex.main=3)
lines(data[start:(start+len),]$datetime, data[start:(start+len),]$TrafficCount, 
	ylab='TrafficCount', xlab='Datetime', ylim=c(min2, max2), col='red')
legend('topright', c("Fitted", 'True'), fill=c("blue", 'red'))
#Plot the ACF of the residuals of the 2nd model
plot(acf(res_2, lag=200, plot=FALSE), ylim=c(-0.4, 0.8),
	main="ACF of Model 2 Residuals", cex.lab=1.5, cex.axis=1.5, cex.main=5)
#Plot the residuals of the 2nd model for the full dataset
plot(data[169:nrow(data),]$datetime, res_2[169:nrow(data)], 
	ylab='GAM Residuals', xlab='Datetime', ylim=c(min(res_2), max(res_2)), 
	col='green', type='l', main="Model 2 Residuals",
	cex.lab=1.5, cex.axis=1.5, cex.main=3)
legend('bottomright', 'Residuals', fill='green')
#Plot the PACF of the residuals of the 2nd model up to 200 lags
plot(pacf(res_2, lag=200, plot=FALSE), ylim=c(-0.4, 0.8), 
	main="PACF of Model 2 Residuals", cex.lab=1.5, cex.axis=1.5, cex.main=5)
dev.off()

auto.arima(mod2$residuals, stationary=TRUE, seasonal=TRUE)
