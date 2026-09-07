
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 05_data_analysis_tabular_stats.R

# Example stats - wind data
###########################################################################

# load packages
library(tidyverse)
library(magrittr)

# import dataset
# import wind file
wind <- read.csv('data/modelled/era5_wind_10m.csv')

# check column names
colnames(wind)

#-----------------------------------------------------
# Function to parse mixed formats
parse_mixed_datetime <- function(x) {
  # Try full timestamp with dash
  dt <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Try full timestamp with colon
  na_idx <- is.na(dt)
  dt[na_idx] <- as.POSIXct(x[na_idx], format = "%Y:%m:%d %H:%M:%S", tz = "UTC")
  
  # Try date-only with dash
  na_idx <- is.na(dt)
  dt[na_idx] <- as.POSIXct(x[na_idx], format = "%Y-%m-%d", tz = "UTC")
  
  # Try date-only with colon
  na_idx <- is.na(dt)
  dt[na_idx] <- as.POSIXct(x[na_idx], format = "%Y:%m:%d", tz = "UTC")
  
  return(dt)
}

# Apply to your dataframe
wind$datetime <- parse_mixed_datetime(wind$datetime)

# Check
any(is.na(wind$datetime))
#---------------------------------------------------------------------------------

#TABULAR STATS
# Create a new column 'months' using the 'timsetamp column as input. 
# This will be needed when we want to generate monthly statistics
wind$month <- factor(format(wind$datetime, "%m"), labels = month.name)

## monthly statistics - wind_speed
monthly_stats <- wind %>%
  group_by(month) %>%
  summarise("mean" = mean(wind_speed),
            "median" = median(wind_speed),
            "sd" = sd(wind_speed),
            "max" = max(wind_speed),
            "min" = min(wind_speed))

## overall statistics - wind_speed
overall_stats <- wind %>%
  summarise("mean" = mean(wind_speed),
            "median" = median(wind_speed),
            "sd" = sd(wind_speed),
            "max" = max(wind_speed),
            "min" = min(wind_speed))
overall_stats$month <- c('Overall') #add a new column to the dataframe
overall_stats <- overall_stats[, c(6, 1:5)] #re-arrange the columns so that 'overall' is first. 

table_wind_speed <- rbind(monthly_stats, overall_stats) #combine the two dataframes
write.csv(table_wind_speed, 'output/wind_speed_stats.csv') #save your output

