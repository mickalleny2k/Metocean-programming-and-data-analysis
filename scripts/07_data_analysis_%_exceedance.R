
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 07_data_analysis_%_exceedance.R

# % Exceedance table
# This is an example piece of code for wind speed.
# For other datasets e.g., wave heights and current speed, change thresholds accordingly
###########################################################################

# import wind file
wind <- read.csv('data/modelled/era5_wind_10m.csv')

# check column names
colnames(wind)

# Define speed thresholds for exceedance calculation
thresholds <- seq(0, 20, by = 2) # defining thresholds from 0 to 20 m/s at 2 m/s intervals

# Calculate % exceedance for each threshold
exceedance <- sapply(thresholds, function(thresh) {
  mean(wind$wind_speed > thresh) * 100
})

# Create a data frame for the table
exceedance_table <- data.frame(
  Threshold = thresholds,
  Percent_Exceedance = exceedance
)

print(exceedance_table)

# Export as CSV file
write.csv(exceedance_table, 'output/exceedance_wind_speed_table.csv')

