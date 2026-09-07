
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 03_data_analysis_tidal_levels.R

# Harmonic analysis on water levels
###########################################################################

library(TideHarmonics)

#1. Use the function ftide() to conduct harmonic analysis on your water level time series
tide.der <- ftide(tide$water_level, tide$time, hc114) # the data frame 'tide' is already in your environment from script 02.

# 2. From your tidal fit, extract tidal datums (from MLWS to MHWS) which are computed based on 
# properties of tidal constituents at your site and stored in 'features1'
SD_tide <- tide.der$features1
SD_tide.df <- as.data.frame(SD_tide)
SD_tide.df <- as.data.frame(t(SD_tide.df))

# 3. Need to calculate HAT and LAT. This is not given in the tide model. 
# Using the tidal constituent properties calculated above, predict a 
# long-term time series over a 18.6 year nodal period

# Define a fine grid of times over the nodal period (e.g. 18.6 years)
start <- min(tide$time)
end   <- start + lubridate::years(19)  # approximate
 

# Use the function predict() from the TideHarmonics R package to predict a 
# long term time series of astronomical water levels from the harmonic analysis results above (stored in tide.der)
predicted <- predict(tide.der, start, end, 1)

# Extract maxima and minima of this record to obtain HAT and LAT
HAT <- max(predicted, na.rm = TRUE)
LAT <- min(predicted, na.rm = TRUE)

#Combine all tidal datums into one dataframe and label appropriately
datums <- cbind(SD_tide.df, HAT, LAT)

#re-arrange in order of high water to low water
datums <- datums[, c('HAT', 'MHWS', 'MHWN', 'MSL', 'MLWN', 'MLWS', 'LAT')]


#Transform dataframe
datums_t <- as.data.frame(t(datums))

# Correct values so that MSL = 0
datums_t$SD_tide <- datums_t$SD_tide - datums_t["MSL", "SD_tide"]

write.csv(datums_t, 'output/tidal_datums.csv')

