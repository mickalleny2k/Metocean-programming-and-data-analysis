
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 04_data_analysis_rose_plot.R

# Example rose plot - wind data
###########################################################################

# load packages
library(openair)

# import wind file
wind <- read.csv('data/modelled/era5_wind_10m.csv')

# check column names
colnames(wind)

# basic wind rose plot
png("output/windrose_plot.png", width = 1500, height = 1500, res = 200)
windRose(wind, ws = "wind_speed", wd = "wind_dir_from", angle = 22.5, breaks= 6, ws.int = 4, calm.thresh = 0, dig.lab= 1, 
         width = 1.5, key.footer = "wind speed (m/s)", annotate = FALSE, main = 'Wind rose (10 m above sea level)', 
         ylab = 'Frequency of counts by wind direction (%)')
dev.off()


