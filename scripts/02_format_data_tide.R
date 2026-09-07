
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 02_format_data_tide.R

# Formats tide data for harmonic analysis
###########################################################################

library(magrittr)
library(dplyr)
library(lubridate)

#Format tide file
#import file
tide <- read.csv('data/modelled/water_level_currents.csv')

tide$time <- as.POSIXct(paste(tide$Time),
                        tz = "GMT",
                        format = "%Y-%m-%d %H:%M:%S")

#wind - arrange time according to earliest to oldest
tide <- tide %>% 
  arrange(time)
tide <- tide[, -c(1)]
tide <- tide[, c(4, 1:3)]
colnames(tide) <- c('time', 'c_dir', 'c_speed', 'water_level')
write.csv(tide, "output/tide_compiled.csv", row.names = FALSE)
