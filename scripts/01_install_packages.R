# First we are creating an object with a list of the packages that we'll need

list.of.packages <- c('readxl', 'ncdf4', 'cmsaf', 'data.table', 
                      'ggplot2', 'magrittr', 'plyr', 'dplyr', 'utils', 'naniar',
                      'lubridate', 'raster', 'tidyverse', 'imputeTS', 'readr', 'ncdf4',
                      'rWind', 'broman', 'car', 'ggplot2', 'tidyr', 'extRemes', 'grDevices',
                      'splitstackshape', 'TideHarmonics', 'viridis', 'DT', 'webshot', 'htmltools', 'knitr') 

# Now we will check to see if any of the packages required are not yet on our system

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

# Any package missing will be added to the ‘new.packages’ object 
# which can then be used to install any missing ones

if(length(new.packages)) install.packages(new.packages)

