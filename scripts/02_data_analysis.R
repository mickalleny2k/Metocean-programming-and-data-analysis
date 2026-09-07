#load packages
library(tidyverse)
library(magrittr)
library(readr)
library(lubridate)

#Import wind timeseries
hs <- read.csv('data/era5_2015_2020_hs.csv')
tp <- read.csv('data/era5_2015_2020_tp.csv')

#format files so they are ready to combine.

# The 'Time' column in the imported CSV file is not yet recognised as time in R. 
# Therefore, we have to create a new 'time' column using the function 'as.POISVct()' with 'time' as input. 
# As our new column is the same name as our existing column it is overwriting it.
hs$time <- as.POSIXct(paste(hs$time),
                      tz = "GMT",
                      format = "%d/%m/%Y %H:%M")

#hs - arrange time according to earliest to oldest
hs <- hs %>% 
  arrange(time)

#Same for mwd & tp but also remove time column as this is retained in hs dataset.
tp$time <- as.POSIXct(paste(tp$time),
                      tz = "GMT",
                      format = "%d/%m/%Y %H:%M")
tp <- tp %>% 
  arrange(time)

#create a new object 'wave' and combine all columns using function 'cbind()'
wave <- cbind(hs, tp)

#Open wave dataframe you just created in your Environment in the top right panel to see what columns you want to remove.
#remove all unwanted columns
wave <- wave[, -c(1, 4:7, 9:10)]

#save formatted dataset
write.csv(wave, 'output/dataset/wave_compiled.csv')


# Bin wave height into categories
wave <- wave %>%
  mutate(height_bin = cut(hs, breaks = seq(0, 14, by = 2), include.lowest = TRUE))

# Scatter Diagram for Incident wave conditions and Kernel Density
#color palette
cols <- rev(rainbow(20)[-c(3, 6, 8,  10, 12, 13, 14, 16, 17, 18, 19, 20)])

plot_kernel <- ggplot(wave, aes(x=tp, y=hs)) + xlab("Peak wave period (s)") +
  ylab("Significant wave height (m)") +
  geom_bin2d(bins = 160, aes(fill = after_stat(count))) +
  labs(fill = "Number of Occurrences") +
  scale_fill_gradientn(colours=cols, breaks = c(1,20, 40, 60, 80, 100))+  # Numbers in this line should be changed for each project. 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,22.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 16)) + 
  ggtitle("Wave period vs height") +
  theme_linedraw()+
  theme(plot.title = element_text(hjust = 0.1), legend.text = element_text(size = 10))+
guides(fill = guide_colorbar(barwidth = 2, barheight = 12))
plot_kernel # print the plot so you can visualise it & check for errors

ggsave(plot_kernel, file="output/plots/hs_tp_kernel_density.jpeg", dpi = 600, height = 3, width = 4.5, units = c("in"))


#TABULAR STATS
# Create a new column 'months' using the 'timsetamp column as input. 
# This will be needed when we want to generate monthly statistics
wave$month <- factor(format(wave$time, "%m"), labels = month.name)

## monthly statistics - Hs
monthly_stats <- wave %>%
  group_by(month) %>%
  summarise("mean" = mean(hs),
            "median" = median(hs),
            "sd" = sd(hs),
            "max" = max(hs),
            "min" = min(hs))

## overall statistics - Hs
overall_stats <- wave %>%
  summarise("mean" = mean(hs),
            "median" = median(hs),
            "sd" = sd(hs),
            "max" = max(hs),
            "min" = min(hs))
overall_stats$month <- c('Overall') #add a new column to the dataframe
overall_stats <- overall_stats[, c(6, 1:5)] #re-arrange the columns so that 'overall' is first. 

table_hs <- rbind(monthly_stats, overall_stats) #combine the two dataframes
write.csv(table_hs, 'output/tables/hs_stats.csv') #save your output


