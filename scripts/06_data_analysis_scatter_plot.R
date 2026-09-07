
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 06_data_analysis_scatter_plot.R

# Example scatter plot - wave data
###########################################################################

#load packages
library(tidyverse)
library(magrittr)

# import data
wave <- read.csv('data/modelled/wave_compiled.csv')
colnames(wave)

# Bin wave height into categories
wave <- wave %>%
  mutate(height_bin = cut(hs, breaks = seq(0, 16, by = 2), include.lowest = TRUE))

# Scatter Diagram for Incident wave conditions and Kernel Density
#color palette
cols <- rev(rainbow(20)[-c(3, 6, 8,  10, 12, 13, 14, 16, 17, 18, 19, 20)])

plot_kernel <- ggplot(wave, aes(x=tp, y=hs)) + xlab("Peak wave period (Tp) (s)") +
  ylab("Significant wave height (Hs) (m)") +
  geom_bin2d(bins = 110, aes(fill = after_stat(count))) +
  labs(fill = "Number of Occurrences") +
  scale_fill_gradientn(colours=cols, breaks = c(1,500, 1000, 1500, 2000, 2500))+  # Numbers in this line should be changed for each project. 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,30)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 16)) + 
  ggtitle("Tp (s) vs Hs (m)") +
  theme_linedraw()+
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 6))
plot_kernel # print the plot so you can visualise it & check for errors

ggsave(plot_kernel, file="output/hs_tp_kernel_density.jpeg", dpi = 600, height = 5, width = 7, units = c("in"))

