###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 11_%occurrence_joint_distribution.R

# Expected dataset:
#  A dataframe called 'data' with columns 'vhub', 'dw', 'hs', 'tp' 
#  where vhub is wind speed at hub height and dw is wind direction

# An example dataset is provided for this
###########################################################################

library(data.table)
library(magrittr)
library(dplyr)

# Import data
data <- read.csv('data/modelled/joint_distribution_example_data.csv')

count <- nrow(data)


# 1. Wind speed bins: 2 m/s
data$category_vhub <- cut(
  data$vhub,
  breaks = seq(0, 49, by = 2)  # start = 0, end = 49, step = 2 m/S
)

# 2. Significant wave height (hs) bins: 0.5 m (centered at 0.25 increments)
data$category_hs <- cut(
  data$hs,
  breaks = seq(0, 15, by = 0.5)  # start=0, end=15, step=0.5
)

# 3. Peak period (tp) bins: 0.5 s (up to max value in dataset)
tp_max <- ceiling(max(data$tp, na.rm = TRUE))  # find max tp and round up
data$category_tp <- cut(
  data$tp,
  breaks = seq(0, tp_max, by = 0.5)
)



#----------------------------------------------------------------------------------
# 1. Generate scatter table per wind speed category with totals and cumulative sums
#----------------------------------------------------------------------------------

cont <- unique(data$category_vhub)
df <- list()  # store all categories

for(i in seq_along(cont)) {
  
  # Subset data for current wind speed category
  temp <- data[data$category_vhub == cont[i], ]
  
  # Contingency table of hs vs tp frequencies
  t <- table(temp$category_hs, temp$category_tp)
  
  # Convert to percentages of total dataset
  t_percent <- (t / count) * 100
  t_percent_matrix <- as.data.frame.matrix(t_percent)
  
  # -------------------------------
  # Row sums and cumulative sum
  t_percent_matrix$sum <- rowSums(t_percent_matrix)
  t_percent_matrix$cusum <- cumsum(t_percent_matrix$sum)
  
  # -------------------------------
  # Column sums and cumulative sum (exclude 'sum' and 'cusum' columns)
  col_sums <- colSums(t_percent_matrix[, 1:(ncol(t_percent_matrix)-2)])
  col_cusum <- cumsum(col_sums)
  
  # Add totals as extra rows
  totals_row <- c(col_sums, sum(col_sums), NA)  # NA for 'cusum' column
  cum_row <- c(col_cusum, NA, NA)               # NA placeholders for 'sum' and 'cusum'
  
  t_percent_matrix <- rbind(
    t_percent_matrix,
    ColTotals = totals_row,
    ColCumSum = cum_row
  )
  
  # Store in list using category name
  name <- paste0(cont[i])
  df[[name]] <- t_percent_matrix
}

# Export as separate CSVs if needed
for(name in names(df)){
  write.csv(df[[name]], paste0('output/joint_prob_bins_', gsub("[^0-9a-zA-Z]", "_", name), '.csv'), row.names = TRUE)
}


#----------------------------------------------------------------------------------
# 2. Generate scatter table for ALL wind speeds with totals and cumulative sums
#----------------------------------------------------------------------------------

n <- table(data$category_hs, data$category_tp)       # contingency table
n_percent <- (n / count) * 100                      # convert counts to percent
n_percent_matrix <- as.data.frame.matrix(n_percent) # convert to data.frame

# --- Row sums and cumulative sum of rows ---
n_percent_matrix$sum <- rowSums(n_percent_matrix)   # sum across each row
n_percent_matrix$cusum <- cumsum(n_percent_matrix$sum) # cumulative sum down rows
# Column sums and cumulative sums (exclude 'sum' and 'cusum' columns)
col_sums <- colSums(n_percent_matrix[, 1:(ncol(n_percent_matrix)-2)])  # only original bins
col_cusum <- cumsum(col_sums)

# Construct row for totals
totals_row <- c(col_sums, sum(col_sums), NA)  # sum(col_sums) for row 'sum', NA for 'cusum' column
cum_row <- c(col_cusum, NA, NA)               # NA placeholders for 'sum' and 'cusum' columns

# rbind with matching column names
n_percent_matrix_totals <- rbind(
  n_percent_matrix,
  ColTotals = totals_row,
  ColCumSum = cum_row
)

# Export CSV
write.csv(n_percent_matrix_totals, 'output/joint_prob_All_wind_speeds.csv', row.names = TRUE)
