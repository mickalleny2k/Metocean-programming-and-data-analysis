
###########################################################################
# GG6581
# Week 9
# Author: Dr. Shauna Creane
# 10_extreme_value_analysis.R

# Example of how to carry out extreme value analysis.
# The following sections are included in this script:

# 1. Simulate hourly wind speeds - synthetic dataset
# 2. Automatic GP Threshold Selection with Mean Residual Life (MRL) Plot
# 3. Declustering for POT (24-hour storm) - create function
# 4. Apply declustering
# 5. Fit distributions, including GEV & Gumbel to AM;GP to declustered peaks (POT)
# 6. Compute Goodness-of-fit statistics (KS & PCC)
# 7. Optional Additional Goodness of fit plots
# 8. Return levels
# 9. Compute parametric bootstrap CI for best-fit distribution (either Gumbel, GEV, or GP)
# 10. Plot return levels for the best-fit distribution with computed confidence intervals
###########################################################################


library(evd)
library(tidyverse)
library(fitdistrplus)
library(extRemes)
library(boot)
library(ggplot2)

# -------------------------------
# 1. Simulate hourly wind speeds - synthetic dataset
# -------------------------------

set.seed(123)

# -------------------------------
# A. Set parameters
# -------------------------------
years <- 1980:2010
days_per_year <- 365
hours_per_day <- 24

# Total hours
n_hours <- length(years) * days_per_year * hours_per_day

# -------------------------------
# B. Gamma parameters for wind speed
# -------------------------------
# Target: mean ≈ 9 m/s, max ≈ 30 m/s
# Gamma mean = shape * scale
# Gamma variance = shape * scale^2

target_mean <- 5
target_max <- 30

# Rough estimate: shape & scale
# Shape controls skewness (larger shape → less skewed)
shape <- 2
scale <- target_mean / shape  # scale = mean / shape

# Simulate hourly wind speeds
wind_data <- data.frame(
  datetime = seq.POSIXt(
    from = as.POSIXct("1980-01-01 00:00"),
    by = "hour",
    length.out = n_hours
  ),
  wind_speed = rgamma(n_hours, shape = shape, scale = scale)
)

# Plot histogram of wind speeds
hist(
  wind_data$wind_speed,
  breaks = 50,
  col = "skyblue",
  main = "Simulated Hourly Wind Speed Distribution",
  xlab = "Wind speed (m/s)",
  ylab = "Frequency",
  xlim = c(0, target_max)
)

# -------------------------------
# C. Check extremes
# -------------------------------
summary(wind_data$wind_speed)
max(wind_data$wind_speed)

# -------------------------------
# E. Compute annual maxima
# -------------------------------
annual_max <- wind_data %>%
  mutate(year = format(datetime, "%Y")) %>%
  group_by(year) %>%
  summarise(max_wind = max(wind_speed)) %>%
  pull(max_wind)

# Inspect
head(annual_max)
summary(annual_max)

# annual max
x <- annual_max

# full dataset
y <- wind_data$wind_speed



# -------------------------------
# 2. Automatic GP Threshold Selection with Mean Residual Life (MRL) Plot
# -------------------------------
select_gp_threshold_plot <- function(data, min_exc = 30) {
  
  # A. Candidate thresholds
  # Generate a sequence of high quantiles as candidate thresholds
  probs <- seq(0.99, 0.9999, length.out = 50) # 99% to 99.99% quantiles
  thresholds <- quantile(data, probs)         # actual data values for candidate thresholds
  
  
  # B. Compute mean excess for each threshold
  mean_excess <- sapply(thresholds, function(u) {
    exc <- data[data > u] - u                     # excesses above threshold u
    if(length(exc) < min_exc) return(NA)          # skip if too few exceedances
    mean(exc)                                     # mean residual life (mean excess)
  })
  
  
  # C. Remove thresholds where mean excess could not be computed
  valid <- !is.na(mean_excess)
  thresholds <- thresholds[valid]
  mean_excess <- mean_excess[valid]
  probs <- probs[valid]
  
  
  # D.Stable threshold
  # Compute relative differences between consecutive mean excesses
  diffs <- diff(mean_excess)
  stable_idx <- which(diffs / mean_excess[-length(mean_excess)] < 0.0001) # near-flat region
  
  if(length(stable_idx) == 0) {
    # No stable threshold found → fallback
    warning("No stable threshold found, using 99th percentile.")
    thresh <- quantile(data, 0.99)
    thresh_prob <- 0.99
  } else {
    # Use the first "stable" threshold
    thresh <- thresholds[stable_idx[1]]
    thresh_prob <- probs[stable_idx[1]]
  }
  
  
  # E. Count exceedances above selected threshold
  n_exceedances <- sum(data > thresh)
  
  
  # F. Plot Mean Residual Life (MRL)
  mrl_df <- data.frame(Threshold = thresholds, MeanExcess = mean_excess)
  p <- ggplot(mrl_df, aes(x = Threshold, y = MeanExcess)) +
    geom_line(color = "blue") +         # MRL curve
    geom_point() +                      # points for each threshold
    geom_vline(xintercept = thresh, color = "red", linetype = "dashed", size = 1) + # selected threshold
    annotate("text", x = thresh, y = max(mean_excess, na.rm = TRUE), 
             label = paste0("Threshold = ", round(thresh, 2),
                            "\nPercentile = ", round(thresh_prob*100, 2), "%",
                            "\nExceedances = ", n_exceedances),
             hjust = -0.1, vjust = 1.2, color = "red") +       # annotate threshold info
    labs(title = "Mean Residual Life (MRL) Plot",
         x = "Threshold (u)",
         y = "Mean Excess above u") +
    theme_minimal()
  
  print(p)  # explicitly print the plot inside the function
  
  
  
  # G. Return all relevant info
  return(list(
    threshold = thresh, # selected threshold value
    percentile = thresh_prob, # percentile of selected threshold
    n_exceedances = n_exceedances, # number of exceedances above threshold
    plot = p # ggplot object of MRL
  ))
}


# -------------------------------
# H. Apply function

# Call the automatic GP threshold selection function on dataset 'y'
gp_info <- select_gp_threshold_plot(y)

# Print the selected threshold value
cat("Selected GP threshold:", gp_info$threshold, "\n")

# Print the percentile of the selected threshold (e.g., 99th percentile)
cat("Corresponding percentile:", gp_info$percentile, "\n")

# Print the number of exceedances above the selected threshold
cat("Number of exceedances above threshold:", gp_info$n_exceedances, "\n")



# -------------------------------
# 3. Declustering for POT (Peaks Over Threshold) using a 48-hour storm window
# -------------------------------
decluster_exceedances <- function(data, threshold, min_gap_hours = 48) {
  
  # Parameters:
  # data: a dataframe with at least 'datetime' and 'wind_speed' columns
  # threshold: the GP threshold used for identifying exceedances
  # min_gap_hours: minimum time gap (in hours) to separate storm events
  #                assuming data is hourly; adjust if using different time resolution
  
  # A. Select only exceedances above the threshold
  exceed <- data %>% filter(wind_speed > threshold)
  
  # B. If there are no exceedances, return an empty dataframe
  if(nrow(exceed) == 0) return(data.frame())  # no exceedances
  
  # C. Order exceedances chronologically
  exceed <- exceed %>% arrange(datetime)
  
  # D. Initialize cluster assignment
  cluster_id <- numeric(nrow(exceed)) # vector to store cluster numbers
  cluster_counter <- 1                # start first cluster
  cluster_id[1] <- cluster_counter    # assign first exceedance to first cluster
  
  # E. Loop through exceedances to assign clusters based on time gap
  for(i in 2:nrow(exceed)) {
    # Calculate time difference from previous exceedance
    diff_hours <- as.numeric(difftime(exceed$datetime[i], exceed$datetime[i-1], units = "hours"))
    
    # If gap exceeds min_gap_hours, start a new cluster
    if(diff_hours > min_gap_hours) {
      cluster_counter <- cluster_counter + 1
    }
    # Assign cluster number
    cluster_id[i] <- cluster_counter
  }
  
  exceed$cluster <- cluster_id # attach cluster IDs to exceedances
  
  
  # F. Keep only the maximum wind_speed from each cluster
  declustered <- exceed %>%
    group_by(cluster) %>%
    summarise(
      datetime = datetime[which.max(wind_speed)],   # datetime of max wind_speed in cluster
              wind_speed = max(wind_speed)          # maximum wind_speed in cluster
      ) %>%
    ungroup()
  
  # G. Return the declustered exceedances
  return(declustered)
}


# -------------------------------
# 4. Apply declustering
# -------------------------------

# Set threshold for GP exceedances using the MRL-selected threshold
threshold <- gp_info$threshold  

# Apply declustering function to wind data
#    - Groups exceedances separated by at least 48 hours into clusters
#    - Keeps only the maximum wind speed per cluster
declustered_exceedances <- decluster_exceedances(wind_data, threshold, min_gap_hours = 48)

# Print summary of exceedances before and after declustering
#    - Original: all points above threshold
cat("Original exceedances:", sum(y > threshold), "\n")

#    - Declustered: only cluster maxima, reduces correlated extremes
cat("Declustered exceedances (48h storm):", nrow(declustered_exceedances), "\n")


# -------------------------------
# 5. Fit distributions, including GEV & Gumbel to AM;GP to declustered peaks (POT)
# -------------------------------

# Note:
# - Fits are performed once per distribution using Maximum Likelihood Estimation (MLE) by default
# - For L-moments based fitting, you could also use the 'lmom' library
# - 'x' is typically the vector of annual maxima (for GEV/Gumbel)
# - 'declustered_exceedances$wind_speed' is used for GP fitting (POT approach)

fit_list <- list(
  # A. Generalized Extreme Value (GEV) for annual maxima
  GEV = fevd(x, type = "GEV"),
  
  # B. Gumbel distribution for annual maxima (special case of GEV with shape=0)
  Gumbel = fevd(x, type = "Gumbel"),
  
  # C. Generalized Pareto (GP) distribution for peaks-over-threshold
  #    Uses declustered exceedances and threshold from MRL
  GP = fevd(declustered_exceedances$wind_speed, type = "GP", threshold = threshold)
)

# Store additional GP info for reporting
fit_gp <- fit_list$GP  # extract the GP fit object
fit_gp$threshold_percentile <- gp_info$percentile  # save percentile used for threshold
fit_gp$n_exceedances <- gp_info$n_exceedances      # save number of exceedances above threshold



# -------------------------------
# 6. Goodness-of-fit statistics
# -------------------------------

# Computes CDF values for goodness-of-fit tests
calc_cdf <- function(obs, fit, dist_name){
  switch(dist_name,
         # A. GEV CDF using fitted location, scale, shape
         GEV = pgev(obs, loc=fit$results$par["location"], scale=fit$results$par["scale"], shape=fit$results$par["shape"]),
         
         # B. Gumbel CDF (special case of GEV with shape=0)
         Gumbel = pgumbel(obs, loc=fit$results$par["location"], scale=fit$results$par["scale"]),
         
         # C. GP CDF for Peaks-Over-Threshold
         GP = {
           thresh <- fit$threshold; sigma <- fit$results$par["scale"]; xi <- fit$results$par["shape"]
           
           # For observations below threshold, CDF=0
           if(length(obs)==0) return(numeric(0))
           
           # Standard GP CDF formula
           p <- ifelse(obs <= thresh, 0, 1 - (1 + xi*(obs - thresh)/sigma)^(-1/xi))
           p
         })
}

# ---------------------------------
# 6a. Kolmogorov-Smirnov statistic
ks_stat <- function(obs, fit, dist_name){
  cdf_vals <- calc_cdf(obs, fit, dist_name)    # compute fitted CDF values
  cdf_vals <- cdf_vals + runif(length(cdf_vals), -1e-8, 1e-8)  # jitter to avoid ties
  ks.test(cdf_vals, "punif")$statistic   # KS statistic vs uniform
}

# ------------------------------------
# 6b. Probability Plot Correlation Coefficient (PPCC)
ppcc_stat <- function(obs, fit, dist_name){
  n <- length(obs); 
  obs_sorted <- sort(obs)   # sort observations
  probs <- (1:n)/(n+1)      # plotting positions
  q_fit <- switch(dist_name,
                  
                  # GEV quantiles
                   GEV = qgev(probs, loc=fit$results$par["location"], scale=fit$results$par["scale"], shape=fit$results$par["shape"]),
                  # Gumbel quantiles
                  Gumbel = qgumbel(probs, loc=fit$results$par["location"], scale=fit$results$par["scale"]),
                  # GP quantiles
                  GP = {
                    thresh <- fit$threshold; sigma <- fit$results$par["scale"]; xi <- fit$results$par["shape"]
                    thresh + sigma * ((1 - probs)^(-xi) - 1)/xi
                  })
  cor(obs_sorted, q_fit)   # correlation between empirical and fitted quantiles
}

#-----------------------------------
# 6c. Compute GOF statistics for all models

gof_table <- data.frame(
  Distribution = names(fit_list),
  
  # KS statistic for each distribution
  KS_Statistic = c(
    ks_stat(annual_max, fit_list$GEV,"GEV"),
    ks_stat(annual_max, fit_list$Gumbel,"Gumbel"),
    ks_stat(declustered_exceedances$wind_speed, fit_list$GP,"GP")
  ),
  
  # PPCC statistic for each distribution
  PPCC = c(
    ppcc_stat(annual_max, fit_list$GEV,"GEV"),
    ppcc_stat(annual_max, fit_list$Gumbel,"Gumbel"),
    ppcc_stat(declustered_exceedances$wind_speed, fit_list$GP,"GP")
  )
)
print(gof_table)  # inspect table of GOF statistics


# ----------------------------------------------------------
# 6d. Extract AIC/BIC (all models)

get_info_criteria <- function(fit) {
  loglik <- fit$results$value        # log-likelihood from fit
  k <- length(fit$results$par)       # number of fitted parameters
  
  # Akaike Information Criterion
  aic <- -2 * loglik + 2 * k
  
  # Bayesian Information Criterion
  bic <- -2 * loglik + log(length(fit$x)) * k
  
  return(c(AIC = aic, BIC = bic))
}

# Build AIC/BIC matrix for all fits
info_mat <- t(sapply(fit_list, get_info_criteria))


# ------------------------------------------------------
# 6e. Combine GOF stats and info criteria

comparison_table <- data.frame(
  Distribution = names(fit_list),
  KS = gof_table$KS_Statistic,
  PPCC = gof_table$PPCC,
  AIC = info_mat[, "AIC"],
  BIC = info_mat[, "BIC"]
)

# ------------------------------------------------------
# 6f. QQ MSE

qq_error <- function(data, fit, dist_name, threshold = NULL) {
  # n = number of observations
  n <- length(data)
  
  # Compute plotting positions for empirical quantiles
  probs <- (1:n)/(n+1)
  
  # Sort the data to get empirical quantiles
  obs <- sort(data)
  
  # Compute theoretical quantiles from the fitted distribution
  theo <- switch(dist_name,
                 # GEV quantiles
                 GEV = qgev(probs, loc=fit$results$par["location"], scale=fit$results$par["scale"], shape=fit$results$par["shape"]),
                 
                 # Gumbel quantiles (special case of GEV with shape=0)
                 Gumbel = qgumbel(probs, loc=fit$results$par["location"], scale=fit$results$par["scale"]),
                 
                 # GP quantiles for Peaks-Over-Threshold
                 GP = {
                   sigma <- fit$results$par["scale"]
                   xi <- fit$results$par["shape"]
                   u <- threshold
                   
                   # Standard GP quantile formula
                    u + sigma * ((1 - probs)^(-xi) - 1)/xi
                 })
  
  # Compute mean squared error between empirical and theoretical quantiles
  mean((obs - theo)^2)  # MSE as a goodness-of-fit measure
}

# Apply QQ MSE calculation for all fitted distributions
comparison_table$QQ_MSE <- c(
  qq_error(annual_max, fit_list$GEV, "GEV"),
  qq_error(annual_max, fit_list$Gumbel, "Gumbel"),
  qq_error(declustered_exceedances$wind_speed, fit_list$GP, "GP", threshold)
)

# -------------------------------
# Ranking system: ONLY apply to GEV and Gumbel
# -------------------------------
rank_models <- function(df) {
  
  # Keep only GEV and Gumbel rows; GP is excluded because POT/GP behaves differently
  df_ranked <- df[df$Distribution %in% c("GEV","Gumbel"), ]  
  
  
  # Rank each goodness-of-fit metric
  df_ranked$KS_rank  <- rank(df_ranked$KS)  # Lower KS statistic → better fit, so rank normally
  df_ranked$AIC_rank <- rank(df_ranked$AIC) # Lower AIC → better fit, rank normally
  df_ranked$BIC_rank <- rank(df_ranked$BIC)  # Lower BIC → better fit, rank normally
  df_ranked$PPCC_rank <- rank(-df_ranked$PPCC) # Higher PPCC → better fit, so rank in descending order (negate values)
  df_ranked$QQ_rank <- rank(df_ranked$QQ_MSE)  # Lower QQ MSE → better fit, rank normally
  
 # Combine all ranks into a total score
  df_ranked$TotalScore <- df_ranked$KS_rank +
    df_ranked$AIC_rank +
    df_ranked$BIC_rank +
    df_ranked$PPCC_rank +
    df_ranked$QQ_rank
  
  # Sort models by total score (lower is better)
  df_ranked <- df_ranked[order(df_ranked$TotalScore), ]
  return(df_ranked)
}

# Apply ranking function to comparison table
ranked_models <- rank_models(comparison_table)
print(ranked_models)

# Automatically select best distribution
# Take the top-ranked distribution among GEV and Gumbel
best_fit_auto <- ranked_models$Distribution[1]
cat("Best distribution (statistical, GEV vs Gumbel only):", best_fit_auto, "\n")

# -------------------------------
# Function to create Q-Q plot for any fitted distribution (visual inspection)
# -------------------------------
plot_qq <- function(data, fit, dist_name, threshold = NULL) {
  n <- length(data) # Number of observations
  probs <- (1:n) / (n + 1) # Empirical plotting positions
  data_sorted <- sort(data) # Sort observed data in ascending order
  
  # Compute theoretical quantiles from fitted distribution
  q_theoretical <- switch(dist_name,
                          
                          # For GEV: use fitted location, scale, shape
                          GEV = qgev(probs, loc = fit$results$par["location"], scale = fit$results$par["scale"], shape = fit$results$par["shape"]),
                          
                          # For Gumbel: shape = 0, only location and scale
                          Gumbel = qgumbel(probs, loc = fit$results$par["location"], scale = fit$results$par["scale"]),
                          
                          # For GP: use threshold-based calculation
                          GP = {
                            sigma <- fit$results$par["scale"]
                            xi <- fit$results$par["shape"]
                            u <- threshold
                            u + sigma * ((1 - probs)^(-xi) - 1)/xi
                          })
  
  # Combine theoretical and observed values into a data frame
  df <- data.frame(Theoretical = q_theoretical, Observed = data_sorted)
  
  # Generate Q-Q plot using ggplot2
  ggplot(df, aes(x = Theoretical, y = Observed)) + # Observed vs theoretical points
    geom_point(color = "darkblue") +               # 45-degree reference line
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Q-Q Plot:", dist_name),  # Add title dynamically
         x = "Theoretical Quantiles",
         y = "Observed Quantiles") +
    theme_minimal()                                # Clean minimal theme
}

# Plot QQ for visual inspection
plot_qq(annual_max, fit_list$GEV, "GEV")
plot_qq(annual_max, fit_list$Gumbel, "Gumbel")
plot_qq(declustered_exceedances$wind_speed, fit_list$GP, "GP", threshold)


# -----------------------------------------------------
# 7. Optional Additional Goodness of fit plots
# ------------------------------------------------
# A. Histograms

# Create a generic function that can be applied to all models
plot_hist_with_fit <- function(data, fit, dist_name, threshold = NULL) {
  
  x_grid <- seq(min(data), max(data), length.out = 500)
  
  density_vals <- switch(dist_name,
                         GEV = dgev(x_grid,
                                    loc = fit$results$par["location"],
                                    scale = fit$results$par["scale"],
                                    shape = fit$results$par["shape"]),
                         
                         Gumbel = dgumbel(x_grid,
                                          loc = fit$results$par["location"],
                                          scale = fit$results$par["scale"]),
                         
                         GP = {
                           sigma <- fit$results$par["scale"]
                           xi <- fit$results$par["shape"]
                           u <- threshold
                           
                           y <- x_grid[x_grid > u]
                           dens <- (1 / sigma) * (1 + xi * (y - u)/sigma)^(-1/xi - 1)
                           
                           full <- rep(NA, length(x_grid))
                           full[x_grid > u] <- dens
                           full
                         }
  )
  
  df <- data.frame(x = x_grid, density = density_vals)
  
  ggplot() +
    geom_histogram(aes(x = data, y = after_stat(density)),
                   bins = 30, fill = "lightgray", color = "black") +
    geom_line(data = df, aes(x = x, y = density),
              color = "blue", linewidth = 1.2) +
    labs(title = paste("Histogram with", dist_name, "fit"),
         x = "Wind Speed", y = "Density") +
    theme_minimal()
}

# Apply the function to all three models
plot_hist_with_fit(x, fit_list$GEV, "GEV")
plot_hist_with_fit(x, fit_list$Gumbel, "Gumbel")
plot_hist_with_fit(declustered_exceedances$wind_speed, fit_list$GP, "GP", threshold)

# ---------------------------------------------
# B. QQ Plots 

# Create a generic function that can be applied to all models
plot_qq <- function(data, fit, dist_name, threshold = NULL) {
  
  n <- length(data)
  probs <- (1:n) / (n + 1)
  data_sorted <- sort(data)
  
  q_theoretical <- switch(dist_name,
                          GEV = qgev(probs,
                                     loc = fit$results$par["location"],
                                     scale = fit$results$par["scale"],
                                     shape = fit$results$par["shape"]),
                          
                          Gumbel = qgumbel(probs,
                                           loc = fit$results$par["location"],
                                           scale = fit$results$par["scale"]),
                          
                          GP = {
                            sigma <- fit$results$par["scale"]
                            xi <- fit$results$par["shape"]
                            u <- threshold
                            
                            u + sigma * ((1 - probs)^(-xi) - 1)/xi
                          }
  )
  
  df <- data.frame(Theoretical = q_theoretical, Observed = data_sorted)
  
  ggplot(df, aes(x = Theoretical, y = Observed)) +
    geom_point(color = "darkblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Q-Q Plot:", dist_name),
         x = "Theoretical Quantiles",
         y = "Observed Quantiles") +
    theme_minimal()
}


# Apply to all models
plot_qq(x, fit_list$GEV, "GEV")
plot_qq(x, fit_list$Gumbel, "Gumbel")
plot_qq(declustered_exceedances$wind_speed, fit_list$GP, "GP", threshold)


# --------------------------------------------------------
# C. PP Plots

# Create a generic function that can be applied to all models
plot_pp <- function(data, fit, dist_name, threshold = NULL) {
  
  data_sorted <- sort(data)
  n <- length(data)
  emp_probs <- (1:n) / (n + 1)
  
  model_probs <- calc_cdf(data_sorted, fit, dist_name)
  
  df <- data.frame(Empirical = emp_probs, Model = model_probs)
  
  ggplot(df, aes(x = Empirical, y = Model)) +
    geom_point(color = "purple") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "red") +
    labs(title = paste("P-P Plot:", dist_name),
         x = "Empirical Probabilities",
         y = "Model Probabilities") +
    theme_minimal()
}

# Apply to all models
plot_pp(x, fit_list$GEV, "GEV")
plot_pp(x, fit_list$Gumbel, "Gumbel")
plot_pp(declustered_exceedances$wind_speed, fit_list$GP, "GP", threshold)

#----------------------------------------------------
# D. Optional - can also plot the figures created by the existing models - shows many plots together

plot(fit_list$GEV)
plot(fit_list$Gumbel)
plot(fit_list$GP)



#-------------------------------------
# 8. Return levels
# -------------------------------

# A. Define return periods of interest (years) and all fitted distributions

return_periods <- c(1.01, 2, 5, 10, 20, 50, 100, 200)
all_distributions <- names(fit_list)

# Store GP threshold percentile inside the fit object for reference
fit_list$GP$threshold_percentile <- gp_info$percentile

# Function to calculate return values for any distribution
calc_return_value_fast <- function(dist_name, T, fit_obj, total_hours = nrow(wind_data)) {
  
  if(dist_name %in% c("GEV", "Gumbel")) {
    # For annual maxima distributions: use built-in return.level function
    rl <- as.numeric(return.level(fit_obj, return.period = T))
    return(list(ReturnValue = rl, Threshold = NA, ThresholdPercentile = NA)) # GP-specific fields NA
    
  } else if(dist_name == "GP") {
    # For Peak-Over-Threshold (GP) distribution
    
    # Extract parameters from fitted object
    thresh <- fit_obj$threshold
    sigma <- fit_obj$results$par["scale"]
    xi <- fit_obj$results$par["shape"]
    
    # Number of declustered exceedances in dataset
    n_exc <- nrow(declustered_exceedances)
    
    # Convert to yearly exceedance rate assuming original data hourly
    lambda_per_year <- n_exc / total_hours * 24 * 365 
    # Note: Adjust formula if dataset is not hourly (e.g., 10-min intervals)
    
    # GP return level formula
    if(xi != 0){
      rl <- thresh + sigma / xi * ((lambda_per_year * T)^xi - 1)
    } else {
      rl <- thresh + sigma * log(lambda_per_year * T)
    }
    
    # Return also threshold info
    percentile <- fit_obj$threshold_percentile
    
    return(list(ReturnValue = rl, Threshold = thresh, ThresholdPercentile = percentile))
  }
}

#--------------------------------
# B. Create a table for all combinations of return periods and distributions
return_values_table <- expand.grid(
  Distribution = names(fit_list), # Distribution names
  ReturnPeriod = return_periods   # Return periods
)

# Compute return values for each distribution and return period combination
return_values_list <- mapply(
  function(dist_name, T) {
    fit_obj <- fit_list[[dist_name]]
    calc_return_value_fast(dist_name, T, fit_obj, total_hours = nrow(wind_data))
  },
  dist_name = as.character(return_values_table$Distribution),
  T = return_values_table$ReturnPeriod,
  SIMPLIFY = FALSE
)

# Extract return values, thresholds, and percentiles into the table
return_values_table$ReturnValue <- sapply(return_values_list, function(x) x$ReturnValue)
return_values_table$Threshold <- sapply(return_values_list, function(x) x$Threshold)
return_values_table$ThresholdPercentile <- sapply(return_values_list, function(x) x$ThresholdPercentile)

# C. Optionally pivot table to wide format (rows = return periods, columns = distributions)
return_values_wide <- return_values_table %>%
  pivot_wider(names_from = Distribution, values_from = ReturnValue)

# E. Inspect and save results
print(return_values_table)  # Long format with all info
print(return_values_wide)   # Wide format for easier comparison

# Export wide table as CSV
write.csv(return_values_wide, 'output/return_values.csv')





# -------------------------------
# 9. Compute PARAMETRIC BOOTSTRAP CI for best-fit distribution (either Gumbel, GEV, or GP)

# Look at lecture slides on canvas - parametric bootstrap
# Read annotations
# -------------------------------

best_fit_name <- "GP"   # options: "Gumbel", "GEV", "GP"
n_boot <- 500            # number of bootstrap samples

# -----------------------------------------------------------------------
#  A. Generic parametric bootstrap (GEV + Gumbel + GP)
bootstrap_return_parametric <- function(data, return_period, fit, dist_name,
                                        threshold = NULL, total_hours = NULL) {
  rl <- NA
  try({
    
    # i) SIMULATE DATA
    if (dist_name %in% c("GEV", "Gumbel")) {
      
      n <- length(data) # Get the number of data points in the original sample
      
      # Extract fitted model parameters for location, scale, and shape
      loc   <- fit$results$par["location"]
      scale <- fit$results$par["scale"]
      
      # For Gumbel, shape parameter is 0; for GEV, use estimated shape
      shape <- ifelse(dist_name == "Gumbel", 0, fit$results$par["shape"])
      
      u <- runif(n) # Generate n random uniform values between 0 and 1 for simulation
      
      # Simulate data from the inverse CDF of the selected distribution:
      # If shape != 0, simulate from GEV; else simulate from Gumbel
      d_sim <- if (shape != 0) {
       
         # GEV inverse CDF formula applied to uniform samples
        loc + (scale / shape) * ((-log(u))^(-shape) - 1)
        # Gumbel inverse CDF formula applied to uniform samples
      } else {
        loc - scale * log(-log(u))   # Gumbel
      }
      
      # Refit the specified extreme value distribution to the simulated data
      fit_sim <- fevd(d_sim, type = dist_name)
      
      # Calculate the return level (quantile) for the given return period from the refit model
      rl <- return.level(fit_sim, return.period = return_period)
      
    } else if (dist_name == "GP") {
      
      # Ensure required parameters are provided
      if (is.null(threshold) || is.null(total_hours)) {
        stop("Threshold and total_hours required for GP")
      }
      
      # Extract fitted GP parameters
      sigma <- fit$results$par["scale"] # scale parameter
      xi    <- fit$results$par["shape"] # shape parameter
      thresh <- threshold  # threshold for exceedances
      
      n_exc <- length(data) # number of exceedances in the observed data
      
      # -------------------------------
      # Frequency simulation (Poisson)
      # -------------------------------
      # Simulate the number of exceedances in the bootstrap sample
      n_exc_sim <- rpois(1, lambda = n_exc)
      
      # If too few simulated exceedances, return NA to avoid unstable fits
      if (n_exc_sim < 5) return(NA)
      
      # -------------------------------
      # Magnitude simulation (GP)
      # -------------------------------
      u <- runif(n_exc_sim) # uniform random values for inverse CDF simulation
      
      # Simulate GP exceedances using the inverse CDF
      d_sim <- if (xi != 0) {
        thresh + (sigma / xi) * (u^(-xi) - 1) # GP formula for xi != 0
      } else {
        thresh - sigma * log(u) # GP formula for xi = 0
      }
      
      # Refit GP to the simulated exceedances
      fit_sim <- fevd(d_sim, type = "GP", threshold = thresh)
      
      # Extract parameters from refit
      sigma_sim <- fit_sim$results$par["scale"]
      xi_sim    <- fit_sim$results$par["shape"]
      
     
      # -------------------------------
      # Update exceedance rate
      # -------------------------------
      # Convert number of exceedances to expected annual rate
      lambda_sim <- n_exc_sim / total_hours * 24 * 365
      
      # -------------------------------
      # Compute return level
      # -------------------------------
      rl <- if (xi_sim != 0) {
        thresh + sigma_sim / xi_sim * ((lambda_sim * return_period)^xi_sim - 1)
      } else {
        thresh + sigma_sim * log(lambda_sim * return_period)
      }
      
    } else {
      stop("Unsupported distribution") # Catch invalid distribution input
    }
    
  }, silent = TRUE) # End of try block; suppress errors to avoid stopping bootstrap

  return(rl) # Return the computed return level from this bootstrap replicate
}

# -----------------------------------------------------------------------

# B. Select the best-fit dataset

# Choose which dataset to use for bootstrapping depending on the selected distribution
data_to_boot <- switch(best_fit_name,
                       GEV = x,        # Use annual maxima for GEV
                       Gumbel = x,     # Use annual maxima for Gumbel
                       GP = declustered_exceedances$wind_speed) # Use declustered exceedances for GP

# For GP, retrieve the threshold used in the original fit; for GEV/Gumbel, threshold not needed
threshold_val <- if(best_fit_name == "GP") gp_info$threshold else NULL

# For GP, calculate total hours of the dataset to convert exceedance counts into a rate per year
# For GEV/Gumbel, this is not required
total_hours_val <- if(best_fit_name == "GP") nrow(wind_data) else NULL


# -----------------------------------------------------------------------


# C. Compute central return levels

# For each return period, compute the "central" return level
# i.e., the estimate from the original fitted model (not from bootstrap yet)
pred_return <- sapply(return_periods, function(T) {
  
  # calc_return_value_fast() is a custom function that calculates
  # the return level for the specified distribution and return period
  # using the parameters from the best-fit model
  calc_return_value_fast(
    best_fit_name,               # distribution type: GEV, Gumbel, or GP
    T,                          # return period of interest (e.g., 10, 50, 100 years)
    fit_list[[best_fit_name]]   # fitted model object for the selected distribution
    )$ReturnValue                # extract the numeric return level value
})



# -----------------------------------------------------------------------

# D. Bootstrap CI - repeat many times

# Initialize a list to store confidence intervals for each return period
ci_results <- list()

# Loop over each return period
for(i in seq_along(return_periods)) {
  T <- return_periods[i]  # current return period (e.g., 10, 50, 100 years)
  
  # Perform parametric bootstrap n_boot times for the current return period
  boot_vals <- replicate(n_boot,
                         bootstrap_return_parametric(
                           data = data_to_boot,             # dataset to simulate from (GEV/Gumbel/GP)
                           return_period = T,               # return period of interest
                           fit = fit_list[[best_fit_name]], # fitted model for the selected distribution
                           dist_name = best_fit_name,       # type of distribution: GEV, Gumbel, GP
                           threshold = threshold_val,       # threshold (for GP only)
                           total_hours = total_hours_val    # total hours (for GP only)
                         )
  )
  
  # Remove any failed bootstrap iterations (NA values)
  boot_vals <- boot_vals[!is.na(boot_vals)]
  
  # Compute the 95% confidence interval from the bootstrap distribution
  # 2.5% and 97.5% quantiles represent the lower and upper bounds
  ci_results[[as.character(T)]] <- quantile(boot_vals, probs = c(0.025, 0.975))
}



# -----------------------------------------------------------------------
# E. Build CI dataframe

# Combine the return periods, central estimates, and bootstrap CIs into a single data frame
ci_df <- data.frame(
  ReturnPeriod = return_periods,  # vector of return periods (e.g., 10, 50, 100 years)
  Estimate = pred_return,         # central return level estimates from the original fitted model
  Lower = sapply(return_periods, function(T) ci_results[[as.character(T)]][1]),  # 2.5% quantile (lower CI)
  Upper = sapply(return_periods, function(T) ci_results[[as.character(T)]][2])   # 97.5% quantile (upper CI)
)

# Display the resulting dataframe in the console
print(ci_df)

# ---------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------
# 10. Plot return levels for the best-fit distribution with computed confidence intervals
# --------------------------------------------------------------------------------------

# Adaptive Return Level Plot - creates plot for one distribution depending on what 
# one you selected as best fit i.e., best_fit_name. Includes Confidence Intervals

# Create function that can be adapted to all three distributions. 
##  Gumbel & GEV were fitted to Annual Max
##  GP was fitted to POT (declustered)

plot_return_levels <- function(best_fit_name, ci_df, annual_max, declustered_exceedances, wind_data, min_return_period = 1) {
  
  # Prepare empirical return periods
  if(best_fit_name == "GP") {
    exc <- declustered_exceedances$wind_speed
    exc_sorted <- sort(exc, decreasing = TRUE)
    n_exc <- length(exc_sorted)
    n_years <- length(unique(format(wind_data$datetime, "%Y")))
    lambda <- n_exc / n_years
    emp_rp <- 1 / ((1:n_exc) / (n_exc + 1) * lambda)
    emp_df <- data.frame(ReturnPeriod = emp_rp, WindSpeed = exc_sorted)
  } else {
    x_sorted <- sort(annual_max, decreasing = TRUE)
    n <- length(x_sorted)
    emp_rp <- (n + 1) / (1:n)
    emp_df <- data.frame(ReturnPeriod = emp_rp, WindSpeed = x_sorted)
  }
  
  # Filter empirical points for min_return_period
  emp_df <- emp_df %>% dplyr::filter(ReturnPeriod >= min_return_period)
  
  # Filter CI data similarly for cleaner plotting (optional)
  ci_df <- ci_df %>% dplyr::filter(ReturnPeriod >= min_return_period)
  
  p <- ggplot() +
    geom_ribbon(data = ci_df,
                aes(x = ReturnPeriod, ymin = Lower, ymax = Upper),
                fill = "lightblue", alpha = 0.3) +
    geom_line(data = ci_df,
              aes(x = ReturnPeriod, y = Estimate),
              color = "blue", linewidth = 1.2) +
    geom_point(data = emp_df,
               aes(x = ReturnPeriod, y = WindSpeed),
               color = "red", shape = 1, size = 2) +
    scale_x_log10() +
    labs(x = "Return Period (Years)",
         y = "Wind Speed (m/s)",
         title = paste("Return Level Plot - Best Fit:", best_fit_name),
         subtitle = "Blue line: fitted return | Light blue: 95% CI | Red: empirical points") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}


# Apply function to plot return levels
plot_return_levels(best_fit_name, ci_df, x, declustered_exceedances, wind_data, min_return_period = 1)

