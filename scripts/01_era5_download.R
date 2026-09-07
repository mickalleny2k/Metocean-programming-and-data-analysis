###########################################################################
# GG6581 Assignment
# 05_era5_download.R
# Download era5 data for Lí Ban (Mermaid Saint) site
###########################################################################

## Aim of script:
# Write a function to automatically download era5 data from CDS, saving variables as monthly files, so to not
# exceed file size limitations. Save all variables into separate subfolders.

# 1. install required package
install.packages("ecmwfr")
install.packages("plyr")
install.packages("sf")

# 2. load required packages
library("ecmwfr")
library("plyr")
library("sf")

# 3. Generate access to CDS API using the ecmwfr package. This line only needs to be run once on your computer.
ecmwfr::wf_set_key(key = "3dd19e02-dad8-42a3-b1fd-dec150098abb", # between "_" insert your key from Climate Data Store
                   user = "95694862@umail.ucc.ie") # between "_" insert your email used for your profile in Climate Data Store e.g., youremail@gmail.com

# 4. Create a function to download multiple ERA5 variables from the dataset 'Reanalysis ERA5 Single Levels'
download_era5_single_levels_nc <- function(shape = NULL,
                                             lat,
                                             lon,
                                             buffer = 0.5,
                                             variable = c("10m_u_component_of_wind",
                                                          "10m_v_component_of_wind",
                                                          "100m_u_component_of_wind",
                                                          "100m_v_component_of_wind",
                                                          "mean_sea_level_pressure",
                                                          "sea_surface_temperature",
                                                          "mean_wave_period",
                                                          "mean_wave_direction",
                                                          "peak_wave_period",
                                                          "significant_height_of_combined_wind_waves_and_swell"),#add more variable as per ERA5 documentation
                                             year = 2022,
                                             month = 1:12,
                                             site  = "test",
                                             user = NULL,
                                             era5_dataset = "reanalysis-era5-single-levels",
                                             path = ".",
                                             subfolder = "era5_downloads") {
  
  # Main download folder
  main_path <- file.path(path, subfolder)
  if (!dir.exists(main_path)) {
    dir.create(main_path, recursive = TRUE)
  }
  
  time <- sprintf("%02d:00", 0:23)
  timestep <- "hourly"
  
  if (!is.null(shape)) {
    bbox <- shape |> sf::st_bbox()
    area <- c(plyr::round_any(bbox["ymax"], 0.5, ceiling),
              plyr::round_any(bbox["xmin"], 0.5, floor),
              plyr::round_any(bbox["ymin"], 0.5, floor),
              plyr::round_any(bbox["xmax"], 0.5, ceiling))
  } else {
    area <- c(plyr::round_any(lat + buffer, 0.5, ceiling),
              plyr::round_any(lon - buffer, 0.5, floor),
              plyr::round_any(lat - buffer, 0.5, floor),
              plyr::round_any(lon + buffer, 0.5, ceiling))
  }
  
  request_list <- list()
  
  for (v in variable) {
    # Create folder for this variable inside the subfolder
    var_folder <- v
    dir.create(file.path(main_path, var_folder), recursive = TRUE, showWarnings = FALSE)
    
    for (y in year) {
      for (m in month) {
        list_name <- paste0(v, "_", y, "_", m)
        
        # TARGET is relative to main_path (no subfolder)
        target_file <- file.path(
          var_folder,
          paste0(era5_dataset, "_", v, "_", timestep, "_", y, "_", m, "_", site, ".nc")
        )
        
        request_list[[list_name]] <- list(
          dataset_short_name = era5_dataset,
          product_type       = "reanalysis",
          download_format    = "unarchived",
          data_format        = "netcdf",
          variable           = v,
          year               = y,
          month              = m,
          day                = sprintf("%02d", 1:31),
          time               = time,
          area               = area,
          target             = target_file
        )
      }
    }
  }
  
  # Check requests
  lapply(request_list, ecmwfr::wf_check_request)
  
  # Batch download
  req <- ecmwfr::wf_request_batch(
    request_list = request_list,
    path = main_path,
    user = user,
    workers = length(request_list)
  )
  
  return(req)
}

### The function has been created above.

# 5. Define your input criteria into your newly created function
lat <- 51.795
lon <- -7.5
year <- 2024
month <- 1:2
variable <- c("100m_u_component_of_wind", "100m_v_component_of_wind")
path <- "data/"
user <- "95694862@umail.ucc.ie" # between "_" insert your specific user email linked to Climate Data Store


# 6. Use the function to download the required NetCDF files.
files <- download_era5_single_levels_nc(
  lat = lat,
  lon = lon,
  year = year,
  month = month,
  variable = variable,
  path = path,
  user = user 
)

files
