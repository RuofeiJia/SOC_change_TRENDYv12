# load the area data for each model and generate masks for area in the study domain
rm(list = ls(all = TRUE))

library(ncdf4)
library(terra)
library(raster)
library(dplyr)
library(tictoc)
library(sp)

# working directory should be the this the folder "code"

#---------1. make raster of area fraction used for data-driven SOC change estimation (cropped) -------

# download information on study domain and area fractions of forests and grasslands
# generated from repository https://github.com/RuofeiJia/SOCTimeSeries

# load area fraction of forests and grasslands in the study domain
load("GlobalSOCChange/results/proportion_decade_df.rda")

# Load lists of pixels within study domain
load("GlobalSOCChange/results/predicted_inputs_brms1_37k12chain.rda")

# make sum of young forest, old forest, grasslands fractions predicted 
# find pixels that are predicted
yf_pred <- input_youngforest_pred %>%
  select(x,y) %>% mutate(yf_predicted=1)
of_pred <- input_oldforest_pred %>%
  select(x,y) %>% mutate(of_predicted=1)
g_pred <- input_grassland_pred %>%
  select(x,y) %>% mutate(g_predicted=1)

# calculate fraction of predicted pixels
predicted_fraction_df <- proportion_decade_df %>%
  left_join(yf_pred, by = c("x", "y")) %>%
  left_join(of_pred, by = c("x", "y")) %>%
  left_join(g_pred, by = c("x", "y")) %>%
  # convert NA to 0s in predicted flag
  dplyr::mutate(across(c(yf_predicted, of_predicted, g_predicted),
                       ~replace(., is.na(.), 0))) %>%
  # calculate predicted fraction
  mutate(
    value = avg_youngforest1992_2020*yf_predicted+avg_oldforest1992_2020*of_predicted+grassland1992_2020*g_predicted
  ) %>%
  dplyr::select(x,y,value)

# # check predicted fraction 
# max(predicted_fraction_df$value) #1
# min(predicted_fraction_df$value) #0
# ggplot(predicted_fraction_df) +
#   geom_histogram(aes(x=value))
# # looks good

# fill the fraction values to full grid
# construct full grid
full_grid <- expand.grid(x=seq(-179.75, 179.75, by=0.5), 
                         y=seq(-89.75, 89.75, by=0.5))
# map fraction vlaues to full grid
predicted_fraction_gridded <- merge(full_grid, predicted_fraction_df, 
                                    by=c("x","y"), all.x=T) %>%
  # convert all NA values into 0
  dplyr::mutate(value = ifelse(is.na(value), 0, value))

# make into spatial dataframe
# convert to SpatialPixelsDataFrame
coordinates(predicted_fraction_gridded) <- ~x+y
gridded(predicted_fraction_gridded) <- TRUE
# convert to raster
predicted_fraction_raster <- raster(predicted_fraction_gridded)
# optional: check
plot(predicted_fraction_raster)
# write raster
writeRaster(predicted_fraction_raster,'../data/prediction_fractions/predicted_fraction_raster.tif',overwrite=TRUE,progress='text')

#---------2. resample natural fraction raster to model resolutions -------------
# load processed area grids for each model
all_files <- list.files("../data/processed_area", pattern = "_area\\.nc$", full.names = TRUE)

# convert each model resolution to rasters
raster_list <- list()
for (dir in all_files) {
  # extract model name
  model <- sub(".*/([^/_]+)_area\\.nc$", "\\1", dir)
  print(model)
  
  # exclude CLM5.0 and CABLEPOP (not included in comparison)
  if (model=="CLM5.0" | model=="CABLEPOP") {
    next
  }
  
  # load raster
  r <- raster(dir)
  
  # reassign value of raster to be 1
  r[] <- 1
  
  # rename raster
  assign(paste0(model, "_area"), r)

  # append to named list
  raster_list[[model]] <- r
}

# make list of unique resolutions
unique_rasters <- unique(raster_list)
names(unique_rasters) <- c(1,2,3,4,5,6,7,8,9)

# group models by the same resolution
group_list <- data.frame()
# for each model
for (i in 1:length(raster_list) ) {
  print(paste0("i ", i))
  # extract model name and resolution raster
  model <- names(raster_list)[i]
  model_raster <- raster_list[[i]]
  
  # compare model raster with each unique raster
  for (j in 1:length(unique_rasters) ) {
    print(paste0("j ", j))
    # extract j-th unique resolution
    unique_raster <- unique_rasters[[j]]

    # when find the matching unique resolution, assign group index and continue to next model
    if ( identical(model_raster, unique_raster) ) {
      # assign group to model in group list
      group_list <- rbind(
        group_list,
        data.frame(model=model, polygon_index=j, stringsAsFactors = FALSE)
        )
      #group_list[[model]] <- j
      print(paste0("appended ", model, " to unique index ", j))
      next
    }
  }
}

# save group_list
write.table(group_list, "../data/prediction_fractions/unique_resolution_list.txt", 
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

# convert predicted_fraction_raster to terra::rast (SpatRaster) for resampling
predicted_fraction_rast <- rast(predicted_fraction_raster)

# resample predcited fraction raster to each resolution so polygons might not be necessary
# for each unique model resolution
for (j in 1:length(unique_rasters)) {
  # extract j-th unique resolution
  unique_rast <- rast(unique_rasters[[j]])
  
  # resample predicted fraction (0.5 res) to unique resolution
  tic()
  resampled_fraction <- resample(predicted_fraction_rast, unique_rast, method="bilinear")
  toc(); print(paste0("resolution ", j ," resampled"))
  
  # save as .nc
  tic()
  save_dir <- paste0("../data/prediction_fractions/fraction_resolution", j, ".nc")
  writeCDF(resampled_fraction, filename = save_dir, overwrite = TRUE)
  toc(); print(paste0("resolution ", j ," saved"))
}
