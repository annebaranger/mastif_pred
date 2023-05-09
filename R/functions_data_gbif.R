#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien BARRERE, Georges KUNSTLER
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 1. Generic functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
#' @author Julien BARRERE
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}



#' Get file from its url and write it on disk, at a specified location. 
#' @param dir.name.in Directory where the file should be written (ex: "data/BWI")
#' @param url.in URL where to download the file.
#' @author Julien BARRERE
get_and_write <- function(dir.name.in, url.in){
  
  # Write directories if they do not exist
  path.in <- strsplit(dir.name.in, "/")[[1]]
  for(i in 1:length(path.in)){
    if(i == 1) path.in_i <- path.in[i]
    else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
    if(!dir.exists(path.in_i)) dir.create(path.in_i)
  }
  
  # Write file on the disk
  url.in_split <- strsplit(url.in, "/")[[1]]
  file.in <- paste(dir.name.in, url.in_split[length(url.in_split)], sep = "/")
  if(!file.exists(file.in)){
    try(GET(url.in, authenticate('guest', ""), write_disk(file.in, overwrite = TRUE)))
    # Specific case of zip file: unzip and delete zip file
    if("zip" %in% strsplit(file.in, split = "\\.")[[1]]){
      unzip(file.in, exdir = dir.name.in, overwrite = T)
      print(paste0("---Getting and unzipping ", file.in))
      unlink(file.in)
    }else{print(paste0("---Getting ", file.in))}
  } 
}


#' Write a table on disk
#' @param table.in dataframe to write on the disk
#' @param file.in Name (and path) of the file on the disk
#' @author Julien BARRERE
write_on_disk <- function(table.in, file.in){
  create_dir_if_needed(file.in)
  write.table(table.in, file = file.in, row.names = F)
  return(file.in)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 2. Functions to get GBIF data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Get the species for which to get GBIF presence data
#' @return A character vector containing species latin name
get_species <- function(){
  c("Abies alba", "Acer nigrum", "Carpinus betulus" , "Cedrus atlantica", "Fagus sylvatica", 
    "Frangula alnus", "Laurus nobilis", "Picea abies", "Pinus pinaster", "Pinus pinea", "Pinus sylvestris",
    "Quercus ilex", "Quercus petraea", "Quercus pubescens", "Quercus robur","Sorbus aucuparia")
  # c("Abies alba", "Acer campestre", "Acer pseudoplatanus", "Alnus glutinosa", "Betula pendula", 
  #   "Carpinus betulus" , "Castanea sativa", "Erica arborea", "Eucalyptus camaldulensis", "Fagus sylvatica", 
  #   "Fraxinus angustifolia", "Fraxinus excelsior", "Juniperus thurifera", "Myrica faya", "Olea europaea", 
  #   "Picea abies", "Pinus canariensis", "Pinus halepensis", "Pinus nigra", "Pinus pinaster", "Pinus pinea", 
  #   "Pinus radiata", "Pinus sylvestris", "Pinus uncinata", "Populus nigra", "Populus tremula", "Prunus avium", 
  #   "Pseudotsuga menziesii", "Quercus faginea", "Quercus ilex", "Quercus petraea", "Quercus pubescens", 
  #   "Quercus pyrenaica", "Quercus robur", "Quercus suber", "Robinia pseudoacacia", "Salix caprea", "Ulmus minor")
}




#' Get the taxon keys from GBIF
#' @param sp_vec vector of species
#' @author Georges Kunstler
get_gbif_taxon_keys <- function(sp_vec){
  # add original name back into data.frame
  gbif_taxon_keys <- get_gbifid_(sp_vec, method="backbone") %>% imap(~ .x %>% mutate(original_sciname = .y)) %>% 
    # combine all data.frames into one
    bind_rows() %T>% 
    # save as side effect for you to inspect if you want
    readr::write_tsv(path = file.path("output", "all_matches.tsv")) %>% 
    # get only accepted and matched names
    filter(matchtype == "EXACT" & status == "ACCEPTED") %>% 
    # remove anything that might have matched to a non-plant
    filter(kingdom == "Plantae") %>% 
    pull(usagekey) 
  return(gbif_taxon_keys)  
}



#' Get the requested data from GBIF
#' @param gbif_taxon_keys List of all taxon keys fo the species for which to get data
#' @param user gbif username
#' @param pwd gbif password
#' @param email email for gbif
#' @author Georges Kunstler, Julien Barrere
get_data_gbif <- function(gbif_taxon_keys, user, pwd, email){
  
  # use matched gbif_taxon_keys from above to send request download
  request_download <- occ_download(
    pred_in("taxonKey", gbif_taxon_keys),
    pred_in("basisOfRecord", c('HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_or(
      pred_not(pred("establishmentMeans","INTRODUCED")),
      pred_not(pred_notnull("establishmentMeans"))
    ),
    pred_or(
      pred_not(pred("establishmentMeans","INVASIVE")),
      pred_not(pred_notnull("establishmentMeans"))
    ),
    pred_or(
      pred_not(pred("establishmentMeans","NATURALISED")),
      pred_not(pred_notnull("establishmentMeans"))
    ),
    pred_or(
      pred_lt("coordinateUncertaintyInMeters",10000),
      pred_not(pred_notnull("coordinateUncertaintyInMeters"))
    ),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email
  )
  
  # Wait for the download to complete
  print("Downloading files")
  occ_download_wait(request_download[[1]])
  
  # In case the output dir is not created, create it
  if(!dir.exists("output")) dir.create("output")
  
  # Once download is complete, get the data
  out <- occ_download_get(request_download[[1]], path = "output", overwrite = TRUE) %>% 
    occ_download_import %>%
    setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
    filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
    filter(!coordinateuncertaintyinmeters %in% c(999, 9999)) %>% 
    filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
    cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
    cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
    cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
    cc_sea() %>% # remove from ocean 
    distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>%
    glimpse()
  
  # Return dataset
  return(out)
}



#' Function to filter gbif observations outside of the species native continent
#' @param data_gbif presence data of gbif with coordinates
filter_data_gbif <- function(data_gbif){
  
  # Create a table with coordinates limits for each species
  coord_lim <- data.frame(species = unique(data_gbif$species)) %>%
    mutate(
      block = case_when(species %in% c("Pseudotsuga menziesii", "Robinia pseudoacacia") ~ "American native block", 
                        species == "Eucalyptus camaldulensis" ~ "Australian native block", 
                        TRUE ~ "Eurasian native block"), 
      coord = case_when(block == "American native block" ~ "0_80_-160_-30", 
                        block == "Australian native block" ~ "-45_-8_110_154", 
                        block == "Eurasian native block" ~ "25_75_-10_160")
    ) %>%
    separate(col = "coord", into = c("lat.min", "lat.max", "long.min", "long.max"), sep = "_") %>%
    mutate(long.min = as.numeric(long.min), long.max = as.numeric(long.max), 
           lat.min = as.numeric(lat.min), lat.max = as.numeric(lat.max))
  
  
  # Filter data gbif
  out <- data_gbif %>%
    left_join(coord_lim, by = "species") %>%
    filter(decimallongitude >= long.min & decimallongitude <= long.max &
             decimallatitude >= lat.min & decimallatitude <= lat.max)
  
  
  # Return the output
  return(out)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 3. Functions to manage CHELSA climatic data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Download chelsa files
#' @param bioclim numeric vector of the bioclim codes to extract
#' @param path character: directory where the data should be stored
#' @return 
#' @author G. Kunstler, J. Barrere
download_CHELSA <- function(bioclim, path){
  
  # Loop on all bioclim codes
  for(b in bioclim){
    if(! b %in% c(1:19)) stop("Argument bioclim needs to be in range 1:19")
    # Name of the file to extract
    filename.ymv <- paste("CHELSA_bio10", b, "land.7z", sep = "_")
    # Name of the url to extract
    url.ymv <- paste("https://zenodo.org/record/4996318/files", 
                     filename.ymv, sep = "/")
    # Download the file
    get_and_write(path, url.ymv)
  }
  
  # Identify archived files
  chelsa_files_archived <- paste(path, list.files(path), sep = "/")
  
  # Loop on all chelsa files archived to unzip them
  for(i in 1:length(chelsa_files_archived)){
    print(paste0("--- Unarchiving file ", chelsa_files_archived[i]))
    # Extract the file(s) from the archive
    archive_extract(chelsa_files_archived[i], dir = path)
    # Remove the archive
    unlink(chelsa_files_archived[i])
  }
  
  
  return(paste(path, list.files(path), sep = "/"))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 4. Functions to get disturbance related climate data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to download GLOBSNOW and unzip snow cover data
#' @param dir.in directory where to save the files
get_GlobSnow <- function(dir.in){
  
  # - Dates for which data are available
  dates.in <- (expand.grid(year = as.character(c(1979:2017)), 
                           month = c(paste0("0", c(1:9)), as.character(c(10:12)))) %>%
                 mutate(year_month = paste0(year, month)))$year_month
  
  # - URL to download the files
  url.in <- paste0(
    "https://www.globsnow.info/swe/archive_v3.0/L3B_monthly_SWE/NetCDF4/", 
    dates.in, 
    "_northern_hemisphere_monthly_swe_0.25grid.nc"
  )
  
  # - Apply function to download all files
  for(j in 1:length(url.in)) get_and_write(dir.in, url.in[j])
  
  # Return the name of the extracted files
  return(paste(dir.in, list.files(dir.in), sep = "/"))
}


#' Function to download wind speed data from the global Wind atlas
#' @param dir.in directory where to save the file
get_GlobalWindAtlas <- function(dir.in){
  
  # - URL to download the file
  url.in <- "https://figshare.com/ndownloader/files/17247017"
  
  # - Name of the file to download
  file.in <- paste(dir.in, "gwa3_250_wind-speed_100m.tif", sep = "/")
  
  # - Create directory if needed
  create_dir_if_needed(file.in)
  
  # - Download the file
  try(GET(url.in, write_disk(file.in, overwrite = TRUE)))
  
  # - Return the name of the extracted file
  return(paste(dir.in, list.files(dir.in), sep = "/"))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 5. Functions to merge all data together ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Extract climate variables for each location 
#' @param chelsa_files character vector containing all chelsa files
#' @param data_gbif presence data of gbif with coordinates
extract_climate_for_gbif <- function(chelsa_files, data_gbif){
  
  out <- data_gbif
  
  # Extract mean annual temperature
  chelsa_file_mat <- grep("bio10_1.tif", chelsa_files, value = TRUE)
  raster_mat <- terra::rast(chelsa_file_mat)
  out$mat <- as.numeric(terra::extract(raster_mat,
                                       cbind(out$decimallongitude,
                                             out$decimallatitude))[, 1])/10
  
  # Extract mean annual precipitation
  chelsa_file_map <- grep("bio10_12.tif", chelsa_files, value = TRUE)
  raster_map <- terra::rast(chelsa_file_map)
  out$map <- as.numeric(terra::extract(raster_map,
                                       cbind(out$decimallongitude,
                                             out$decimallatitude))[, 1])
  
  # Extract min temperature
  chelsa_file_tmin <- grep("bio10_6.tif", chelsa_files, value = TRUE)
  raster_tmin <- terra::rast(chelsa_file_tmin)
  out$tmin <- as.numeric(terra::extract(raster_tmin,
                                        cbind(out$decimallongitude,
                                             out$decimallatitude))[, 1])/10
  
  
  # - Finish formatting
  out <- out %>%
    group_by(species) %>%
    summarize(mat.low = quantile(mat, probs = 0.025, na.rm = TRUE), 
              mat.high = quantile(mat, probs = 0.975, na.rm = TRUE), 
              mat = mean(mat, na.rm = TRUE), 
              map.low = quantile(map, probs = 0.025, na.rm = TRUE), 
              map.high = quantile(map, probs = 0.975, na.rm = TRUE), 
              map = mean(map, na.rm = TRUE), 
              tmin.low = quantile(tmin, probs = 0.025, na.rm = TRUE), 
              tmin.high = quantile(tmin, probs = 0.975, na.rm = TRUE), 
              tmin = mean(tmin, na.rm = TRUE))
  
  return(out)
  
}



#' Extract climate variables for each location 
#' @param chelsa_files character vector containing all chelsa files
#' @param data_gbif presence data of gbif with coordinates
extract_disturbance_index_for_gbif <- function(
  fireweatherindex_files, globalwindatlas_files, globsnow_files, data_gbif){
  
  # - Initialize output
  out <- data_gbif
  
  
  # - Extract wind data
  
  # --- Convert file into raster
  raster_windspeed <- terra::rast(globalwindatlas_files)
  # --- Extract data from raster for each gbif point
  out$windspeed <- as.numeric(
    terra::extract(raster_windspeed, cbind(out$decimallongitude, out$decimallatitude))[, 1])
  
  
  
  # - Compute mean fire weather index
  
  # --- Identify the years for which we have FWI data
  years <- unique(substr(gsub(".+FWI\\_", "", fireweatherindex_files), 1, 4))
  # --- Loop to calculate the mean FWI for each year
  for(i in 1:length(years)){
    # -- List of files for year i
    files.i <- fireweatherindex_files[grep(paste0("FWI_", years[i]), fireweatherindex_files)]
    # Initialize the raster for year i with the first day
    raster_year.i <- terra::rast(files.i[1])
    # Loop on all other days
    for(j in 2:length(files.i)) raster_year.i <- c(raster_year.i, terra::rast(files.i[j]))
    # Average the fire weather index over all days
    raster_year.i <- app(raster_year.i, mean)
    # Convert longitude from 0-360 to -180-180 scale
    raster_year.i <- rotate(raster_year.i)
    # If i == 1, initialize the final raster, otherwise add it to final raster
    if(i == 1) raster_fire <- raster_year.i
    else raster_fire <- c(raster_fire, raster_year.i)
  }
  # --- Combine raster of all years into one raster
  raster_fire <- app(raster_fire, mean)
  # --- Extract data from raster for each gbif point
  out$fwi <- as.numeric(
    terra::extract(raster_fire, cbind(out$decimallongitude, out$decimallatitude))[, 1])
  
  
  # - Compute mean snow water equivalent
  
  # --- Initialize with the raster of the first month
  raster_snow <- terra::rast(globsnow_files[1])
  # --- Assemble all raster files into one multi-layer raster
  for(i in 2:length(globsnow_files)){
    # Check that the file is large enough (there are empty rasters in the directory)
    if(file.info(globsnow_files[i])$size > 1000){
      # Check that the file corresponds to a winter month
      if(as.numeric(substr(globsnow_files[i], 19, 20)) %in% c(12, 1, 2, 3)){
        raster_snow <- c(raster_snow, terra::rast(globsnow_files[i]))
      }
    }
  }
  # --- Average the snow water equivalent over all layers
  raster_snow <- terra::app(raster_snow, mean)
  raster_snow <- project(x = raster_snow, y = "EPSG:4326", method = "bilinear")
  # --- Extract data from raster for each gbif point
  out$swe <- as.numeric(
    terra::extract(raster_snow, cbind(out$decimallongitude, out$decimallatitude))[, 1])
  
  
  # - Finish formatting
  out <- out %>%
    group_by(species) %>%
    summarize(windspeed.low = quantile(windspeed, probs = 0.025, na.rm = TRUE), 
              windspeed.high = quantile(windspeed, probs = 0.975, na.rm = TRUE), 
              windspeed = mean(windspeed, na.rm = TRUE), 
              fwi.low = quantile(fwi, probs = 0.025, na.rm = TRUE), 
              fwi.high = quantile(fwi, probs = 0.975, na.rm = TRUE), 
              fwi = mean(fwi, na.rm = TRUE), 
              swe.low = quantile(swe, probs = 0.025, na.rm = TRUE), 
              swe.high = quantile(swe, probs = 0.975, na.rm = TRUE), 
              swe = mean(swe, na.rm = TRUE))
  
  return(out)
  
}



