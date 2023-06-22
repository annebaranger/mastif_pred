# _targets.R file
library(targets)
# lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
source("R/functions_data_gbif.R")
packages.in <- c("dplyr", "ggplot2", "readxl", "purrr", "readr", "magrittr", "rgbif", "taxize",
                 "CoordinateCleaner", "RCurl", "httr", "archive", "terra", "cowplot", "R.utils", 
                 "sf", "ggspatial", "rnaturalearth", "rnaturalearthdata", "tidyr", "scales")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
options(tidyverse.quiet = TRUE)
sf::sf_use_s2(FALSE)
tar_option_set(packages = packages.in)

species.eu=tar_read(species.eu,store = "target_data")
species.am=tar_read(species.am,store = "target_data")

list(
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - GBIF SPECIES DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Give user and password to access gbif data
  tar_target(user, "anbaranger"), 
  tar_target(pwd, "ooDait9C.choco"), 
  tar_target(email, "annebarang@gmail.com"), 
  
  # Make and send request to download data from GBIF
  tar_target(df.species,get_species(species.am,species.eu)),
  tar_target(sp_vec, df.species$TaxonName), 
  tar_target(gbif_taxon_keys, get_gbif_taxon_keys(sp_vec)), 
  tar_target(data_gbif, get_data_gbif(gbif_taxon_keys, user, pwd, email)),
  
  # Remove observations out of native species range
  tar_target(data_gbif_filtered, filter_data_gbif(data_gbif,df.species)),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - CHELSA CLIMATE DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Download CHELSA climatic data
  tar_target(chelsa_files, download_CHELSA(bioclim = c(1, 6, 12), path = "data/CHELSA"), format = "file"),

  # Compute optimal climatic conditions for each species
  tar_target(meanClimate_species, extract_climate_for_gbif(chelsa_files, data_gbif_filtered)), 
  
  # Export file
  tar_target(meanClimate_species_file, write_on_disk(
    meanClimate_species, "output/sp_gbif_climate.csv"), format = "file"), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - DISTURBANCE INDEX DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # # Download GlobSnow data on snow cover
  # tar_target(globsnow_files, get_GlobSnow("data/GlobSnow"), format = "file"),
  # 
  # # Download Gloabal Wind atlas data on wind speed
  # tar_target(globalwindatlas_files, list.files("data/GlobalWindAtlas", full.names = TRUE), format = "file"),
  # 
  # # Fire weather index files at global scale (need to download it before running the script)
  # tar_target(fireweatherindex_files, list.files("data/FireWeatherIndex", full.names = TRUE), format = "file"),
  # 
  # # Compute mean disturbance index for each species
  # tar_target(meanDistIndex_species, extract_disturbance_index_for_gbif(
  #   fireweatherindex_files, globalwindatlas_files, globsnow_files, data_gbif_filtered)), 
  # 
  # # Export file
  # tar_target(meanDistIndex_species_file, write_on_disk(
  #   meanDistIndex_species, "output/sp_gbif_disturbance.csv"), format = "file"), 
  # 
  # 
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - PLOTS ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # tar_target(fig_sp_distribution, plot_species_distribution(data_gbif_filtered, "fig/species_distribution.jpg"), 
  #            format = "file"), 
#   tar_target(fig_sp_distribution_density, plot_species_distribution_density(
#     data_gbif_filtered, "fig/species_distribution_density.jpg"), format = "file"), 
#   tar_target(fig_native_blocks, plot_native_distribution_blocks("fig/native_distribution_blocks.jpg"), 
#               format = "file")
 NULL
)
