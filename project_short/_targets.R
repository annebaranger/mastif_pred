# _targets.R file
library(targets)
# lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
source("R/functions_data.R")
source("R/functions_fit.R")
packages.in <- c("rworldmap","stringr","ggplot2","tidyr","dplyr","terra","factoextra","modi","tibble","mastif")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
options(tidyverse.quiet = TRUE)
sf::sf_use_s2(FALSE)
tar_option_set(packages = packages.in)
mastif.am=readRDS("_targets/objects/mastif.am")
mastif.eu=readRDS("_targets/objects/mastif.eu")
fecundity.am_clim=readRDS("_targets/objects/fecundity.am_clim")
fecundity.eu_clim=readRDS("_targets/objects/fecundity.eu_clim")
species.meta=readRDS("_targets/objects/species.meta")
list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - BLUR DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # was done using complete dataset
  # tar_target(mastif.am,
  #            blur_data(mastif.file="data/not_public/mastif.am",
  #                      species.meta.file="data/not_public/species.meta.csv")),
  # tar_target(mastif.eu,
  #            blur_data(mastif.file="data/not_public/mastif.eu",
  #                      species.meta.file="data/not_public/species.meta.csv")),
  # tar_target(fecundity.am_clim,
  #            blur_nfi(nfi.file="data/not_public/fecundity.am_clim",
  #                     species.meta.file="data/not_public/species.meta.csv")),
  # tar_target(fecundity.eu_clim,
  #            blur_nfi(nfi.file="data/not_public/fecundity.eu_clim",
  #                     species.meta.file="data/not_public/species.meta.csv")),
  # 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - GBIF SPECIES DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(gbif_niche,
             read.csv("data/sp_gbif_climate.csv",header = TRUE)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - SPECIES SELECTION ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # tar_target(
  #   species.meta, # phylogeny of Mastif species
  #   read.csv("data/not_public/species.meta.csv") |> select(-species_name,-X)
  # ),
  tar_target(
    species_selection_nfi_margin,
    select_species_nfi_margin(fecundity.eu_clim,
                              fecundity.am_clim,
                              mastif.am,
                              mastif.eu)),
  tar_target(
    species_selection_gbif_margin,
    select_species_gbif_margin(gbif_niche,
                               mastif.am,
                               mastif.eu)),
  tar_target(
    species_selection_rank,
    species_selection_nfi_margin$species_rank |> 
      filter(species%in%species_selection_gbif_margin$species_rank$species)
  ),
  tar_target(
    species_selection_quant,
    species_selection_nfi_margin$species_quant |> 
      filter(species%in%species_selection_gbif_margin$species_quant$species)
  ),
  tar_target(
    species.select,
    species.meta |>  # keep only species present in nfi and mastif
      left_join(species_selection_quant) |>
      filter(!is.na(cold_valid)) |> 
      filter(!species %in% c(15,39,70))),
  tar_target(
    species_selection,
    species.select$species
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - PREPARE DATA FOR FIT---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    fecundity.fit.1, # margin set with 10th percentiles
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    species.select,
                    species_selection,
                    thresh=0.1)
  ),
  tar_target(
    fecundity.fit.25, # margin set with 25th percentiles
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    species.select,
                    species_selection,
                    thresh=0.25)
  ),
  tar_target(
    species.biome,
    get_biome(fecundity.fit.1)
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - FIT MODEL ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    fitcontinentdiscrete_25,
    fit.continent.discrete(fecundity.fit.25,
                           folder="model/continent_discrete_25_margin")
  ),
   tar_target(
    fitspecies_25,
    fit.species(fecundity.fit.25,
                folder="model/species_continuous_25_margin")
  ),
  tar_target(
    fitbiomecontinuous_25,
    fit.biome.continuous(subset(fecundity.fit.25,dh_valid),
                         species.biome,
                         folder="model/biome_continous_25_margin")
  ),
  tar_target(
    fitbiomediscrete_25,
    fit.biome.discrete(fecundity.fit.25,
                       species.biome,
                       folder="model/biome_discrete_25_margin")
  ),
  NULL
)
