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


list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - BLUR DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # was done using complete dataset
  # tar_target(mastif.am,
  #            blur_data(mastif.file="data/not_public/mastif.am")),
  # tar_target(mastif.eu,
  #            blur_data(mastif.file="data/not_public/mastif.eu")),
  # tar_target(fecundity.am_clim,
  #            blur_nfi(nfi.file="data/not_public/fecundity.am_clim")),
  # tar_target(fecundity.eu_clim,
  #            blur_nfi(nfi.file="data/not_public/fecundity.eu_clim")),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - GBIF SPECIES DATA ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(gbif_niche,
             read.table("data/sp_gbif_climate.csv",header = TRUE)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - SPECIES SELECTION ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    species.phylo, # phylogeny of Mastif species
    readRDS("data/species.phylo")
  ),
  tar_target(
    species_class, #classification of nfi species
    readRDS("data/species_class")
  ),
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
    phylo.select,
    species.phylo |> # species in mastif
      inner_join(species_class) |>  # keep only species present in nfi and mastif
      left_join(species_selection_quant) |>
      filter(!is.na(cold_valid)) |> 
      separate(species_l,into=c("genus","species_only"),remove=FALSE) |> 
      mutate(s_p=str_replace(species_l," ","_"),
             sp_little=paste0(tolower(substr(genus,1,4)),substr(species_only,1,4)),
             dir.sp=case_when(block=="america"~paste0("data/USTreeAtlas-main/shp/",sp_little),
                              block=="europe"~paste0("data/chorological_maps_dataset/",species_l,"/shapefiles")),
             direx=dir.exists(dir.sp)) |> 
      filter(direx) |> 
      mutate(
        file.sp=case_when(block=="america"~sp_little,
                          block=="europe"~paste0(s_p,"_plg")),
        file.sp=case_when(species=="fagusSylvatic"~"Fagus_sylvatica_sylvatica_plg",
                          species=="quercusIlex"~"Quercus_ilex_ilex_plg",
                          TRUE~file.sp),
        filex=file.exists(file.path(dir.sp,paste0(file.sp,".shp"))) 
      ) |>
      filter(filex)),
  tar_target(
    species_selection,
    phylo.select$species
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - PREPARE DATA FOR FIT---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    phylo.zone, # correction of zone attribution
    phylo.select|> 
      mutate(zone=case_when(species=="abiesCephalon"~"europe",
                            species=="acerCircinat"~"west",
                            species=="acerSpicatum"~"east",
                            species=="amelanchArborea"~"east",
                            species=="asiminaTriloba"~"east",
                            species=="cedrusAtlantic"~"europe",
                            species=="frangulaAlnus"~"europe",
                            species=="linderaBenzoin"~"east",
                            TRUE~zone)) 
  ),
  tar_target(
    fecundity.fit.1, # margin set with 10th percentiles
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    phylo.zone,
                    species_selection,
                    thresh=0.1)
  ),
  tar_target(
    fecundity.fit.25, # margin set with 25th percentiles
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    phylo.zone,
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
