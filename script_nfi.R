# _targets.R file

#library
library(targets)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","tidyr","dplyr","terra","factoextra","modi","tibble","mastif"),
               error = "continue") 

mastif.eu=tar_read(mastif.eu,store = "target_data_2")
mastif.am=tar_read(mastif.am,store = "target_data_2")
# meanClimate_species=tar_read(meanClimate_species,store="target_gbif")

#Targets
list(
  # nfi data
  tar_target(
    clim_list,
    get_climate()
  ),
  tar_target(
    df.nfi.eu,
    get_nfi(clim_list,
            fit="fit2024",
            continent="europe")
  ),
  tar_target(
    df.nfi.am,
    get_nfi(clim_list,
            fit="fit2024",
            continent="america")
  ),
  
  # # margins
  # tar_target(
  #   nfi.am_clim,
  #   compute_acp(df.nfi.am,
  #                  clim_var=c("sgdd","wai"))
  # ),
  # tar_target(
  #   nfi.eu_clim,
  #   compute_acp(df.nfi.eu,
  #                  clim_var=c("sgdd","wai"))
  #   ),
  
  # fecundity
  tar_target(
    fec.eu,
    get_fecundity(continent="europe",
                  fit="fit2024",
                  mastif.eu$df.species.select$species)
  ),
  tar_target(
    fec.am,
    get_fecundity(continent="america",
                  fit="fit2024",
                  mastif.am$df.species.select$species)
  ),
  
  # fecundity cv
  # tar_target(
  #   cv.eu,
  #   get_cv(continent="europe",
  #          sp.select.eu)
  # ),
  # tar_target(
  #   cv.am,
  #   get_cv(continent="america",
  #          sp.select.am)
  # ),
  # 
  
  # fecundity x plot
  tar_target(
    fecundity.eu_clim,
    fec.eu |> 
      # left_join(cv.eu,by=c("plot","species")) |> 
      left_join(df.nfi.eu)
  ),
  tar_target(
    fecundity.am_clim,
    fec.am |> 
      # left_join(cv.am,by=c("plot","species")) |> 
      left_join(df.nfi.am)
  ),
  
  # species selection #
  tar_target(
    species_selection,
    select_species(fecundity.eu_clim,
                   fecundity.am_clim,
                   mastif.am,
                   mastif.eu)),
  # tar_target(
  #   sp.select.am,
  #   species_selection$df.select[species_selection$df.select$block=="america"&
  #                                 species_selection$df.select$select.quant==TRUE,"species"][[1]] # extract only selection from am
  # ),
  # tar_target(
  #   sp.select.eu,
  #   species_selection$df.select[species_selection$df.select$block=="europe"&
  #                                 species_selection$df.select$select.quant==TRUE,"species"][[1]]
  # ),
  ## create a dataframe with species selection and phylogeny
  tar_target(
    species.phylo.file,
    "target_data/objects/species.phylo",
    format="file"
  ),
  tar_target(
    species.phylo, # phylogeny of Mastif species
    readRDS(species.phylo.file)
  ),
  tar_target(
    species_class, #classification of nfi species
    class_species(fecundity.am_clim,
                  fecundity.eu_clim)
  ),
  tar_target(
    phylo.select,
    species.phylo |>
      inner_join(species_class) |>  # keep only species present in nfi and mastif
      filter(species %in% species_selection$select.quant) # filter selection
      
  ),
  NULL
)
