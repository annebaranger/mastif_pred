# _targets.R file

#library
library(targets)
library(dplyr)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_fit.R")
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = c("stringr","ggplot2","tidyr","dplyr","terra",
                            "modi","tibble","rstan",
                            "future"),
               error = "continue",
               memory = "transient") 
future::plan(future::multisession, workers = 6)
fecundity.eu_clim=tar_read(fecundity.eu_clim,store="target_nfi")
fecundity.am_clim=tar_read(fecundity.am_clim,store="target_nfi")
# species_selection_narrow=tar_read(species_selection_narrow,store="target_nfi")
phylo.select=tar_read(phylo.select_margin,store="target_nfi") #|>
# filter(species%in%species_selection_narrow)
species_selection=phylo.select$species

#Targets
list(
  tar_target(
    phylo.zone,
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
    fecundity.fit.1,
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    phylo.zone,
                    species_selection,
                    thresh=0.1)
  ),
  tar_target(
    fecundity.fit.25,
    raw_data_margin(fecundity.eu_clim,
                    fecundity.am_clim,
                    phylo.zone,
                    species_selection,
                    thresh=0.25)
  ),
  # tar_target(
  #   fecundity.fit.05,
  #   raw_data(fecundity.eu_clim,
  #            fecundity.am_clim,
  #            species_selection$select.quant,
  #            phylo.zone,
  #            thresh=0.05)
  # ),
  tar_target(
    species.biome,
    get_biome(fecundity.fit.1)
  ),
  # tar_target(
  #   fitcontinentdiscrete,
  #   fit.continent.discrete(fecundity.fit.1,
  #                          folder="model/continent_discrete")
  # ),
  # tar_target(
  #   fitcontinentdiscreteexcl,
  #   fit.continent.discrete.excl(fecundity.fit.1,
  #                               folder="model/continent_discrete")
  # ),
  # tar_target(
  #   fitcontinentcontinous,
  #   fit.continent.continous(fecundity.fit.1,
  #                           excluded=TRUE,
  #                           folder="model/continent_continous")
  # ),
  # tar_target(
  #   fitspecies,
  #   fit.species(fecundity.fit.1,
  #               folder="model/species_continuous")
  # ),
  # tar_target(
  #   fitbiomecontinuous,
  #   fit.biome.continuous(fecundity.fit.1,
  #                        species.biome,
  #                        folder="model/biome_continous")
  # ),
  # tar_target(
  #   fitbiomediscrete,
  #   fit.biome.discrete(fecundity.fit.1,
  #                      species.biome,
  #                      folder="model/biome_discrete")
  # ),
  ## margin 0.25
  tar_target(
    fitcontinentdiscrete_25,
    fit.continent.discrete(fecundity.fit.25,
                           folder="model/continent_discrete_25_margin")
  ),
  tar_target(
    fitcontinentdiscreteexcl_25,
    fit.continent.discrete.excl(fecundity.fit.25,
                                folder="model/continent_discrete_25_margin")
  ),
  tar_target(
    fitcontinentcontinous_25,
    fit.continent.continous(subset(fecundity.fit.25,dh_valid),
                            excluded=TRUE,
                            folder="model/continent_continous_25_margin")
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
    fitbiomediscrete_dh_25,
    fit.biome.discrete.dh(subset(fecundity.fit.25,dh_valid),
                         species.biome,
                         folder="model/biome_continous_dh_25_margin")
  ),
  tar_target(
    fitbiomediscrete_25,
    fit.biome.discrete(fecundity.fit.25,
                       species.biome,
                       folder="model/biome_discrete_25_margin")
  ),
  NULL
)
