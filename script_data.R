# _targets.R file

#library
library(targets)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","tidyr","dplyr","terra","factoextra","taxize"),
               error = "continue") 

#Targets
list(
  tar_target(
    mastif.eu,
    get_mastif(dir.data="data/mastif_plots/fittedFecundityMastifEurope.rdata")
  ),
  tar_target(
    mastif.am,
    get_mastif(dir.data="data/mastif_plots/fittedFecundityMastifNorthAmerica.rdata")
  ),
  tar_target(
    species.eu,
    get_specieslist(mastif.eu$df.species.select$species,
                    block="europe")
  ),
  tar_target(
    species.am,
    get_specieslist(mastif.am$df.species.select$species,
                    block="america")
  ),
  tar_target(
    species.phylo,
    rbind(species.am,species.eu)
  ),
  NULL
)