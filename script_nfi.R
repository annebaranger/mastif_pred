# _targets.R file

#library
library(targets)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("stringr","ggplot2","tidyr","dplyr","terra","factoextra","modi","tibble"),
               error = "continue") 

mastif.eu=tar_read(mastif.eu,store = "target_data")
mastif.am=tar_read(mastif.am,store = "target_data")
meanClimate_species=tar_read(meanClimate_species,store="target_gbif")

#Targets
list(
  # species selection
  tar_target(
    species_selection,
    select_species(mastif.am,
                   mastif.eu,
                   meanClimate_species)
  ),
  tar_target(
    sp.select.am,
    species_selection$sp.select.am
  ),
  tar_target(
    sp.select.eu,
    species_selection$sp.select.eu
  ),
  
  # nfi data
  tar_target(
    clim_list,
    get_climate()
  ),
  tar_target(
    df.nfi.eu,
    get_nfi(clim_list,
            continent="europe")
  ),
  tar_target(
    df.nfi.am,
    get_nfi(clim_list,
            continent="america")
  ),
  
  # margins
  tar_target(
    nfi.am_clim,
    compute_acp(df.nfi.am,
                   clim_var=c("sgdd","wai"))
  ),
  tar_target(
    nfi.eu_clim,
    compute_acp(df.nfi.eu,
                   clim_var=c("sgdd","wai"))
    ),
  
  # fecundity
  tar_target(
    fec.eu,
    get_fecundity(continent="europe",
                  sp.select.eu)
  ),
  tar_target(
    fec.am,
    get_fecundity(continent="america",
                  sp.select.am)
  ),
  
  # fecundity margin
  tar_target(
    sp.margin.eu,
    compute_margin(nfi.eu_clim$df.nfi.clim,
                   fec.eu)
  ),
  tar_target(
    sp.margin.am,
    compute_margin(nfi.am_clim$df.nfi.clim,
                   fec.am)
  ),
  
  # fecundity x plot
  tar_target(
    fecundity.eu_clim,
    fec.eu |> 
      left_join(nfi.eu_clim$df.nfi.clim)|> 
      left_join(sp.margin.eu) |> 
      mutate(margin=case_when(PC1<quant05~"low",
                              PC1>quant475&PC1<quant525~"opt",
                              PC1>quant95~"up",
                              TRUE~NA))
  ),
  tar_target(
    fecundity.am_clim,
    fec.am |> 
      left_join(nfi.am_clim$df.nfi.clim)|> 
      left_join(sp.margin.am) |> 
      mutate(margin=case_when(PC1<quant05~"low",
                              PC1>quant475&PC1<quant525~"opt",
                              PC1>quant95~"up",
                              TRUE~NA))
  ),
  NULL
)