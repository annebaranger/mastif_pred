# _targets.R file

#library
library(targets)
# lapply(c("stringr","ggplot2","tidyr","dplyr"),require,character.only=TRUE)


#Options
source("R/functions_data.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("rworldmap","stringr","ggplot2","tidyr","dplyr","terra","factoextra","modi","tibble","mastif"),
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
  # tar_target(
  #   df.nfi.eu,
  #   get_nfi(clim_list,
  #           fit="fit2024",
  #           continent="europe")
  # ),
  # tar_target(
  #   df.nfi.am,
  #   get_nfi(clim_list,
  #           fit="fit2024",
  #           continent="america")
  # ),
  # 
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
  # tar_target(
  #   fec.eu,
  #   get_fecundity(continent="europe",
  #                 fit="fit2024",
  #                 mastif.eu$df.species.select$species)
  # ),
  # tar_target(
  #   fec.am,
  #   get_fecundity(continent="america",
  #                 fit="fit2024",
  #                 mastif.am$df.species.select$species)
  # ),
  
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
  # tar_target(
  #   fecundity.eu_clim,
  #   fec.eu |> 
  #     # left_join(cv.eu,by=c("plot","species")) |> 
  #     left_join(df.nfi.eu)
  # ),
  # tar_target(
  #   fecundity.am_clim,
  #   fec.am |> 
  #     # left_join(cv.am,by=c("plot","species")) |> 
  #     left_join(df.nfi.am)
  # ),
  # 
  # get nfi new function
  tar_target(
    fecundity.am_clim,
    get_nfipred_plot(continent="america",
                     fit="fit2023",
                     chelsa_files=list(mat="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                                       map="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                                       pet="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif"),
                     sp.select=mastif.am$df.species.select$species)
  ),
  tar_target(
    fecundity.eu_clim,
    get_nfipred_plot(continent="europe",
                     fit="fit2023",
                     chelsa_files=list(mat="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                                       map="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                                       pet="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif"),
                     sp.select=mastif.eu$df.species.select$species)
  ),
  # species selection #
  tar_target(
    species_selection_nfi,
    select_species_nfi(fecundity.eu_clim,
                       fecundity.am_clim,
                       mastif.am,
                       mastif.eu)),
  tar_target(gbif_niche_file,
             "output/sp_gbif_climate.csv",
             format = "file"),
  tar_target(gbif_niche,
             read.table(gbif_niche_file,header = TRUE)),
  tar_target(
    species_selection_gbif,
    select_species_gbif(gbif_niche,
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
    "target_data_2/objects/species.phylo",
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
    species_selection_narrow,
    intersect(species_selection_gbif$select.rank,species_selection_nfi$select.rank)
  ),
  tar_target(
    species_selection_medium,
    intersect(species_selection_gbif$select.quant,species_selection_nfi$select.quant)
  ),
  tar_target(
    species_selection_large,
    unique(c(species_selection_gbif$select.quant,species_selection_nfi$select.quant))
  ),
  # tar_target(
  #   phylo.select,
  #   species.phylo |>
  #     inner_join(species_class) |>  # keep only species present in nfi and mastif
  #     filter(species %in% species_selection_narrow) # filter selection
  #     
  # ),
  tar_target(
    phylo.select,
    species.phylo |>
      inner_join(species_class) |>  # keep only species present in nfi and mastif
      filter(species %in% species_selection_narrow) |> 
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
      filter(filex)
    
  ),
  tar_target(
    selection_expertmaps,
    select_species_expert(map_file="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                          pet_file="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif",
                          mat_file="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                          species_selection_gbif,
                          species_selection_nfi,
                          phylo.select
    )
  ),
  NULL
)
