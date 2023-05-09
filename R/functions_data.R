#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analyses.R  
#' @description R script containing all functions relative to data
#               extraction and formating
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Data collection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract mastif data over europe and extract climate
#' @param dir.data directory of mastif rdata
#' @param files.climate list of files of climate variables to extract
#' "mat","tmin","map"
get_mastif_europe <- function(dir.data="data/fecundityMastif.rdata",
                              files.climate=list(mat="data/CHELSA/CHELSA_bio10_1.tif",
                                                 tmin="data/CHELSA/CHELSA_bio10_6.tif",
                                                 map="data/CHELSA/CHELSA_bio10_12.tif")){
  clim=rast(lapply(names(files.climate),
                   function(z){r=rast(files.climate[[z]])
                               names(r)=z
                               return(r)}
                   )
            )
  load(dir.data)
  fecundityPred<-fecundityPred |> 
    tibble::rownames_to_column(var="ID") |> 
    separate(ID,into=c("treeID","year"),sep="_",remove=FALSE) |> 
    separate(treeID,into=c("plotID",NA),sep="-",remove=FALSE)
  
  fecundityPred<-cbind(fecundityPred,
                       terra::extract(clim,
                               y=data.frame(x=fecundityPred$lon,
                                            y=fecundityPred$lat))[,-1]
  )

  #build df with only plots
  df.plot<-fecundityPred |> 
    select(plotID,lon,lat,year) |> 
    distinct() |> 
    group_by(plotID,lon,lat) |> 
    summarise(year_serie=n()) |> 
    ungroup()
  
  # build a df with species selection, by restricting only to species with more 
  # than 200 individuals
  df.species <- fecundityPred |> 
    select(treeID,species,lon,lat) |> 
    distinct() |> 
    group_by(species) |> 
    filter(n()>200) |> 
    summarise(n=n(),
              minlat=min(lat),
              maxlat=max(lat),
              minlon=min(lon),
              maxlon=max(lon)) |> 
    ungroup()
  
  #select trees
  df.tree<-fecundityPred |> 
    filter(species %in% df.species$species)
  return(list(df.alltree=fecundityPred,
              df.allplot=df.plot,
              df.species.select=df.species,
              df.tree=df.tree
              ))
}