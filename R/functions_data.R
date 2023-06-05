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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - MASTIF data collection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract mastif data over europe and extract climate
#' @param dir.data directory of mastif rdata
#' @param files.climate list of files of climate variables to extract
#' "mat","tmin","map"
get_mastif <- function(dir.data="data/fecundityMastif.rdata",
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


get_specieslist <- function(species.list,block){
  species.world=read.csv2("data/global_tree_search_trees_1_7.csv")[,1:2] |> 
    separate(TaxonName,into=c("genus","sp"),remove=FALSE) |> 
    mutate(species=paste0(substr(tolower(genus),1,8),substr(str_to_title(sp),1,8))) 
  species.code=species.world |> 
    filter(species%in%species.list) |> 
    mutate(block=block)
  return(species.code)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2- Species selection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select_species <- function(mastif.am,
                           mastif.eu,
                           meanClimate_species){
  df.tree=rbind(mastif.am$df.tree,mastif.eu$df.tree) 
  df.species=rbind(mastif.am$df.species.select,mastif.eu$df.species.select)
  
  # format gbif output
  df.niche= meanClimate_species |> 
    separate(species,into=c("genus","sp"),remove=FALSE) |> 
    mutate(species_l=species,
           species=paste0(substr(tolower(genus),1,8),substr(str_to_title(sp),1,8))) |> 
    rename_with(.cols=c("mat","map","tmin"),
                ~paste0(.x,".opt")) |> 
    select(-genus,-sp) |> 
    mutate(mat.range=mat.high-mat.low,
           map.range=map.high-map.low,
           tmin.range=tmin.high-tmin.low) |> 
    pivot_longer(cols=-c("species","species_l"),
                 names_to = "clim",
                 values_to = "clim_niche")
  
  df.range<-df.niche |> 
    select(species,clim,clim_niche) |> 
    filter(grepl("range",clim)) |> 
    separate(clim,into=c("var","drop")) |> 
    select(-drop) |> 
    rename(range="clim_niche")
  
  # compute weight for trees
  df.tree.weight=df.tree |>
    # replace 0 estimations by NA and compute, rescale mat and tmin
    mutate(fecEstSe=na_if(fecEstSe,0),
           fecEstMu=na_if(fecEstMu,0),
           across(c("mat","tmin"),
                  ~.x/10)
    ) |> 
    filter(!is.na(fecEstMu)) |> filter(!is.na(fecEstSe)) |> 
    # compute numbers of observations per plot
    group_by(plotID,species,mat,map,tmin) |> 
    summarise(n_plot=n(),
              fecEstMu_plot=mean(fecEstMu,na.rm=TRUE),
              se_plot=sqrt(sum(fecEstSe^2)/n_plot),
              cv_plot=se_plot/fecEstMu_plot) |> 
    # Compute mean of se per plot
    ungroup() |> 
    pivot_longer(cols=c("mat","tmin","map"),
                 names_to = "clim",
                 values_to = "clim_val")  |> 
    filter(!is.na(clim_val)) |> 
    group_by(species,clim) |>
    mutate(opt=weighted.quantile(clim_val,1/cv_plot,prob=0.5)[1],
           low=weighted.quantile(clim_val,1/cv_plot,prob=0.025)[1],#,probs=0.1)[1],
           high=weighted.quantile(clim_val,1/cv_plot,prob=0.975)[1]) |> 
    ungroup()
  
  (df.tree.weight |>
    # format to match df.niche 
    select(species,clim,opt,low,high) |> unique() |> 
    pivot_wider(names_from = "clim",
                values_from = c("opt","high","low"),
                names_sep = ".") |> 
    pivot_longer(cols=-species,
                 names_to = "clim",
                 values_to = "clim_obs") |> 
    mutate(clim=paste0(str_split_i(clim,"\\.",2),".",str_split_i(clim,"\\.",1))) |> 
    #merge with gbif climate
    left_join(df.niche,by=c("species","clim")) |> 
    relocate(species_l,.after=species)|> 
    separate(clim,into=c("var","quant"),remove=FALSE) |> 
    relocate(var,.after=species_l) |> 
  # merge climate range
    left_join(df.range,by=c("species","var")) |> 
    filter(!is.na(species_l)) |> 
    mutate(dif=case_when(grepl("opt",clim)~100*abs(clim_obs-clim_niche)/range,
                         grepl("high",clim)~100*(clim_obs-clim_niche)/range,
                         grepl("low",clim)~100*(clim_niche-clim_obs)/range)) |> 
    filter(dif<(-50)) |> 
    select(species) |> unique())$species->sp.deleted
  sp.select=setdiff(unique(df.tree$species),sp.deleted)
  sp.select.eu=sp.select[sp.select %in% mastif.eu$df.species.select$species]
  sp.select.am=sp.select[sp.select %in% mastif.am$df.species.select$species]
  
  return(list(df.gbif=df.niche,
              df.tree_w=df.tree.weight,
              sp.select.am=sp.select.am,
              sp.select.eu=sp.select.eu))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - NFI data collection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get a list of climate files
#' @param clim_list list of climatic variables of interest
get_climate<-function(clim_list=list(mat="data/CHELSA/CHELSA_bio10_1.tif",
                                     tmin="data/CHELSA/CHELSA_bio10_6.tif",
                                     map="data/CHELSA/CHELSA_bio10_12.tif",
                                     sgdd="data/CHELSA/CHELSA_gdd5_1981-2010_V.2.1.tif",
                                     pet="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif")){
  return(clim_list)
}

#' Load NFI plots for a targetted continent with associated climatic variables
#' and plot size
#' @param clim_list climatic variable list
#' @param continent wether "europe" or "america"
get_nfi <- function(clim_list,
                    continent){
  # load clim
  clim=rast(lapply(names(clim_list),
                   function(z){
                     r=rast(clim_list[[z]])
                     names(r)=z
                     return(r)
                   }
  )
  )
  #load plotData.rdata
  nfi.file=paste0("data/",continent,"/plotData.rdata")
  size.file=paste0("data/",continent,"/SIZE.rdata")
  load(nfi.file)
  load(size.file)
  SIZE=as.data.frame(SIZE) |> 
    rownames_to_column(var="plot")
  plotData<-plotData |> 
    left_join(SIZE,by="plot")
  df.nfi<-cbind(plotData,
                extract(clim,
                        y=data.frame(x=plotData$lon,
                                     y=plotData$lat))[,-1]) |> 
    mutate(wai=(map-12*pet)/(12*pet))
  
  return(df.nfi)
}


#' Compute margin class for each plot of nfi
#' @param df.nfi
compute_margin<-function(df.nfi,
                         clim_var=c("sgdd","wai")){
  clim_pca<-prcomp(df.nfi[,clim_var] |> drop_na(),
                   scale. = TRUE) 
  clim_pred <- cbind(df.nfi,
                     predict(clim_pca, newdata = df.nfi[,clim_var])) |> 
    rownames_to_column() |> 
    select(-rowname)
  
  return(list(clim_pca=clim_pca,
              df.nfi.clim=clim_pred))
}


#' Load fecundity for each plot, including mean and sd estimation, ISP a,d basal area
#' @param continent
#' @param sp.select
get_fecundity <- function(continent,
                          sp.select){
  # load BA
  load(paste0("data/",continent,"/BA.rdata"))
  BA<-as.data.frame(BA)|> 
    rownames_to_column(var="plot") |> 
    select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "BA")
  
  #load ISP
  load(paste0("data/",continent,"/ISP.rdata"))
  ISP<-as.data.frame(ISP)|> 
    rownames_to_column(var="plot") |> 
    select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "ISP")
  
  #load fecgmmu&sd
  load(paste0("data/",continent,"/fecGmSd.rdata"))
  fecGmSd<-as.data.frame(summary)|> 
    rownames_to_column(var="plot") |> 
    select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmSd")
  load(paste0("data/",continent,"/fecGmMu.rdata"))
  fecGmMu<-as.data.frame(summary)|> 
    rownames_to_column(var="plot") |> 
    select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmMu")
  rm(summary)
  
  fecundity<-BA |> 
    left_join(ISP,by=c("plot","species")) |> 
    left_join(fecGmMu,by=c("plot","species")) |> 
    left_join(fecGmSd,by=c("plot","species")) |> 
    mutate(fecGmMu=na_if(fecGmMu,0),
           fecGmSd=na_if(fecGmSd,0),
           ISP=na_if(ISP,0)) |> 
    group_by(plot,species) |> 
    mutate(n_plot=n(),
           fecGmMu_plot=mean(fecGmMu,na.rm=TRUE),
           se_plot=sqrt(sum(fecGmSd^2)/n_plot),
           cv_plot=se_plot/fecGmMu_plot) |> 
    # Compute mean of se per plot
    ungroup()
  
  return(fecundity)
}