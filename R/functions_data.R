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


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - MASTIF data collection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract mastif data over europe and extract climate
#' @param dir.data directory of mastif rdata
#' @param files.climate list of files of climate variables to extract
#' "mat","tmin","map"
get_mastif <- function(dir.data=dir.data,
                       # c("data/mastif_plots/fittedFecundityMastifNorthAmerica.rdata",
                       #            "data/mastif_plots/fittedFecundityMastifEurope.rdata"),
                       files.climate=list(mat="data/CHELSA/CHELSA_bio10_1.tif",
                                          tmin="data/CHELSA/CHELSA_bio10_6.tif",
                                          map="data/CHELSA/CHELSA_bio10_12.tif")){
  clim=rast(lapply(names(files.climate),
                   function(z){r=rast(files.climate[[z]])
                               names(r)=z
                               return(r)}
                   )
            )
  # fecundityPred_i=data.frame()
  # for (i in 1:length(dir.data)){
  #   print(i)
  #   load(dir.data[i])
  #   fecundityPred_i=bind_rows(fecundityPred_i,
  #                             fecundityPred)
  # }
  fecundityPred=loadRData(dir.data)
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
    dplyr::select(treeID,species,lon,lat) |> 
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
  # load a list of world species
  species.world=read.csv2("data/global_tree_search_trees_1_7.csv")[,1:2] |> 
    separate(TaxonName,into=c("genus","sp"),remove=FALSE) |> 
    mutate(species=paste0(substr(tolower(genus),1,8),substr(str_to_title(sp),1,8))) # recreate ID of mastif
  # create correspondence between ID and full species name
  species.code=species.world |> 
    filter(species%in%species.list) 
  
  # add phylogeny
  classification(unique(species.code$TaxonName),db="gbif",rows=1)->out # get phylogeny from species names
  class2tree(out,check=FALSE)->out2 # random command
  
  species.phylo <- out2$classification |> # extract phylogeny from taxize object
    tibble::rownames_to_column( var = "species_init") |>
    # reformat long species name
    mutate(species_l=case_when(is.na(species)~species_init,
                             TRUE~species)) |> 
    # reformat taxa
    mutate(taxa=case_when(class=="Pinopsida"~"gymnosperm",
                          class=="Magnoliopsida"~"angiosperm")) |>
    dplyr::select(species_l,genus,family,order,taxa) |> # select relevant fields
    left_join(species.code[,c("TaxonName","species")],by=c("species_l"="TaxonName")) |>
    unique()
  all.species<-data.frame(species=species.list) |> 
    left_join(species.phylo)|> 
    mutate(block=block)
  
  rm(out,out2)
  
  return(all.species)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2- Species selection ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Select species based on NFI
#' @description Compare climatic distribution of NFI plots and MASTIF plots. 
#' Select species for which only MASTIF coverage is large enough compared to NFI.
#' @param fecundity.eu_clim NFI plots for European species
#' @param fecundity.am_clim NFI plots for American species
#' @param mastif.eu MASTIF data for European species
#' @param mastif.am MASTIF data for American species
select_species_nfi<-function(fecundity.eu_clim,
                             fecundity.am_clim,
                             mastif.eu,
                             mastif.am){
  
  mastif.species.niche<-rbind(mastif.am$df.tree, 
                              mastif.eu$df.tree) |> 
    mutate(dh=pet-(map/12)) |> 
    mutate(fecEstSe=na_if(fecEstSe,0),
           fecEstMu=na_if(fecEstMu,0)) |> 
    filter(!is.na(fecEstMu)) |> filter(!is.na(fecEstSe)) |> 
    # compute numbers of observations per plot
    group_by(plotID,species,mat,dh) |> 
    summarise(n_plot=n(),
              fecEstMu_plot=mean(fecEstMu,na.rm=TRUE),
              se_plot=sqrt(sum(fecEstSe^2)/n_plot),
              cv_plot=se_plot/fecEstMu_plot) |> 
    # Compute mean of se per plot
    ungroup() |> 
    pivot_longer(cols=c("mat","dh"),
                 names_to = "clim",
                 values_to = "clim_val")  |> 
    filter(!is.na(clim_val)) |> 
    group_by(species,clim) |>
    summarise(#opt=weighted.quantile(clim_val,1/cv_plot,prob=0.5)[1],
      qlow=weighted.quantile(clim_val,1/cv_plot,prob=0.025)[1],#,probs=0.1)[1],
      qhigh=weighted.quantile(clim_val,1/cv_plot,prob=0.975)[1],
      qrange.mastif=qhigh-qlow,
      rhigh=weighted.quantile(clim_val,1/cv_plot, prob = 1- min(5, n())/n())[1],
      rlow=weighted.quantile(clim_val,1/cv_plot, prob = min(5, n())/n())[1],
      rrange.mastif=rhigh-rlow) |> 
    ungroup() |> 
    pivot_longer(cols=c("qlow","qhigh","rhigh","rlow"),
                 names_to = "limit",
                 values_to = "value_mastif") |> 
    mutate(range.mastif=case_when(grepl("q",limit)~qrange.mastif,
                                  grepl("r",limit)~rrange.mastif)) |> 
    select(-qrange.mastif,-rrange.mastif)
  
  nfi.species.niche<-rbind(fecundity.eu_clim,
                           fecundity.am_clim) |> 
    # mutate(dh=12*pet-map) |> 
    filter(BA!=0) |> 
    filter(dh>(-2500)) |>
    pivot_longer(cols=c("mat","dh"),
                 names_to = "clim",
                 values_to = "clim_val")  |> 
    filter(!is.na(clim_val)) |> 
    group_by(species,clim) |>
    summarise(#opt=weighted.quantile(clim_val,BA,prob=0.5)[1],
      qlow=weighted.quantile(clim_val,BA,prob=0.025)[1],#,probs=0.1)[1],
      qhigh=weighted.quantile(clim_val,BA,prob=0.975)[1],
      qrange.nfi=qhigh-qlow,
      rhigh=weighted.quantile(clim_val,BA, prob = 1- min(5, n())/n())[1],
      rlow=weighted.quantile(clim_val,BA, prob = min(5, n())/n())[1],
      rrange.nfi=rhigh-rlow) |> 
    ungroup() |> 
    pivot_longer(cols=c("qlow","qhigh","rhigh","rlow"),
                 names_to = "limit",
                 values_to = "value_nfi") |> 
    mutate(range.nfi=case_when(grepl("q",limit)~qrange.nfi,
                               grepl("r",limit)~rrange.nfi)) |> 
    select(-qrange.nfi,-rrange.nfi)
  
  
  mastif.species.niche |> 
    left_join(nfi.species.niche,by=c("species","clim","limit")) |> 
    filter(!is.na(range.nfi)) |> 
    mutate(dev=case_when(grepl("high",limit)~100*(value_mastif-value_nfi)/range.nfi,
                         grepl("low",limit)~100*(value_nfi-value_mastif)/range.nfi)) |> 
    filter(grepl("r",limit)) |> 
    group_by(species) |> 
    mutate(ind_neg=sum(dev>(-50))) |> ungroup() |> 
    filter(ind_neg==4) |> 
    select(species) |> unique() |> pull() -> select.rank
  mastif.species.niche |> 
    left_join(nfi.species.niche,by=c("species","clim","limit")) |> 
    filter(!is.na(range.nfi)) |> 
    mutate(dev=case_when(grepl("high",limit)~100*(value_mastif-value_nfi)/range.nfi,
                         grepl("low",limit)~100*(value_nfi-value_mastif)/range.nfi)) |> 
    filter(grepl("q",limit)) |> 
    group_by(species) |> 
    mutate(ind_neg=sum(dev>(-50))) |> ungroup() |> 
    filter(ind_neg==4) |> 
    select(species) |> unique() |> pull() -> select.quant
  
  return(list(df.quant=nfi.species.niche,
              select.rank=select.rank,
              select.quant=select.quant))
}




#' Select species based on GBIF
#' @description Compare climatic distribution of GBIF plots and MASTIF plots. 
#' Select species for which only MASTIF coverage is large enough compared to GBIF
#' @note  NFI plots do not always cover the whole distribution of species 
#' @param fecundity.eu_clim NFI plots for European species
#' @param fecundity.am_clim NFI plots for American species
#' @param mastif.eu MASTIF data for European species
#' @param mastif.am MASTIF data for American species
select_species_gbif<-function(gbif_niche,
                              mastif.am,
                              mastif.eu){
  
  mastif.species.niche<-rbind(mastif.am$df.tree, 
                              mastif.eu$df.tree) |> 
    mutate(dh=pet-(map/12)) |> 
    mutate(fecEstSe=na_if(fecEstSe,0),
           fecEstMu=na_if(fecEstMu,0)) |> 
    filter(!is.na(fecEstMu)) |> filter(!is.na(fecEstSe)) |> 
    # compute numbers of observations per plot
    group_by(plotID,species,mat,dh) |> 
    summarise(n_plot=n(),
              fecEstMu_plot=mean(fecEstMu,na.rm=TRUE),
              se_plot=sqrt(sum(fecEstSe^2)/n_plot),
              cv_plot=se_plot/fecEstMu_plot) |> 
    # Compute mean of se per plot
    ungroup() |> 
    pivot_longer(cols=c("mat","dh"),
                 names_to = "clim",
                 values_to = "clim_val")  |> 
    filter(!is.na(clim_val)) |> 
    group_by(species,clim) |>
    summarise(#opt=weighted.quantile(clim_val,1/cv_plot,prob=0.5)[1],
      qlow=weighted.quantile(clim_val,1/cv_plot,prob=0.025)[1],#,probs=0.1)[1],
      qhigh=weighted.quantile(clim_val,1/cv_plot,prob=0.975)[1],
      qrange.mastif=qhigh-qlow,
      rhigh=weighted.quantile(clim_val,1/cv_plot, prob = 1- min(5, n())/n())[1],
      rlow=weighted.quantile(clim_val,1/cv_plot, prob = min(5, n())/n())[1],
      rrange.mastif=rhigh-rlow) |> 
    ungroup() |> 
    pivot_longer(cols=c("qlow","qhigh","rhigh","rlow"),
                 names_to = "limit",
                 values_to = "value_mastif") |> 
    mutate(range.mastif=case_when(grepl("q",limit)~qrange.mastif,
                                  grepl("r",limit)~rrange.mastif)) |> 
    select(-qrange.mastif,-rrange.mastif)
  
  gbif_niche<-gbif_niche[,c("species","mat.qhigh","mat.qlow","mat.rhigh","mat.rlow",
                 "dh.qhigh","dh.qlow","dh.rhigh","dh.rlow")]|> 
    separate(species,into=c("genus","sp")) |> 
    mutate(species=paste0(substr(tolower(genus),1,8),substr(str_to_title(sp),1,8))) |> 
    select(-genus,-sp) |> 
    pivot_longer(cols=c("mat.qhigh","mat.qlow","mat.rhigh","mat.rlow",
                        "dh.qhigh","dh.qlow","dh.rhigh","dh.rlow"),
                 names_to = "limit",
                 values_to = "value_gbif") |> 
    separate(limit,into=c("clim","limit")) |> 
    mutate(type=str_sub(limit,1,1)) |> 
    group_by(species,clim,type) |> 
    mutate(range_gbif=max(value_gbif)-min(value_gbif)) 
  
  mastif.species.niche |> 
    left_join(gbif_niche,by=c("species","clim","limit")) |> 
    filter(!is.na(range_gbif)) |> 
    mutate(dev=case_when(grepl("high",limit)~100*(value_mastif-value_gbif)/range_gbif,
                         grepl("low",limit)~100*(value_gbif-value_mastif)/range_gbif)) |> 
    filter(grepl("r",limit)) |> 
    group_by(species) |> 
    mutate(ind_neg=sum(dev>(-50))) |> ungroup() |> 
    filter(ind_neg==4) |> 
    select(species) |> unique() |> pull() -> select.rank
  mastif.species.niche |> 
    left_join(gbif_niche,by=c("species","clim","limit")) |> 
    filter(!is.na(range_gbif)) |> 
    mutate(dev=case_when(grepl("high",limit)~100*(value_mastif-value_gbif)/range_gbif,
                         grepl("low",limit)~100*(value_gbif-value_mastif)/range_gbif)) |> 
    filter(grepl("q",limit)) |> 
    group_by(species) |> 
    mutate(ind_neg=sum(dev>(-50))) |> ungroup() |> 
    filter(ind_neg==4) |> 
    select(species) |> unique() |> pull() -> select.quant
  
  return(list(df.quant=gbif_niche,
              select.rank=select.rank,
              select.quant=select.quant))
}



#' Assess margin reliability based on expert maps
#' @description Compare climatic distribution of GBIF plots and MASTIF plots. 
#' Select species for which only MASTIF coverage is large enough compared to GBIF
#' @param map_file
#' @param pet_file
#' @param mat_file
#' @param species_selection_gbif
#' @param species_selection_nfi description
#' @param phylo.select 
select_species_expert<-function(map_file="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                                pet_file="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif",
                                mat_file="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                                species_selection_gbif,
                                species_selection_nfi,
                                phylo.select
                              ){
  map_raster=terra::rast(map_file)
  pet_raster=terra::rast(pet_file)
  mat_raster=terra::rast(mat_file)
  dh_raster=pet_raster-map_raster/12
  list_rast_min <- lapply(list(mat=mat_raster,
                               dh=dh_raster), 
                          function(r) {
                            aggregate(r, fact = 10,fun="min")
                            })
  list_rast_max <- lapply(list(mat=mat_raster,
                               dh=dh_raster), 
                          function(r) {
                            aggregate(r, fact = 10,fun="max")
                          })
  sf::sf_use_s2(FALSE)
  usmap <- sf::st_as_sf(getMap(resolution = "high")) |> 
    sf::st_crop(xmin=-170,xmax=-55,ymin=025,ymax=90)
  eumap <- sf::st_as_sf(getMap(resolution = "high")) |> 
    sf::st_crop(xmin=-10,xmax=35,ymin=35,ymax=70)
  
  df_extrema<-data.frame(species=phylo.select$species) |> 
    crossing(clim=c("mat","dh")) |> group_by(species,clim) |> 
    crossing(limit=c("qhigh","qlow")) |> 
    left_join(species_selection_gbif$df.quant) |> 
    left_join(species_selection_nfi$df.quant) |> 
    select(species,clim,limit,value_gbif,range_gbif,value_nfi,range.nfi) |> 
    mutate(value_exp=NA,range_exp=NA)
  
  for(sp in unique(phylo.select$species)){
    print(sp)
    if(phylo.select[phylo.select$species==sp,"block"]=="america"){
      continent.map=usmap$geometry
    }
    if(phylo.select[phylo.select$species==sp,"block"]=="europe"){
      continent.map=eumap$geometry
    }
    
    sp_map=sf::read_sf(dsn=phylo.select[phylo.select$species==sp,"dir.sp"],
                       layer=phylo.select[phylo.select$species==sp,"file.sp"]) |>
      sf::st_set_crs(sf::st_crs(continent.map)) |> 
      sf::st_intersection(continent.map)
    
    for(i in 1:2){
      terra::mask(list_rast_min[[i]],vect(sp_map))->mask_min
      terra::mask(list_rast_max[[i]],vect(sp_map))->mask_max
      q05=global(mask_min, quantile, probs=c(0.05), na.rm=TRUE)[["X5."]]
      q95=global(mask_min, quantile, probs=c(0.95), na.rm=TRUE)[["X95."]]
      df_extrema[df_extrema$species==sp &
                   df_extrema$clim==names(list_rast_min)[i]&
                   df_extrema$limit=="qhigh",c("value_exp","range_exp")]=list(q95,
                                                                              q95-q05)
      df_extrema[df_extrema$species==sp &
                   df_extrema$clim==names(list_rast_min)[i]&
                   df_extrema$limit=="qlow",c("value_exp","range_exp")]=list(q05,
                                                                             q95-q05)
    }
    
  }
  
  df_extrema<-df_extrema |> 
    mutate(dif=case_when(limit=="qhigh"~100*(value_exp-value_nfi)/range_exp,
                         limit=="qlow"~100*(value_nfi-value_exp)/range_exp)) 
  return(df_extrema)
  
}


check_selection <- function(phylo.select,
                            data_gbif,
                            fecundity.am_clim,
                            fecundity.eu_clim,
                            mastif.am,
                            mastif.eu){
  mastif_data<-rbind(mastif.am$df.tree, 
                     mastif.eu$df.tree) |> 
    mutate(dh=pet-(map/12)) |> 
    mutate(fecEstSe=na_if(fecEstSe,0),
           fecEstMu=na_if(fecEstMu,0)) |> 
    filter(!is.na(fecEstMu)) |> filter(!is.na(fecEstSe)) |> 
    # compute numbers of observations per plot
    group_by(plotID,lon,lat,species,mat,dh) |> 
    summarise(n_plot=n(),
              fecEstMu_plot=mean(fecEstMu,na.rm=TRUE),
              se_plot=sqrt(sum(fecEstSe^2)/n_plot),
              cv_plot=se_plot/fecEstMu_plot) |> 
    # Compute mean of se per plot
    ungroup() 
  
  nfi_data<-rbind(fecundity.eu_clim,
                  fecundity.am_clim) |> 
    # mutate(dh=12*pet-map) |> 
    filter(BA!=0) |> 
    filter(dh>(-2500)) |>
    group_by(plot,lon,lat,species,mat,dh) |> 
    summarise(n=n())
  sf::sf_use_s2(FALSE)
  usmap <- sf::st_as_sf(getMap(resolution = "high")) |> 
    sf::st_crop(xmin=-135,xmax=-55,ymin=18,ymax=65)
  eumap <- sf::st_as_sf(getMap(resolution = "high")) |> 
    sf::st_crop(xmin=-10,xmax=35,ymin=35,ymax=70)
  
  littlemap<-read.csv("data/USTreeAtlas-main/Little_datatable.csv")
  
  for (sp in phylo.select$species){
    tryCatch(
      {
        print(sp)
        sp_l<-phylo.select[phylo.select$species==sp,"species_l"]
        s_p<-str_replace(sp_l," ","_")
        splittle<-littlemap |> filter(Latin.Name==sp_l) |> pull(SHP..)
        if(phylo.select[phylo.select$species==sp,"block"]=="america"){
          continent.map=usmap
          species_map=sf::read_sf(dsn=paste0("data/USTreeAtlas-main/shp/",splittle,"/"),
                                  layer=splittle) |> 
            sf::st_set_crs(sf::st_crs(eumap))
        }
        if(phylo.select[phylo.select$species==sp,"block"]=="europe"){
          continent.map=eumap
          if(file.exists(paste0("data/chorological_maps_dataset/",sp_l,"/shapefiles/",s_p,"_syn_plg.dbf"))){
            file.sp=paste0(s_p,"_syn_plg")
          }else{
            file.sp=paste0(s_p,"_plg")
          }
          species_map=sf::read_sf(dsn=paste0("data/chorological_maps_dataset/",sp_l,"/shapefiles/"),
                                  layer=file.sp) |> 
            sf::st_set_crs(sf::st_crs(eumap))
        }
        contbox=sf::st_bbox(continent.map)
        mastif.sp<-mastif_data[mastif_data$species==sp,]
        gbif.sp<-data_gbif[data_gbif$species==sp_l,] |> 
          filter(decimallongitude<contbox[["xmax"]],
                 decimallongitude>contbox[["xmin"]],
                 decimallatitude<contbox[["ymax"]],
                 decimallatitude>contbox[["ymin"]])
        nfi.sp<-nfi_data[nfi_data$species==sp,]
        ggplot()+
          geom_sf(data=continent.map,fill="tan",alpha=0.4,color=NA)+
          geom_point(data=gbif.sp,aes(x=decimallongitude,y=decimallatitude),alpha=0.2,color="darkgrey",size=0.2)+
          geom_point(data=mastif.sp,aes(x=lon,y=lat),color="red",size=0.4)+
          geom_sf(data=species_map,fill=NA,color="blue")+
          theme(axis.title = element_blank(),
                axis.text = element_blank())->plot.gbif
        ggplot()+
          geom_sf(data=continent.map,fill="tan",alpha=0.4,color=NA)+
          geom_point(data=nfi.sp,aes(x=lon,y=lat),alpha=0.2,color="darkgrey",size=0.2)+
          geom_point(data=mastif.sp,aes(x=lon,y=lat),color="red",size=0.4)+
          geom_sf(data=species_map,fill=NA,color="blue")+
          theme(axis.title = element_blank(),
                axis.text = element_blank())->plot.nfi
        
        gbif.sp |> select(species,mat,dh)  |> 
          pivot_longer(cols=c("mat","dh"),names_to="clim",values_to = "quantile") |>
          group_by(clim) |> summarise(q05_gbif=quantile(quantile,probs=0.05,na.rm=TRUE),
                                      q95_gbif=quantile(quantile,probs=0.95,na.rm=TRUE)) -> quant.gbif
        nfi.sp |> select(species,mat,dh)  |> 
          pivot_longer(cols=c("mat","dh"),names_to="clim",values_to = "quantile") |>
          group_by(clim) |> summarise(q05_nfi=quantile(quantile,probs=0.05,na.rm=TRUE),
                                      q95_nfi=quantile(quantile,probs=0.95,na.rm=TRUE)) -> quant.nfi
        mastif.sp |> select(species,mat,dh) |> 
          pivot_longer(cols=c("mat","dh"),names_to="clim",values_to = "clim_val") |> 
          left_join(quant.gbif) |>
          left_join(quant.nfi) |> 
          ggplot(aes(clim_val))+
          geom_histogram(alpha=0.5)+
          geom_vline(aes(xintercept=q05_gbif,color="gbif"))+
          geom_vline(aes(xintercept=q95_gbif,color="gbif"))+
          geom_vline(aes(xintercept=q05_nfi,color="nfi"))+
          geom_vline(aes(xintercept=q95_nfi,color="nfi"))+
          facet_wrap(~clim,scales = "free")+
          theme(axis.title = element_blank(),
                axis.text.y = element_blank(),
                legend.position = "bottom")+
          labs(color="Quantiles from:")-> plot.var
        ggplot()+
          geom_hex(data=gbif.sp,aes(x=mat,y=dh))+
          scale_fill_continuous(type = "viridis") +
          geom_point(data=mastif.sp,aes(x=mat,y=dh),color="red")
        
        plot2<-cowplot::plot_grid(plot.gbif,plot.nfi, ncol = 2, labels = c("GBIF x MASTIF", "NFI x MASTIF"))
        plot3 <-cowplot::plot_grid(plot2,plot.var, nrow = 2,rel_heights =  c(1,0.5))
        ggsave(filename = paste0("species_map/",sp,".png"),width=8,height = 8.5, plot = plot3)
      },
      error=function(e)simpleError("No maps")
    )
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - NFI data fecundity ####
#' @authors Anne Baranger (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get a list of climate files
#' @param clim_list list of climatic variables of interest
get_climate<-function(clim_list=list(mat="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                                     tmin="data/CHELSA/CHELSA_bio6_1981-2010_V.2.1.tif",
                                     map="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                                     pmax="data/CHELSA/CHELSA_bio13_1981-2010_V.2.1.tif",
                                     sgdd="data/CHELSA/CHELSA_gdd5_1981-2010_V.2.1.tif",
                                     pet="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif",
                                     cmi_min="data/CHELSA/CHELSA_cmi_min_1981-2010_V.2.1.tif",
                                     cmi_mean="data/CHELSA/CHELSA_cmi_mean_1981-2010_V.2.1.tif",
                                     cmi_max="data/CHELSA/CHELSA_cmi_max_1981-2010_V.2.1.tif")){
  return(clim_list)
}

#' Load NFI plots for a targetted continent with associated climatic variables
#' and plot size
#' @param clim_list climatic variable list
#' @param continent wether "europe" or "america"
get_nfi <- function(clim_list,
                    fit,
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
  nfi.file=paste0("data/",continent,"/",fit,"/plotData.rdata")
  size.file=paste0("data/",continent,"/",fit,"/SIZE.rdata")
  load(nfi.file)
  load(size.file)
  SIZE=as.data.frame(SIZE) |> 
    rownames_to_column(var="plot")
  plotData<-plotData |> 
    rownames_to_column(var="plot") |> 
    left_join(SIZE,by="plot")
  df.nfi<-cbind(plotData,
                extract(clim,
                        y=data.frame(x=plotData$lon,
                                     y=plotData$lat))[,-1]) |> 
    mutate(wai=(map-12*pet)/(12*pet))
  
  return(df.nfi)
}

#' Compute acp prediction for each plot of nfi
#' @param df.nfi
compute_acp<-function(df.nfi,
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

#' Compute margin class for each plot of nfi
#' @param df.nfi.clim
compute_margin<-function(df.nfi.clim,
                         fecundity){
  species.margin <- fecundity |> 
    filter(BA!=0) |> 
    left_join(df.nfi.clim,by="plot") |> 
    group_by(species) |> 
    summarise(quant95=quantile(PC1,probs=0.95,na.rm=TRUE)[[1]],
              quant525=quantile(PC1,probs=0.525,na.rm=TRUE)[[1]],
              quant475=quantile(PC1,probs=0.475,na.rm=TRUE)[[1]],
              quant05=quantile(PC1,probs=0.05,na.rm=TRUE)[[1]])
  return(species.margin)
}


#' Load fecundity for each plot, including mean and sd estimation, ISP and basal area
#' @param continent
#' @param sp.select
get_fecundity <- function(continent,
                          fit,
                          sp.select){
  # load BA
  load(paste0("data/",continent,"/",fit,"/BA.rdata"))
  BA<-as.data.frame(BA)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "BA")
  
  #load ISP
  load(paste0("data/",continent,"/",fit,"/ISP.rdata"))
  ISP<-as.data.frame(ISP)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "ISP")
  
  #load fecgmmu&sd
  load(paste0("data/",continent,"/",fit,"/fecGmSd.rdata"))
  fecGmSd<-as.data.frame(summary)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmSd")
  load(paste0("data/",continent,"/",fit,"/fecGmMu.rdata"))
  fecGmMu<-as.data.frame(summary)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
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
    group_by(species) |> 
    mutate(n=sum(BA>0,na.rm = TRUE)) |> filter(n>200) |> 
    group_by(plot,species) |> 
    mutate(n_plot=n(),
           fecGmMu_plot=mean(fecGmMu,na.rm=TRUE),
           se_plot=sqrt(sum(fecGmSd^2)/n_plot),
           cv_plot=se_plot/fecGmMu_plot) |> 
    # Compute mean of se per plot
    ungroup()
  
  return(fecundity)
}

#' Separate american block
#'@description based on mean longitude and a threshold, it classifies wether a 
#'species is part of the easter or western species block in N America
#'@param fecundity.am_clim
#'@param fecundity.eu_clim

class_species <- function(fecundity.am_clim,
                          fecundity.eu_clim){
  species_cat=rbind(fecundity.eu_clim |> mutate(block="europe"),
        fecundity.am_clim|> mutate(block="america")) |>
    filter(BA!=0) |> #rm absences
    # filter(ISP<quantile(ISP,probs=0.97)) |>  # remove extreme values
    group_by(species) |> 
    summarize(mean_lon=mean(lon)) |> 
    mutate(zone=case_when(mean_lon > (-10) ~ "europe",
                          mean_lon<(-97)~"west",
                          TRUE~"east"))
  return(species_cat)
}


get_nfipred_plot<-function(continent="europe",
                           fit="fit2024",
                           chelsa_files=list(mat="data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif",
                                          map="data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif",
                                          pet="data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif"),
                           sp.select){
  
  if(fit=="fit2024"){
    ## plotdata
    plotData<-loadRData(file.path("data",continent,"inventory","plotData.rdata")) 
    if(sum(grepl("plot",colnames(plotData)))==0){
      plotData<-plotData |> tibble::rownames_to_column(var="plot")
    }
    rownames(plotData)=NULL
    ## pet, deficit, temperature
    pet <- apply(loadRData(file.path("data",continent,"inventory","pet.rdata")),
                 MARGIN = 1,
                 mean) |> 
      as.data.frame() |>  tibble::rownames_to_column()
    colnames(pet)=c("plot","pet")
    mat <- apply(loadRData(file.path("data",continent,"inventory","temp.rdata")),
                 MARGIN = 1,
                 mean) |> 
      as.data.frame() |>  tibble::rownames_to_column()
    colnames(mat)=c("plot","mat")
    map <- apply(loadRData(file.path("data",continent,"inventory","prec.rdata")),
                  MARGIN = 1,
                  mean) |> 
      as.data.frame() |>  tibble::rownames_to_column()
    colnames(map)=c("plot","map")
    dh <- apply(loadRData(file.path("data",continent,"inventory","def.rdata")),
                 MARGIN = 1,
                 mean) |> 
      as.data.frame()|>  tibble::rownames_to_column()
    colnames(dh)=c("plot","dh")
    size<-loadRData(file.path("data",continent,fit,"SIZE.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    colnames(size)=c("plot","size","cv_size")
    # gather
    list_tab=list(pet,mat,map,dh,size)
    names(list_tab)=c("pet","mat","map","dh","size")
    for(df in names(list_tab)){
      if(sum(list_tab[[df]][["plot"]]!=plotData[["plot"]])==0){
        plotData=cbind(plotData,list_tab[[df]][df])
      }else{
        plotData<-plotData |> left_join(list_tab[[df]],by="plot")
      }
    }
  }
  
  if(fit=="fit2023"){
    plotData<-loadRData(file.path("data",continent,fit,"plotData.rdata")) 
    plotData$pet <- as.numeric(terra::extract(rast(chelsa_files$pet),
                                       cbind(plotData$lon,
                                             plotData$lat))[, 1])
    plotData$map <- as.numeric(terra::extract(rast(chelsa_files$map),
                                       cbind(plotData$lon,
                                             plotData$lat))[, 1])
    plotData$mat <- as.numeric(terra::extract(rast(chelsa_files$mat),
                                       cbind(plotData$lon,
                                             plotData$lat))[, 1])
    plotData$dh <- plotData$pet - plotData$map/12
    size<-loadRData(file.path("data",continent,fit,"SIZE.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    if(sum(plotData[["plot"]]!=size[["plot"]])==0){
      plotData=cbind(plotData,size[,c("SIZE","CV")])
    }else{
      plotData<-plotData |> left_join(size,by="plot")
    }
  }
    # structure data
    BA=loadRData(file.path("data",continent,fit,"BA.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    fecGmMu=loadRData(file.path("data",continent,fit,"fecGmMu.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    fecGmSd=loadRData(file.path("data",continent,fit,"fecGmSd.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    ISP=loadRData(file.path("data",continent,fit,"ISP.rdata")) |> 
      as.data.frame() |>
      tibble::rownames_to_column(var="plot")
    list_tab=list(BA,fecGmMu,fecGmSd,ISP)
    names(list_tab)=c("BA","fecGmMu","fecGmSd","ISP")
    plotData<-plotData |> crossing(species=sp.select)
    for(df in c("BA","fecGmMu","fecGmSd","ISP")){
        df_long<-list_tab[[df]] |> 
          dplyr::select(plot,matches(sp.select)) |> 
          pivot_longer(cols = -plot,
                       names_to = "species",
                       values_to = "value") %>%
          filter(!is.na(value)) %>%
          rename_with(~df, .cols = value)
        plotData<-left_join(plotData,df_long,by=c("plot","species"))
    }
    
    plotData<-plotData |>
      mutate(fecGmMu=na_if(fecGmMu,0),
             fecGmSd=na_if(fecGmSd,0),
             ISP=na_if(ISP,0)) |>  
      filter(!is.na(ISP)) |> 
      group_by(species) |> 
      mutate(n=sum(BA>0,na.rm = TRUE)) |> filter(n>200) |> 
      group_by(plot,species) |> 
      mutate(n_plot=n(),
             fecGmMu_plot=mean(fecGmMu,na.rm=TRUE),
             se_plot=sqrt(sum(fecGmSd^2)/n_plot),
             cv_plot=se_plot/fecGmMu_plot) |> 
      # Compute mean of se per plot
      ungroup()
    
    return(plotData)
}

#' Get temporal variability
#'@description uses timeseries of fecundity to compute wieghted mean, weighted variance
#'and subsequent coefficient of variation of time series
#'@param continent
#'@param sp.select

get_cv <- function(continent,
                   sp.select){
  # get years
  filename=list.files(paste0("data/",continent))
  year=as.numeric(sub(".*-(\\d+)\\.rdata", "\\1", filename))
  year=unique(year[!is.na(year)])
  
  # get weights, mu ,sd
  for (y in year){
    load(paste0("data/",continent,"/fecGmSd-",y,".rdata"))
    sd=summary
    rm(summary)
    load(paste0("data/",continent,"/fecGmMu-",y,".rdata"))
    mu=summary
    rm(summary)
    assign(paste0("w_",y),
           mu/sd)
    assign(paste0("sd_",y),
           sd)
    assign(paste0("mu_",y),
           mu)
  }
  rm(mu,sd)
  w=array(data=unlist(lapply(ls()[grepl( "w_",ls())],
                             FUN = function(x)eval(parse(text=x)))),
          dim = c(nrow(w_2010),
                  ncol(w_2010),
                  length(ls()[grepl( "w_",ls())])))
  rm(list=ls()[grepl( "w_",ls())])
  gc()
  mu=array(data=unlist(lapply(ls()[grepl( "mu_",ls())],
                              FUN = function(x)eval(parse(text=x)))),
           dim = c(nrow(mu_2010),
                   ncol(mu_2010),
                   length(ls()[grepl( "mu_",ls())])))
  rm(list=ls()[grepl( "mu_",ls())])
  gc()
  sd=array(data=unlist(lapply(ls()[grepl( "sd_",ls())],
                              FUN = function(x)eval(parse(text=x)))),
           dim = c(nrow(sd_2010),
                   ncol(sd_2010),
                   length(ls()[grepl( "sd_",ls())])))
  rm(list=ls()[grepl( "sd_",ls())])
  gc()
  # components of weighted mean and var
  sum_w=apply(w,
              MARGIN=c(1,2),
              sum,
              na.rm=TRUE)
  sum_w2=apply(w*w,
               MARGIN=c(1,2),
               sum,
               na.rm=TRUE)
  var_mu=apply(mu,
               MARGIN=c(1,2),
               var,
               na.rm=TRUE)
  mu_w=apply(mu*w,
             MARGIN=c(1,2),
             sum,
             na.rm=TRUE)
  
  # weighted mean and var
  mu_weigthed=mu_w/sum_w
  var_weighted=(var_mu*sum_w2)/(sum_w^2)
  
  # coef of var
  cv=sqrt(var_weighted)/mu_weigthed
  
  # volatility
  load(paste0("data/",continent,"/fecGmSd-",y,".rdata"))
  mu_0=mu[,colnames(summary)%in%sp.select,]
  mu_0[is.na(mu_0)]=0
  volatility<-apply(mu_0,
                    MARGIN = c(1,2),
                    function(y){
                      if(sum(y)==0){return(0)}else{
                      mastSpectralDensity(log(y+1),maxPeriod = 6,PLOT=FALSE)$volatility}
                      }
                    )
  colnames(volatility)=colnames(summary)[colnames(summary)%in%sp.select]
  rownames(volatility)<-rownames(summary)
  

  
  # get data together
  colnames(mu_weigthed)<- colnames(summary)
  rownames(mu_weigthed)<-rownames(summary)
  colnames(var_weighted)<- colnames(summary)
  rownames(var_weighted)<-rownames(summary)
  colnames(cv)<- colnames(summary)
  rownames(cv)<-rownames(summary)
  
  mu_weigthed<-as.data.frame(mu_weigthed)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmMuW")
  var_weighted<-as.data.frame(var_weighted)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmVarW")
  cv<-as.data.frame(cv)|> 
    rownames_to_column(var="plot") |> 
    dplyr::select(plot,matches(sp.select)) |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "fecGmCv")
  volatility<-as.data.frame(volatility)|> 
    rownames_to_column(var="plot") |> 
    pivot_longer(cols=-plot,
                 names_to = "species",
                 values_to = "volatility")
  
  
  
  fecundity_cv<-mu_weigthed |> 
    left_join(var_weighted,by=c("plot","species")) |> 
    left_join(cv,by=c("plot","species")) |> 
    left_join(volatility,by=c("plot","species")) |> 
    mutate(fecGmMuW=na_if(fecGmMuW,0),
           fecGmMuW=na_if(fecGmMuW,NaN),
           fecGmVarW=na_if(fecGmVarW,0),
           fecGmVarW=na_if(fecGmVarW,NaN),
           fecGmCv=na_if(fecGmCv,0),
           fecGmCv=na_if(fecGmCv,NaN),
           volatility=na_if(volatility,0),
           volatility=na_if(volatility,NaN)) |> 
    group_by(species) |> 
    mutate(n=sum(fecGmMuW>0,na.rm=TRUE)) |> filter(n>200) |> 
    filter(species %in% sp.select) |> 
    dplyr::select(-n)
  
  return(fecundity_cv)
  }


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - Traits compilation ####
#' @authors Anne Baranger, Julien Barrère (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compile all traits data
#' @author Julien Barrère
#' @param TRY_file file containing TRY request
#' @param species.in character vector of all the species for which to extract the data
compile_traits_TRY <- function(TRY_file, NFI.traits_file){
  
  # -- Species name
  species.in <- fread(NFI.traits_file)$species
  
  # -- Translate TRY traits code into abbreviated traits name
  TRY.traits.name <- data.frame(
    TraitID = c(24, 3117, 146, 14, 56, 15, 46, 65, 1111, 2809, 2807, 2808, 159, 30, 318, 
                31, 719, 59, 819, 45, 773, 413, 324, 1229, 153, 865, 837, 3446), 
    trait = c("bark.thickness", "leaf.sla", "leaf.CN.ratio", "leaf.N.mass", "leaf.NP.ratio", 
              "leaf.P.mass", "leaf.thickness", "root.type", "seedbank.density", 
              "seedbank.duration", "seedbank.n.layers", "seedbank.thickness.toplayer", 
              "seedbank.type", "tolerance.drought", "tolerance.fire", "tolerance.frost", 
              "xylem.hydraulic.vulnerability", "plant.lifespan", "plant.resprouting.capacity", 
              "stomata.conductance", "crown.height", "leaf.Chl.content", "crown.length", 
              "wood.Nmass", "budbank.height.distribution", "budbank.seasonality", 
              "bark.structure", "plant.biomass")
  )
  
  
  
  ## -- Compile numeric traits data
  traits.TRY <-data.table::fread(TRY_file) %>%
    filter(!is.na(TraitID)) %>%
    filter(!is.na(StdValue)) %>%
    filter(AccSpeciesName %in% species.in) %>%
    left_join(TRY.traits.name, by = "TraitID") %>%
    mutate(trait = paste("TRY", trait, gsub("\\ ", "", UnitName), sep = "_"), 
           trait = gsub("\\/", ".", trait)) %>%
    rename("species" = "AccSpeciesName") %>%
    group_by(species, trait) %>%
    summarize(trait.value = mean(StdValue, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(trait) %>%
    mutate(n.species.per.trait = n()) %>%
    filter(n.species.per.trait >= 10) %>% 
    dplyr::select(-n.species.per.trait) %>%
    spread(key = trait, value = trait.value)  
  
  return(traits.TRY)
}



