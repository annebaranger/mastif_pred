#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_fit.R  
#' @description R script containing all functions relative to model 
#' fitting
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 1. Format dataset for models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Format dataset 
#' @param fecundity.eu_clim fecundity pred on nfi data, europe
#' @param fecundity.am_clim fecundity pred on nfi data, america
#' @param phylo.select information on species zone
raw_data <- function(fecundity.eu_clim, 
                     fecundity.am_clim,
                     phylo.select,
                     species_selection,
                     thresh=0.1){
  fec_tot=rbind(fecundity.eu_clim |> mutate(block="europe"),
                fecundity.am_clim|> mutate(block="america")) |>
    left_join(phylo.select) |>
    dplyr::select(plot,lon,lat,taxa,genus,species,BA,ISP,fecGmMu,fecGmSd,dh,pet,map,mat,block,zone) |>
    filter(BA!=0) |> #rm absences
    # mutate(dh=12*pet-map) |> 
    filter(!is.na(ISP)) |> 
    filter(species %in% species_selection) |> 
    group_by(species) |> 
    # compute weighted quantiles of sgdd and wai, and correlation between wai and sgdd
    mutate(margin.temp=case_when(mat<=weighted.quantile(mat,w=BA, prob=thresh)[[1]]~"cold",
                                 mat<weighted.quantile(mat,w=BA, prob=0.5+thresh/2)[[1]]&
                                   mat>weighted.quantile(mat,w=BA, prob=0.5-thresh/2)[[1]]~"midtemp",
                                 mat>=weighted.quantile(mat,w=BA, prob=1-thresh)[[1]]~"hot"),
           margin.temp=factor(margin.temp,levels=c("midtemp","cold","hot")),
           margin.deficit=case_when(dh<weighted.quantile(dh,w=BA, prob=thresh)[[1]]~"humid",
                                    dh<weighted.quantile(dh,w=BA, prob=0.5+thresh/2)[[1]]&
                                      dh>weighted.quantile(dh,w=BA, prob=0.5-thresh/2)[[1]]~"midhum",
                                    dh>weighted.quantile(dh,w=BA, prob=1-thresh)[[1]]~"arid"),
           margin.deficit=factor(margin.deficit,levels=c("midhum","humid","arid"))) |> 
    ungroup() #|> 
    # mutate(ISP=fecGmMu/BA,
    #        s_ISP=fecGmSd/(BA^2))
  ecoregions=sf::read_sf(dsn="data/WWF/official", #read ecoregions
                         layer="wwf_terr_ecos") |>
    select(BIOME,ECO_NAME,geometry) |> 
    filter(!BIOME%in%c(98,99))
  sf::sf_use_s2(FALSE)
  points=sf::st_as_sf(fec_tot[,c("lon","lat")],coords=c("lon","lat"),crs=sf::st_crs(ecoregions))
  points_biome=sf::st_join(points,ecoregions) 
  
  fec_tot= cbind(fec_tot,
        biome=points_biome$BIOME) 
  
  return(fec_tot)
}


#' Associate each species to one biome
#' @param fecundity.fit fec data wtih extracted biome 
#' @param threshold threshold for prevalence below which biome are not selected 
#' for a species
get_biome <- function(fecundity.fit,
                      threshold=0.1){
  sp.biome <- fecundity.fit |> 
    group_by(species) |> mutate(n_sp=n()) |> 
    group_by(n_sp,species,biome) |>
    summarise(n=n()/n_sp) |> unique() |> ungroup() |> 
    group_by(species) |> 
    arrange(desc(n), .by_group = TRUE) %>%
    slice_head(n = 2) |> 
    filter(n>threshold) |> 
    ungroup() |> 
    select(species,biome,n) |> 
    mutate(biome=case_when(biome==13 ~ 12,
                           TRUE ~ biome))
  return(sp.biome)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 2. Fit continent models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit.continent.discrete<-function(data_fit,
                                 folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_cont<-data_fit |> 
    filter(!is.na(ISP)) |> 
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]]) |>
    ungroup() |>
    filter(!is.na(margin.temp)) |> 
    mutate(dISP=ISP/mid_isp) |> 
    filter(!is.na(dISP))

  X=model.matrix(~margin.temp*zone, fec_cont)
  
  if(file.exists(file.path(folder,"anova_allsp_margintemp_zone.RData"))){
    load(file.path(folder,"anova_allsp_margintemp_zone.RData"))
  }else{
    data_list=list(N=dim(fec_cont)[1],
                   NX=ncol(X),
                   S=nlevels(as.factor(fec_cont$species)),
                   species=as.numeric(as.factor(fec_cont$species)),
                   ISP=fec_cont$dISP,
                   X=X)
    
    fit <- stan(file = "stan/lmm_dif.stan",
                   data=data_list,
                   iter=1000,
                   chains=3,
                   core=3)
    save(fit,file=file.path(folder,"anova_allsp_margintemp_zone.RData"))
    
  }
  
  # launch_shinystan(fit.lm)
  posteriors_fec<-as.data.frame(fit) |> 
    dplyr::select(!matches("beta_")) |> 
    dplyr::select(matches("beta"))
  colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
  
  
  # contrast.max=
  contrast_cor=as.data.frame(X) |> 
    unique() |>
    mutate(margin.temp=case_when(margin.tempcold==1~"cold",
                                 margin.temphot==1~"hot",
                                 TRUE~"mid"),
           zone=case_when(zoneeurope==1~"europe",
                          zonewest==1~"west",
                          TRUE~"east")) |> 
    mutate(margin.temp=factor(margin.temp,levels=c("cold","mid","hot"))) |> 
    tibble::rownames_to_column(var = "code") |> 
    dplyr::select(code,margin.temp,zone) |> 
    mutate(contrast=paste0(margin.temp,"_",zone))
  
  contrast_mat= X |> unique() |>  t()
  contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat
  
  post_distrib=as.data.frame(contrast_post) |> 
    pivot_longer(cols=everything(),
                 names_to = "code",
                 values_to = "posterior") |>
    left_join(contrast_cor,by="code") |> 
    mutate(zone=factor(zone,levels=c("west","east","europe"))) |> 
    select(-code)
  
  return(post_distrib)

}



fit.continent.discrete.excl<-function(data_fit,
                                      folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_cont<-data_fit |> 
    filter(!is.na(ISP)) |> 
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]]) |>
    ungroup() |>
    filter(!is.na(margin.temp)) |> 
    mutate(dISP=ISP/mid_isp) |> 
    filter(!is.na(dISP))
  
  fec_cont |> 
    filter(!is.na(dh)) |> 
    mutate(dh=12*(pet-map)) |> 
    group_by(species) |> 
    mutate(dh.cat=cut(dh,breaks = seq(min(dh),max(dh),by=100)),
           cormatdh=cor(mat,dh,use="complete")) |>
    group_by(species,cormatdh,margin.temp,dh.cat) |> 
    summarise(n=n()) |> 
    filter(margin.temp!="midtemp") |> 
    pivot_wider(names_from = margin.temp,
                values_from = n) |> 
    filter(!is.na(dh.cat)) |> 
    mutate(across(c("cold","hot"),
                  ~replace_na(.,0))) |> 
    mutate(d=abs(cold-hot)/(hot+cold)) |> 
    summarise(n=mean(d)) |> arrange(n) |> 
    filter(n>0.7) |> pull(species)->speciescor
  

  fec_cont<-fec_cont[!fec_cont$species%in% speciescor,]
  
  X=model.matrix(~margin.temp*zone, fec_cont)
  
  
  if(file.exists(file.path(folder,"anova_allsp_margintemp_zone_spselect.RData"))){
    load(file.path(folder,"anova_allsp_margintemp_zone_spselect.RData"))
  }else{
    data_list=list(N=dim(fec_cont)[1],
                   NX=ncol(X),
                   S=nlevels(as.factor(fec_cont$species)),
                   species=as.numeric(as.factor(fec_cont$species)),
                   ISP=fec_cont$dISP,
                   X=X)
    
    fit <- stan(file = "stan/lmm_dif.stan",
                data=data_list,
                iter=1000,
                chains=3,
                core=3)
    save(fit,file=file.path(folder,"anova_allsp_margintemp_zone_spselect.RData"))
  }
  fit=fit.lm
  
  
  # launch_shinystan(fit.lm)
  posteriors_fec<-as.data.frame(fit) |> 
    dplyr::select(!matches("beta_")) |> 
    dplyr::select(matches("beta"))
  colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
  
  
  # contrast.max=
  contrast_cor=as.data.frame(X) |> 
    unique() |>
    mutate(margin.temp=case_when(margin.tempcold==1~"cold",
                                 margin.temphot==1~"hot",
                                 TRUE~"mid"),
           zone=case_when(zoneeurope==1~"europe",
                          zonewest==1~"west",
                          TRUE~"east")) |> 
    mutate(margin.temp=factor(margin.temp,levels=c("cold","mid","hot"))) |> 
    tibble::rownames_to_column(var = "code") |> 
    dplyr::select(code,margin.temp,zone) |> 
    mutate(contrast=paste0(margin.temp,"_",zone))
  
  contrast_mat= X |> unique() |>  t()
  contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat

  
  post_distrib=as.data.frame(contrast_post) |> 
    pivot_longer(cols=everything(),
                 names_to = "code",
                 values_to = "posterior") |>
    left_join(contrast_cor,by="code") |> 
    mutate(zone=factor(zone,levels=c("west","east","europe"))) |> 
    select(-code)
  
  return(post_distrib)
  
}


fit.continent.continous<-function(data_fit,
                                  excluded=TRUE,
                                  folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_cont<-data_fit |> 
    filter(!is.na(ISP)) |> 
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]]) |>
    ungroup() |>
    filter(!is.na(margin.temp)) |> 
    mutate(dISP=ISP/mid_isp) |> 
    filter(!is.na(dISP)) |> 
    filter(dh>(-2000)) |>
    mutate(dh_old=dh,
           dh=scale(dh,center=TRUE,scale=FALSE),
           dh2=dh^2) |> 
    filter(!is.na(dh))
    
  out_ancova_cont=fec_cont |> 
    mutate(fit=NA,
           quant05=NA,
           quant95=NA,
           mean=NA) |> 
    slice(0)
  speciesexcl=c()
  for (z in c("europe","east","west")){
    print(z)
    # filter data
    sp.zone=fec_cont[fec_cont$zone==z,] 
    
    if(excluded==TRUE){
      sp.zone |> 
        mutate(dh_old=12*dh_old) |> #because expressed along year
        group_by(species) |> 
        mutate(dh.cat=cut(dh_old,breaks = seq(min(dh_old),max(dh_old),by=100)),
               cormatdh=cor(mat,dh_old,use="complete")) |>
        group_by(species,cormatdh,margin.temp,dh.cat) |> 
        summarise(n=n()) |> 
        filter(margin.temp!="midtemp") |> 
        pivot_wider(names_from = margin.temp,
                    values_from = n) |> 
        filter(!is.na(dh.cat)) |> 
        mutate(across(c("cold","hot"),
                      ~replace_na(.,0))) |> 
        mutate(d=abs(cold-hot)/(hot+cold)) |> 
        summarise(n=mean(d)) |> arrange(n) |> 
        filter(n>0.7) |> pull(species)->speciescor
      # fit model
      sp.zone=sp.zone |> 
        filter(!species %in% speciescor)
      speciesexcl=c(speciesexcl,speciescor)
    }
    
    # # plot of excluded species
    # sp.zone |> 
    #   filter(species %in% speciescor) |>
    #   filter(margin.temp!="midtemp") |> 
    #   ggplot(aes(dh_old,log(ISP),color=margin.temp))+
    #   geom_point(alpha=0.4)+
    #   geom_boxplot(aes(dh_old,0,color=margin.temp),outlier.colour = NA)+
    #   geom_smooth()+
    #   scale_color_manual(values=c("steelblue1","firebrick1"))+
    #   facet_wrap(~species)->plot
    # ggsave(filename = paste0(folder,"spExcl_",z,".png"),
    #        plot)
    
    X=model.matrix(~margin.temp*dh, sp.zone)
    model=paste0("ancova_",z,"_lin.RData")

   
    if(file.exists(file.path(folder,model))){
      load(file.path(folder,model))
    }else{
      data_zone=list(N=dim(sp.zone)[1],
                     NX=ncol(X),
                     S=nlevels(as.factor(sp.zone$species)),
                     species=as.numeric(as.factor(sp.zone$species)),
                     ISP=sp.zone$dISP,
                     X=X,
                     NULL)
      fit.ancova_lin <- stan(file = "stan/lmm_ancova.stan",
                             data=data_zone,
                             iter=1000,
                             chains=3,
                             core=3,
                             include=FALSE,
                             pars=c("beta_raw"))
      save(fit.ancova_lin,file=file.path(folder,model))    
    }
    fit=fit.ancova_lin
    
    posteriors_fec<-as.data.frame(fit) |>
      dplyr::select(!matches("beta_")) |>
      dplyr::select(matches("beta"))
    colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
    speciescode=sp.zone |> ungroup()|> dplyr::select(species) |> mutate(code=as.numeric(as.factor(species))) |> unique()
    speciesranef=as.data.frame(summary(fit)$summary) |> 
      rownames_to_column(var="parameter") |> 
      filter(grepl("alpha",parameter)) |> 
      dplyr::select(parameter,mean) |> 
      mutate(code=as.numeric(gsub(".*?\\[([^]]*)\\].*", "\\1", parameter))) |> 
      left_join(speciescode) |> 
      dplyr::select(species,mean)
    
    contrast_mat= X |>  t()
    contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat
    fit=apply(contrast_post,2,median)
    quant05=apply(contrast_post,2,quantile,probs=0.05)
    quant95=apply(contrast_post,2,quantile,probs=0.95)
    
    output=cbind(sp.zone,
                 fit=fit,
                 quant05=quant05,
                 quant95=quant95) |>
      left_join(speciesranef) 
    
    out_ancova_cont=rbind(out_ancova_cont,
                          output)
  }
  
  return(list(posterior=out_ancova_cont,
              speciesexcl=speciesexcl))
  
}








#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 3. Fit biomes models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit.biome.continuous<-function(data_fit,
                               species.biome,
                               folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_biome<-data_fit |> 
    filter(!is.na(ISP)) |> 
    filter(!is.na(dh)) |>
    select(- biome) |> 
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]]) |>
    ungroup() |>
    filter(!is.na(margin.temp)) |> 
    left_join(species.biome) |> 
    filter(biome!=8) |> 
    mutate(dISP=ISP/mid_isp) |> 
    filter(!is.na(dISP))
  
  # output file
  out_anova_biome=fec_biome |> 
    mutate(fit=NA,
           quant05=NA,
           quant95=NA,
           mean=NA) |> 
    slice(0)
  
  for(b in unique(fec_biome$biome)){
    print(b)
    sub_biome=fec_biome |> 
      filter(biome==b) |> 
      filter(!is.na(species)) |> 
      ungroup()
    
    X=model.matrix(~margin.temp*dh, sub_biome)
    model=paste0("anova_",b,".RData")
    
    if(file.exists(file.path(folder,model))){
      load(file.path(folder,model))
    }else{
      data_biome=list(N=dim(sub_biome)[1],
                      NX=ncol(X),
                      S=nlevels(as.factor(sub_biome$species)),
                      species=as.numeric(as.factor(sub_biome$species)),
                      ISP=sub_biome$dISP,
                      X=X,
                      NULL)
      fit <- stan(file = "stan/lmm_ancova.stan",
                  data=data_biome,
                  iter=1000,
                  chains=3,
                  core=3,
                  include=FALSE,
                  pars=c("beta_raw"))
      save(fit,file=file.path(folder,model))    
    }
    posteriors_fec<-as.data.frame(fit) |>
      dplyr::select(!matches("beta_")) |>
      dplyr::select(matches("beta"))
    colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
    speciescode=sub_biome |> ungroup()|> dplyr::select(species) |> mutate(code=as.numeric(as.factor(species))) |> unique()
    speciesranef=as.data.frame(summary(fit)$summary) |> 
      rownames_to_column(var="parameter") |> 
      filter(grepl("alpha",parameter)) |> 
      dplyr::select(parameter,mean) |> 
      mutate(code=as.numeric(gsub(".*?\\[([^]]*)\\].*", "\\1", parameter))) |> 
      left_join(speciescode) |> 
      dplyr::select(species,mean)
    
    contrast_mat= X |>  t()
    contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat
    fit=apply(contrast_post,2,median)
    quant05=apply(contrast_post,2,quantile,probs=0.05)
    quant95=apply(contrast_post,2,quantile,probs=0.95)
    
    output=cbind(sub_biome,
                 fit=fit,
                 quant05=quant05,
                 quant95=quant95) |>
      left_join(speciesranef,by="species") 
    
    out_anova_biome=rbind(out_anova_biome,
                          output)
  }
  
  return(out_anova_biome)
  }


fit.biome.discrete<-function(data_fit,
                             species.biome,
                             folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_biome<-data_fit |> 
    filter(!is.na(ISP)) |> 
    filter(!is.na(dh)) |>
    select(- biome) |> 
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]]) |>
    ungroup() |>
    filter(!is.na(margin.temp)) |> 
    filter(!is.na(margin.deficit)) |>
    left_join(species.biome) |> 
    filter(biome!=8) |> 
    mutate(dISP=ISP/mid_isp) |> 
    filter(!is.na(dISP))
  
  # output file
  out_anova_biome=fec_biome |> 
    mutate(fit=NA,
           quant05=NA,
           quant95=NA,
           mean=NA) |> 
    slice(0)
  
  for(b in unique(fec_biome$biome)){
    print(b)
    sub_biome=fec_biome |> 
      filter(biome==b) |> 
      filter(!is.na(species)) |> 
      ungroup()
    
    X=model.matrix(~margin.temp*margin.deficit, sub_biome)
    model=paste0("anova_",b,".RData")
    if(file.exists(file.path(folder,model))){
      load(file.path(folder,model))
    }else{
      data_biome=list(N=dim(sub_biome)[1],
                      NX=ncol(X),
                      S=nlevels(as.factor(sub_biome$species)),
                      species=as.numeric(as.factor(sub_biome$species)),
                      ISP=sub_biome$dISP,
                      X=X,
                      NULL)
      fit <- stan(file = "stan/lmm_ancova.stan",
                  data=data_biome,
                  iter=1000,
                  chains=3,
                  core=3,
                  include=FALSE,
                  pars=c("beta_raw"))
      save(fit,file=file.path(folder,model))    
    }
    posteriors_fec<-as.data.frame(fit) |>
      dplyr::select(!matches("beta_")) |>
      dplyr::select(matches("beta"))
    colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
    
    
    speciescode=sub_biome |> ungroup()|> dplyr::select(species) |> mutate(code=as.numeric(as.factor(species))) |> unique()
    speciesranef=as.data.frame(summary(fit)$summary) |> 
      rownames_to_column(var="parameter") |> 
      filter(grepl("alpha",parameter)) |> 
      dplyr::select(parameter,mean) |> 
      mutate(code=as.numeric(gsub(".*?\\[([^]]*)\\].*", "\\1", parameter))) |> 
      left_join(speciescode) |> 
      dplyr::select(species,mean)
    
    contrast_mat= X |>  t()
    contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat
    fit=apply(contrast_post,2,median)
    quant05=apply(contrast_post,2,quantile,probs=0.05)
    quant95=apply(contrast_post,2,quantile,probs=0.95)
    
    output=cbind(sub_biome,
                 fit=fit,
                 quant05=quant05,
                 quant95=quant95) |>
      left_join(speciesranef,by="species") 
    
    out_anova_biome=rbind(out_anova_biome,
                          output)
  }
  
  return(out_anova_biome)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 4. Fit species models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fit.species <- function(data_fit,
                        # model.type,
                        folder){
  if(!dir.exists(folder)) dir.create(folder)
  fec_species<-data_fit |>  
    group_by(species) |>
    mutate(mid_isp=median(ISP)[[1]],
           dISP=ISP/mid_isp,
           mean.dh=mean(dh,na.rm=TRUE),
           mean.mat=mean(mat,na.rm=TRUE),
           cor.mat.dh=cor(dh,mat,use="complete.obs")) |>
    ungroup() |> 
    filter(dh>quantile(dh,probs = 0.01,na.rm=TRUE)) |>
    filter(dh<quantile(dh,probs=0.99,na.rm=TRUE)) |> 
    # filter(dh>(-2000)) |> 
    filter(!is.na(margin.temp)) |> 
    filter(!is.na(dISP)) |> 
    filter(!is.na(dh))
  
  out_ancova_species=fec_species |> 
    mutate(fit=NA,
           quant05=NA,
           quant95=NA) |> 
    slice(0)
  # sp=unique(out_ancova_species$species)[1]
  for (sp in unique(fec_species$species)){
    print(sp)
    sp.mod<-fec_species[fec_species$species==sp,]
    model=paste0("ancova_",sp,"_lin.RData")
    X=model.matrix(~margin.temp*dh, sp.mod) 
    if(file.exists(file.path(folder,model))){
      load(file.path(folder,model))
    }else{
      data_sp=list(N=dim(sp.mod)[1],
                     NX=ncol(X),
                     ISP=sp.mod$dISP,
                     X=X,
                     NULL)
      fit <- stan(file = "stan/lmm_species.stan",
                  data=data_sp,
                  iter=1000,
                  chains=3,
                  core=3,
                  include=FALSE,
                  pars=c("beta_raw"))
      save(fit,file=file.path(folder,model))    
    }
    posteriors_fec<-as.data.frame(fit) |>
      dplyr::select(!matches("beta_")) |>
      dplyr::select(matches("beta"))
    colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
    
    contrast_mat= X |>  t()
    contrast_post=as.matrix(posteriors_fec)[,1:ncol(X)] %*% contrast_mat
    fit=apply(contrast_post,2,median)
    quant05=apply(contrast_post,2,quantile,probs=0.05)
    quant95=apply(contrast_post,2,quantile,probs=0.95)
    
    output=cbind(sp.mod,
                 fit=fit,
                 quant05=quant05,
                 quant95=quant95) 
    
    out_ancova_species=rbind(out_ancova_species,
                             output)
  }
  return(out_ancova_species)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- 5. Fit taxa models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fit.taxa<-function(data_fit,
                   folder){
  if(!dir.exists(folder)) dir.create(folder)
  
}