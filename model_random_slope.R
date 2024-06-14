# script for model with part of the data, random slope
fecundity_data<-rbind(tar_read(fecundity.am_clim,store="target_nfi"),
                      tar_read(fecundity.eu_clim,store="target_nfi")) |> 
  left_join(tar_read(phylo.select,store="target_nfi")) |>
  dplyr::select(plot,lon,lat,taxa,genus,species,BA,ISP,fecGmMu,fecGmSd,pet,map,mat,block,zone) |>
  filter(BA!=0) |> #rm absences
  # mutate(dh=12*pet-map) |> 
  filter(!is.na(ISP)) |> 
  filter(species %in% tar_read(species_selection,store="target_nfi")$select.quant) |> 
  group_by(species) |> 
  # compute weighted quantiles of sgdd and wai, and correlation between wai and sgdd
  mutate(margin.temp=case_when(mat<=weighted.quantile(mat,w=BA, prob=0.25)[[1]]~"cold",
                               mat<weighted.quantile(mat,w=BA, prob=0.5+0.25/2)[[1]]&
                                 mat>weighted.quantile(mat,w=BA, prob=0.5-0.25/2)[[1]]~"midtemp",
                               mat>=weighted.quantile(mat,w=BA, prob=1-0.25)[[1]]~"hot"),
         margin.temp=factor(margin.temp,levels=c("midtemp","cold","hot"))) |> 
  ungroup()
# fecundity_data<-tar_read(fecundity.fit.25,store="target_fit_2")

fecundity_data |> 
  group_by(species) |> summarise(n=n()) |> 
  arrange(n)

species_select=c("abiesBalsamea","fagusSylvatic","quercusIlex","abiesConcolor")

fecundity_data |> 
  # filter(species%in%species_select) |> 
  filter(zone=="east") |> 
  group_by(species) |> 
  filter(ISP<quantile(ISP,probs=0.99),ISP>quantile(ISP,probs=0.01)) |>
  mutate(lisp=log(ISP/median(ISP)),
         margin.temp = factor(margin.temp, levels = c( "midtemp","cold", "hot"))) |> 
  ungroup() |> 
  filter(margin.temp%in%c("cold","hot")) -> fecundity_select

fecundity_select |> group_by(species,margin.temp) |> summarise(misp=median(lisp),mmisp=mean(lisp))

fecundity_select |> 
  # group_by(species) |> 
  # filter(ISP<quantile(ISP,probs=0.99)) |>
  ggplot(aes(lisp,color=species))+geom_density()

ggplot(fecundity_select,aes(margin.temp,lisp))+
  geom_boxplot()+
  facet_wrap(~species)

lme_mod <- lmer(lisp ~ margin.temp + (1 | species), data = fecundity_select)
lme_mod_rand <- lmer(lisp ~ margin.temp+(margin.temp | species), data = fecundity_select)

summary(lme_mod)
summary(lme_mod_rand)


## check nfi
library(terra)
pet=rast("data/CHELSA/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif")
mat=rast("data/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif")
map=rast("data/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif")
def=pet-map/12
tar_load(mastif.am,store="target_data_2")
tar_load(fecundity.am_clim,store="target_nfi_2")
phylo.select=tar_read(phylo.select,store="target_nfi_2") |> filter(block=="america")
littlemap=read.csv("data/USTreeAtlas-main/Little_datatable.csv")
sf::sf_use_s2(FALSE)
worldmap <- sf::st_as_sf(getMap(resolution = "high")) |> 
  sf::st_crop(xmin=-130,xmax=-60,ymin=25,ymax=70)
plot(worldmap)
for(sp in unique(phylo.select$species)){
  print(sp)
  tryCatch({
    sp_long=phylo.select[phylo.select$species==sp,"species_l"]
    if(sum(grepl(pattern = sp_long,littlemap$Latin.Name))==1){
      littleshp=littlemap |> filter(Latin.Name==sp_long) |> pull(SHP..)
      map=sf::st_read(dsn=file.path("data/USTreeAtlas-main/shp",littleshp),
                      layer = littleshp) |>
        st_set_crs(st_crs(worldmap))
      fecundity.am_clim |>
        filter(species==sp) |>
        ggplot()+
        geom_sf(data=worldmap,fill="grey")+
        geom_sf(data=map,color="red",fill="darkgrey",size=0.2)+
        geom_point(aes(x=lon,y=lat),color="blue",size=0.1,alpha=0.4)->plot
      ggsave(plot=plot,filename=paste0("species_map/",sp,".png"))}
  },
  error = function(e) simpleError("test error"), finally = print("Hello"))
  
  
}


for(sp in ){
  print(sp)
  tryCatch({
    sp_long=phylo.select[phylo.select$species==sp,"species_l"]
    if(sum(grepl(pattern = sp_long,littlemap$Latin.Name))==1){
      littleshp=littlemap |> filter(Latin.Name==sp_long) |> pull(SHP..)
      map=sf::st_read(dsn=file.path("data/USTreeAtlas-main/shp",littleshp),
                      layer = littleshp) |>
        st_set_crs(st_crs(worldmap))
      fecundity.am_clim |>
        filter(species==sp) |>
        ggplot()+
        geom_sf(data=worldmap,fill="grey")+
        geom_sf(data=map,color="red",fill="darkgrey",size=0.2)+
        geom_point(aes(x=lon,y=lat),color="blue",size=0.1,alpha=0.4)->plot
      ggsave(plot=plot,filename=paste0("species_map/",sp,".png"))}
  },
  error = function(e) simpleError("test error"), finally = print("Hello"))
  
  
}
