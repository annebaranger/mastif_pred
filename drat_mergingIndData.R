# get european inventory
library(RPostgreSQL) 
source("R/Fundiv/Get_SamsaraDataset_from_PgAdmin.R")
data_plots <- Get_SamsaraDataset_from_PgAdmin("calib_plots", "natheo_beauchamp", "natheo")
data_tree <- Get_SamsaraDataset_from_PgAdmin("calib_trees", "natheo_beauchamp", "natheo")


#1 faire une liste des genres/espèces/continents
#2 load recursively chaque fichier et select les bonnes espèces

phylo.select.genus<-phylo.select |> 
  select(genus,taxa,block) |> filter(block=="europe") |> unique()
for (i in 1:dim(phylo.select.genus)[1]){
  sp.list<-phylo.select |> 
    filter(block==phylo.select.genus[i,"block"]&
             genus==phylo.select.genus[i,"genus"]) |> 
    select(species) |> pull()
  load(paste0("data/",phylo.select.genus[i,"block"],"/annual/",
              tolower(phylo.select.genus[i,"genus"]),".rdata"))
  inventoryFec<-inventoryFec |> 
    filter(species %in% sp.list)

}
# rm(sp.list)

# check plots in mastif pred and ifn
plots_ifn=abies |> select(plotcode) |> 
  unique() |> pull(plotcode)
plots_mastif=inventoryFec |> select(plot) |> 
  unique() |> 
  pull(plot)

# try to get back to IFN id
inventoryFec |> 
  mutate(country=str_sub(plot,1,2)) |> 
  select(country) |> unique() |> pull()

mastifspred_abies=inventoryFec |> 
  filter(year==2020) |> 
  select(plot) |> unique()

View(mastifspred_abies |> 
       group_by(plot) |> 
       mutate(ifn=sum(grepl(pattern=plot,x=unique(abies$plotcode)))))

View(mastifspred_abies |> 
       filter(!grepl("FG",plot)) |> 
  mutate(country=str_sub(plot,1,2),
         plot_new= gsub("\\.","_",plot),
         plotcode=case_when(country=="ES"~paste0("1_",str_sub(plot_new,3,-1)),
                            country=="WA"~paste0("2_",str_sub(plot_new,3,-1)),
                            country=="DE"~paste0("3_",str_sub(plot_new,3,-1)),
                            country=="SW"~paste0("4_",str_sub(plot_new,3,-1)),
                            country=="FI"~paste0("5_",str_sub(plot_new,3,-1)),
                            country=="FR"~paste0("6_",str_sub(plot_new,3,-1)),
                            country=="PL"~paste0("7_",str_sub(plot_new,3,-1)),
                            country=="SI"~paste0("8_",str_sub(plot_new,3,-1)),
                            country=="SK"~paste0("9_",str_sub(plot_new,3,-1)),
                            country=="CZ"~paste0("10_",str_sub(plot_new,3,-1)),
                            TRUE~paste0("1_",plot))) |> 
    left_join(data_plots,by="plotcode")
    ) 


## match per species
# load data_tree
load("data/IFN/data_plots_europe.rdata")
load("data/IFN/data_tree_europe.rdata")
tar_load(phylo.select,store="target_nfi")

# load data_plots
data_tree <- data_tree |> 
  mutate(treecode=gsub(" ","_",treecode)) |> 
  left_join(data_plots,by="plotcode")

data_tree |> 
  filter(is.na(plotcode)) # check if all trees were matched to a plot

# fit growth model
data_growth<-data_tree |> 
  filter(species %in% phylo.select$species_l) |> 
  filter(treestatus==2) |> 
  mutate(diftime=lubridate::time_length(difftime(time1 = surveydate2,time2=surveydate1,units="days"),"years"),
         growth=(dbh2-dbh1)/diftime) |> 
  group_by(species) |> 
  mutate(growth05=quantile(growth,probs = 0.05)) |> 
  filter(growth>quantile(growth,probs = 0.05)) |> 
  mutate(G=log(growth+growth05)) |> 
  filter(G>-10000) |> 
  ungroup() 


data_tree |> 
  filter(species %in% phylo.select$species_l) |> 
  filter(treestatus==2) |> 
  mutate(diftime=lubridate::time_length(difftime(time1 = surveydate2,time2=surveydate1,units="days"),"years"),
         growth=(dbh2-dbh1)/diftime) |> 
  group_by(species) |> 
  mutate(growth05=quantile(growth,probs = 0.05)) |> ungroup() |> 
  filter(species=="Abies alba") |>
  # filter(country=="SI") |> 
  ggplot(aes(growth))+
  geom_density() +
  xlim(-5,5)+
  geom_vline(aes(xintercept=growth05))
  
data_tree |> 
  filter(species %in% phylo.select$species_l) |> 
  filter(treestatus==2) |> 
  mutate(diftime=lubridate::time_length(difftime(time1 = surveydate2,time2=surveydate1,units="days"),"years"),
         growth=(dbh2-dbh1)/diftime) |> 
  group_by(species) |> 
  # filter(country!="SI") |> 
  summarise(growth05=quantile(growth,probs = 0.1,na.rm = TRUE)) 


species.sampling=sample(unique(data_growth$species),10)

data_growth |> 
  filter(species %in% species.sampling) |> 
  ggplot(aes(log(growth+growth05),color=species))+
  geom_density()+
  theme_minimal()

data_mod<-data_growth |> 
  filter(species %in% species.sampling) |> 
  mutate(G=scale(G))
growth.mod<-lme4::lmer(G ~ species + species:dbh1 + species:log(dbh1) + (1 | country) + (1 | country:plotcode) ,
                       data = data_mod)

cbind(data_mod,
      pred=unname(predict(growth.mod,data_mod))) |> 
  ggplot(aes(G,pred))+
  geom_point()


data_predict<-data_growth |> 
  select(plotcode,species,country,growth05) |> 
  mutate(dbh1=250) |> 
  unique()
data_predict$G=unname(predict(growth.mod,data_predict))





# fit survival model
data_survival<-data_tree |> 
  filter(species %in% phylo.select$species_l) |> 
  filter(treestatus==2) |> 
  mutate(diftime=lubridate::time_length(difftime(time1 = surveydate2,time2=surveydate1,units="days"),"years"),
         growth=(dbh2-dbh1)/diftime) |> 
  group_by(species) |> 
  mutate(growth05=quantile(growth,probs = 0.05)) |> 
  filter(growth>quantile(growth,probs = 0.05)) |> 
  mutate(G=log(growth+growth05)) |> 
  filter(G>-10000) |> 
  ungroup()

inventoryAll <- data.frame()
plots_ifn <- data_plots$plotcode
for (i in 1:dim(phylo.select.genus)[1]){
  print(i)
  sp.list<-phylo.select |> 
    filter(block==phylo.select.genus[i,"block"]&
             genus==phylo.select.genus[i,"genus"]) |> 
    select(species) |> pull()
  load(paste0("data/",phylo.select.genus[i,"block"],"/annual/",
              tolower(phylo.select.genus[i,"genus"]),".rdata"))
  inventoryFec<-inventoryFec |> 
    filter(species %in% sp.list)  
  inventoryFec <- inventoryFec |> 
    # filter(!grepl("_FG",plot)) |> 
    filter(!grepl("IT",plot)) |>
    filter(!grepl("INRA",plot)) |> 
    filter(!grepl("MASTIF",plot)) |> 
    # filter(year==2020) |>
    select(plot,treeID,lon,lat,species,diam,ecoReg,year,fecMu,fecSd,fecGmMu,fecGmSd) |> 
    mutate(country=ifelse(str_sub(plot,-2,-1)=="FG","FR",str_sub(plot,1,2)),
           plot_new= gsub("\\.","_",plot),
           plotcode=case_when(country=="ES"~paste0("1_",str_sub(plot_new,3,-1)),
                              country=="WA"~paste0("2_",str_sub(plot_new,3,-1)),
                              country=="DE"~paste0("3_",str_sub(plot_new,3,-1)),
                              country%in% c("SW","83","84","85")~paste0("4_",plot_new),
                              country=="FI"~paste0("5_",str_sub(plot_new,3,-1)),
                              country=="FR"~paste0("6_",str_sub(plot_new,1,-4)),
                              country=="PL"~paste0("7_",str_sub(plot_new,3,-1)),
                              country=="SI"~paste0("8_",str_sub(plot_new,3,-1)),
                              country=="SK"~paste0("9_",str_sub(plot_new,3,-1)),
                              country=="CZ"~paste0("10_",str_sub(plot_new,3,-1)),
                              TRUE~paste0("1_",plot))) |> 
    left_join(data_plots,by=c("plotcode"))
  not_found <- inventoryFec |> 
    filter(is.na(longitude)) |> 
    filter(year==2020)
  inventoryFec<-inventoryFec |> 
    filter(!is.na(country.y)) |>  # filter out unmatched plots (mostly french plots without remeasures 
    # and czech permanent plots)
    select(-country.x) |> 
    rename(country=country.y,
           mastif_plot=plot,
           mastif_tree=treeID)
  inventoryFR<-inventoryFec |> 
    filter(country=="FR")
  inventoryFec<-inventoryFec |> 
    filter(country!="FR") |> 
    mutate(treecode=case_when(country=="ES"~paste0("1_",gsub("\\.","_",
                                                           gsub(".+_(.*)$", "\\1",
                                                                mastif_tree))),
                              country=="WA"~paste0("2_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="DE"~paste0("3_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="SW"~paste0("4_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="FI"~paste0("5_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="FR"~paste0("6_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="PL"~paste0("7_",gsub("-","_",
                                                             str_sub(gsub(".+_(.*)$", "\\1",
                                                                          mastif_tree),
                                                                     3,-1)
                                                             )
                                                   ),
                              country=="SI"~paste0("8_",gsub("\\.","_",
                                                             gsub(".+_(.*)$", "\\1",
                                                                  mastif_tree))),
                              country=="SK"~paste0("9_",gsub("-","_",
                                                             str_sub(gsub(".+_(.*)$", "\\1",
                                                                          mastif_tree),
                                                                     3,-1)
                                                             )
                                                   ),
                              country=="CZ"~paste0("10_",gsub("-","_",
                                                              str_sub(gsub(".+_(.*)$", "\\1",
                                                                           mastif_tree),
                                                                      3,-1)
                                                              )
                                                   ),
                              TRUE~"OUT")) |> 
    # mutate(treeID=gsub("-","_",treeID),
    #        treeID=gsub(".+_(.*)$", "\\1", treeID),
    #        treeID=gsub("\\.","_",treeID),
    #        treecode=paste0(gsub("(.+?)(\\_.*)", "\\1", plotcode),"_",
    #                        treeID)) |> 
    left_join(data_tree,by=c("plotcode","treecode")) |> 
    rename(species_l=species.y,
           species=species.x)
  # pb de merge pour la france car id tree très bizarre
  # merge par id plot et diamètre, delete trees with different possible combination
  inventoryFR=inventoryFR |> 
    left_join(phylo.select[,c("species_l","species")]) |>  # to get the whol species name
    mutate(dbh1=diam*10) |> 
    left_join(data_tree,by=c("plotcode","dbh1","species_l"="species")) |> 
    mutate(idtree_mastif=str_sub(gsub("\\.","_",
                                      gsub(".+_(.*)$", "\\1",
                                           mastif_tree)),3,-1)) |> 
    group_by(treecode) |> 
    mutate(n=n()) |> 
    filter(n<=15) # delete individual for which multiple combination are possible
  inventoryFec<-inventoryFec |> 
    filter(!is.na(treestatus)) |> 
    bind_rows(inventoryFR)
  
  inventoryFec<-inventoryFec |> 
    filter(!treestatus %in% c(1,6)) |> 
    mutate(diftime=lubridate::time_length(difftime(time1 = surveydate2,time2=surveydate1,units="days"),"years"),
           growth=case_when(treestatus==2~(dbh2-dbh1)/diftime,
                            TRUE~-10),
           mortality=(treestatus!=2))
  
  save(inventoryFec,
       file=paste0("data/",phylo.select.genus[i,"block"],"/merge_ifn/",
              tolower(phylo.select.genus[i,"genus"]),".rdata"))
}


## USA
