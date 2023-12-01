### loop for fitting models with continuous
sp.ancova=rbind(fecundity.eu_clim |> mutate(block="europe"),
                fecundity.am_clim|> mutate(block="america")) |>
  left_join(phylo.select) |> 
  dplyr::select(plot,lon,lat,species,BA,fecGmMu,fecGmSd,sgdd,wai,pet,map,mat,block,zone) |> 
  filter(BA!=0) |> #rm absences
  mutate(dh=12*pet-map) |> 
  group_by(species) |> 
  # compute weighted quantiles of sgdd and wai, and correlation between wai and sgdd
  mutate(margin.temp=case_when(mat<=weighted.quantile(mat,w=BA, prob=0.1)[[1]]~"cold",
                               mat<weighted.quantile(mat,w=BA, prob=0.55)[[1]]&
                                 mat>weighted.quantile(mat,w=BA, prob=0.45)[[1]]~"midtemp",
                               mat>=weighted.quantile(mat,w=BA, prob=0.9)[[1]]~"hot"),
         margin.temp=factor(margin.temp,levels=c("midtemp","cold","hot"))
         ) |>
  filter(dh>quantile(dh,probs = 0.01,na.rm=TRUE)) |> filter(dh<quantile(dh,probs=0.99,na.rm=TRUE)) |> 
  filter(!is.na(margin.temp)) |> 
  # filter(!is.na(margin.temp)&!is.na(margin.deficit)) |>
  # ungroup() |> group_by(species,margin.temp) |> 
  mutate(dh_old=dh,
         dh=scale(dh),
         dh2=dh^2) |> 
  ungroup() |> 
  mutate(ISP=fecGmMu/BA,
         s_ISP=fecGmSd/(BA^2))|> 
  filter(!is.na(ISP)) |> 
  filter(!is.na(dh))

folder="mod_ancova_scale/"
out_ancova_cont=sp.ancova |> 
  mutate(fit=NA,
         quant05=NA,
         quant95=NA,
         mean=NA) |> 
  slice(0)
for (z in c("europe","east","west")){
  print(z)
  # filter data
  sp.zone=sp.ancova |> 
    filter(zone==z)
  
  sp.zone |> 
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
  
  # plot of excluded species
  sp.zone |> 
    filter(species %in% speciescor) |>
    filter(margin.temp!="midtemp") |> 
    ggplot(aes(dh_old,log(ISP),color=margin.temp))+
    geom_point(alpha=0.4)+
    geom_boxplot(aes(dh_old,0,color=margin.temp),outlier.colour = NA)+
    geom_smooth()+
    scale_color_manual(values=c("steelblue1","firebrick1"))+
    facet_wrap(~species)->plot
  ggsave(filename = paste0(folder,"spExcl_",z,".png"),
         plot)
  
  # fit model
  sp.mod=sp.zone |> 
    filter(!species %in% speciescor)
  

  
  if(file.exists(paste0(folder,"ancova_",z,".RData"))){
    load(paste0(folder,"ancova_",z,".RData"))
  }else{
    X=model.matrix(~margin.temp*dh+margin.temp*dh2, sp.mod)
    data_zone=list(N=dim(sp.mod)[1],
                   NX=ncol(X),
                   S=nlevels(as.factor(sp.mod$species)),
                   species=as.numeric(as.factor(sp.mod$species)),
                   ISP=sp.mod$ISP,
                   X=X,
                   NULL)
    fit.ancova <- stan(file = "lmm_ancova.stan",
                       data=data_zone,
                       iter=1000,
                       chains=3,
                       core=3,
                       include=FALSE,
                       pars=c("beta_raw"))
    save(fit.ancova,file=paste0(folder,"ancova_",z,".RData"))    
  }

  if(file.exists(paste0(folder,"ancova_",z,"_lin.RData"))){
    load(paste0(folder,"ancova_",z,"_lin.RData"))
  }else{
    X=model.matrix(~margin.temp*dh, sp.mod)
    data_zone=list(N=dim(sp.mod)[1],
                   NX=ncol(X),
                   S=nlevels(as.factor(sp.mod$species)),
                   species=as.numeric(as.factor(sp.mod$species)),
                   ISP=sp.mod$ISP,
                   X=X,
                   NULL)
    fit.ancova_lin <- stan(file = "lmm_ancova.stan",
                       data=data_zone,
                       iter=1000,
                       chains=3,
                       core=3,
                       include=FALSE,
                       pars=c("beta_raw"))
    save(fit.ancova_lin,file=paste0(folder,"ancova_",z,"_lin.RData"))    
  }
  bic_quad=-2*summary(fit.ancova)$summary["lp__","mean"]+8*dim(sp.mod)[1]
  bic_lin=-2*summary(fit.ancova_lin)$summary["lp__","mean"]+5*dim(sp.mod)[1]
  
  if(bic_quad<bic_lin){
    fit=fit.ancova
    X=model.matrix(~margin.temp*dh+margin.temp*dh2, sp.mod)
  }else{
    fit=fit.ancova_lin
    X=model.matrix(~margin.temp*dh, sp.mod)
  }
  
  posteriors_fec<-as.data.frame(fit) |>
    dplyr::select(!matches("beta_")) |>
    dplyr::select(matches("beta"))
  colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
  speciescode=sp.mod |> ungroup()|> dplyr::select(species) |> mutate(code=as.numeric(as.factor(species))) |> unique()
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
  
  output=cbind(sp.mod,
               fit=fit,
               quant05=quant05,
               quant95=quant95) |>
    left_join(speciesranef) 

  out_ancova_cont=rbind(out_ancova_cont,
                   output)
  rm(output,quant05,quant95,fit,contrast_post,contrast_mat,speciesranef,
     speciescode,posteriors_fec,data_zone,sp.mod,sp.zone)
}

out_ancova_cont |> 
  ggplot()+
  geom_line(aes(dh,fit,color=margin.temp),size=0.6)+
  geom_point(aes(dh,log(ISP)-mean,color=margin.temp),size=0.6,alpha=0.6)+
  geom_smooth(aes(dh,log(ISP)-mean,color=margin.temp),method="gam")+
  geom_ribbon(aes(x=dh,ymin=quant05,ymax=quant95,fill=margin.temp),alpha=0.8)+
  scale_fill_manual(values=c("bisque2","steelblue1","firebrick1"))+
  scale_color_manual(values=c("bisque2","steelblue1","firebrick1"))+
  facet_wrap(~zone)
# out_ancova_cont |> 
#   group_by(zone,margin.temp,plot) |> summarise(BAtot=sum(BA)) |> ungroup() |> 
#   ggplot(aes(BAtot,color=margin.temp))+
#   geom_density()+
#   scale_color_manual(values=c("bisque2","steelblue1","firebrick1"))+
#   facet_wrap(~zone)+
#   scale_x_log10()
out_ancova_cont |> 
  mutate(margin.temp=factor(margin.temp,levels=c("cold","midtemp","hot"))) |> 
  ggplot()+
  geom_line(aes(dh,fit,color=zone),size=0.6)+
  geom_point(aes(dh,log(ISP)-mean,color=zone),size=0.6,alpha=0.1)+
  geom_ribbon(aes(x=dh,ymin=quant05,ymax=quant95,fill=zone),alpha=0.8)+
  geom_boxplot(aes(dh,y=-10,color=zone),position=position_dodge(width=1),outlier.colour = NA)+
  geom_boxplot(aes(x=4,y=log(ISP)-mean,color=zone),position=position_dodge(width=1),outlier.colour = NA)+
  scale_fill_manual(values=c("firebrick1","slateblue","firebrick4"))+
  scale_color_manual(values=c("firebrick1","slateblue","firebrick4"))+
  theme(axis.text.x = element_text(colour = inferno(6)))+
  facet_wrap(~margin.temp)




### Analysis of common genuses
sp.data=rbind(fecundity.eu_clim |> mutate(block="europe"),
              fecundity.am_clim|> mutate(block="america")) |>
  filter(BA!=0) |> dplyr::select(species) |> unique() |> pull()

genus.com<-phylo.select |> filter(species %in% sp.data) |> dplyr::select(genus,zone) |> unique() |> 
  group_by(genus) |> 
  summarise(nsp=n()) |> 
  filter(nsp==3) |> 
  pull(genus)


sp.genus=rbind(fecundity.eu_clim |> mutate(block="europe"),
                fecundity.am_clim|> mutate(block="america")) |>
  left_join(phylo.select) |> 
  dplyr::select(plot,lon,lat,genus,species,BA,fecGmMu,fecGmSd,sgdd,wai,pet,map,mat,block,zone) |> 
  filter(genus %in% genus.com) |> 
  filter(BA!=0) |> #rm absences
  mutate(dh=12*pet-map) |> 
  group_by(species) |> 
  # compute weighted quantiles of sgdd and wai, and correlation between wai and sgdd
  mutate(margin.temp=case_when(mat<=weighted.quantile(mat,w=BA, prob=0.1)[[1]]~"cold",
                               mat<weighted.quantile(mat,w=BA, prob=0.55)[[1]]&
                                 mat>weighted.quantile(mat,w=BA, prob=0.45)[[1]]~"midtemp",
                               mat>=weighted.quantile(mat,w=BA, prob=0.9)[[1]]~"hot"),
         margin.temp=factor(margin.temp,levels=c("midtemp","cold","hot"))
  ) |>
  filter(dh>quantile(dh,probs = 0.01,na.rm=TRUE)) |> filter(dh<quantile(dh,probs=0.99,na.rm=TRUE)) |> 
  filter(!is.na(margin.temp)) |> 
  # filter(!is.na(margin.temp)&!is.na(margin.deficit)) |>
  ungroup() |> group_by(genus) |> 
  mutate(dh_old=dh,
         dh=scale(dh),
         dh2=dh^2) |> 
  ungroup() |> 
  mutate(ISP=fecGmMu/BA,
         s_ISP=fecGmSd/(BA^2))|> 
  filter(!is.na(ISP)) |> 
  filter(!is.na(dh))

folder="mod_ancova_genus_2/"
out_ancova=sp.genus |> 
  mutate(fit=NA,
         quant05=NA,
         quant95=NA,
         mean=NA) |> 
  slice(0)
for(g in genus.com){
  print(g)
  for (z in c("europe","east","west")){
    print(z)
    # filter data
    sp.mod=sp.genus |> 
      filter(genus==g) |> 
      filter(zone==z)
      
    if(file.exists(paste0(folder,"ancova_",z,"_",g,".RData"))){
      load(paste0(folder,"ancova_",z,"_",g,".RData"))
    }else{
      X=model.matrix(~margin.temp*dh+margin.temp*dh2, sp.mod)
      data_zone=list(N=dim(sp.mod)[1],
                     NX=ncol(X),
                     S=nlevels(as.factor(sp.mod$species)),
                     species=as.numeric(as.factor(sp.mod$species)),
                     ISP=sp.mod$ISP,
                     X=X,
                     NULL)
      fit.ancova <- stan(file = "lmm_ancova.stan",
                         data=data_zone,
                         iter=1000,
                         chains=3,
                         core=3,
                         include=FALSE,
                         pars=c("beta_raw"))
      save(fit.ancova,file=paste0(folder,"ancova_",z,"_",g,".RData"))
    }
    if(file.exists(paste0(folder,"ancova_",z,"_",g,"_lin.RData"))){
      load(paste0(folder,"ancova_",z,"_",g,"_lin.RData"))
    }else{
      X=model.matrix(~margin.temp*dh, sp.mod)
      data_zone=list(N=dim(sp.mod)[1],
                     NX=ncol(X),
                     S=nlevels(as.factor(sp.mod$species)),
                     species=as.numeric(as.factor(sp.mod$species)),
                     ISP=sp.mod$ISP,
                     X=X,
                     NULL)
      fit.ancova_lin <- stan(file = "lmm_ancova.stan",
                             data=data_zone,
                             iter=1000,
                             chains=3,
                             core=3,
                             include=FALSE,
                             pars=c("beta_raw"))
      save(fit.ancova_lin,file=paste0(folder,"ancova_",z,"_",g,"_lin.RData"))    
    }
    
    bic_quad=-2*summary(fit.ancova)$summary["lp__","mean"]+8*dim(sp.mod)[1]
    bic_lin=-2*summary(fit.ancova_lin)$summary["lp__","mean"]+5*dim(sp.mod)[1]
    
    if(bic_quad<bic_lin){
      fit=fit.ancova
      X=model.matrix(~margin.temp*dh+margin.temp*dh2, sp.mod)
    }else{
      fit=fit.ancova_lin
      X=model.matrix(~margin.temp*dh, sp.mod)
    }
    
    posteriors_fec<-as.data.frame(fit.ancova) |>
      dplyr::select(!matches("beta_")) |>
      dplyr::select(matches("beta"))
    colnames(posteriors_fec)[1:ncol(X)]=colnames(X)
    speciescode=sp.mod |> ungroup()|> dplyr::select(species) |> mutate(code=as.numeric(as.factor(species))) |> unique()
    speciesranef=as.data.frame(summary(fit.ancova)$summary) |> 
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
    
    output=cbind(sp.mod,
                 fit=fit,
                 quant05=quant05,
                 quant95=quant95) |>
      left_join(speciesranef) 
    
    out_ancova=rbind(out_ancova,
                     output)
    rm(output,quant05,quant95,fit,contrast_post,contrast_mat,speciesranef,speciescode,posteriors_fec,data_zone,sp.mod,sp.zone)
  }
}

out_ancova |> 
  ggplot()+
  geom_line(aes(dh_old,fit,color=margin.temp),size=0.6)+
  geom_point(aes(dh_old,log(ISP)-mean,color=margin.temp),size=0.6,alpha=0.6)+
  geom_ribbon(aes(x=dh_old,ymin=quant05,ymax=quant95,fill=margin.temp),alpha=0.8)+
  scale_fill_manual(values=c("bisque2","steelblue1","firebrick1"))+
  scale_color_manual(values=c("bisque2","steelblue1","firebrick1"))+
  facet_grid(genus~zone)
out_ancova |> 
  mutate(margin.temp=factor(margin.temp,levels=c("cold","midtemp","hot"))) |> 
  ggplot()+
  geom_line(aes(dh_old,fit,color=zone),size=0.6)+
  geom_point(aes(dh_old,log(ISP)-mean,color=zone),size=0.6,alpha=0.1)+
  geom_ribbon(aes(x=dh_old,ymin=quant05,ymax=quant95,fill=zone),alpha=0.8)+
  geom_boxplot(aes(dh_old,y=-10,color=zone),position=position_dodge(width=1),outlier.colour = NA)+
  # geom_boxplot(aes(x=-2000,y=log(ISP)-mean,color=zone),position=position_dodge(width=100))+
  scale_fill_manual(values=c("firebrick1","slateblue","firebrick4"))+
  scale_color_manual(values=c("firebrick1","slateblue","firebrick4"))+
  theme(axis.text.x = element_text(colour = inferno(6)))+
  facet_grid(genus~margin.temp)
