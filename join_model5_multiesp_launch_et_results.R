## Mise en forma donnee pour modele Stan
## Modele joint de croissance et de mortalite

###Résumé
## Etape A : simplification optionnelle des donnée
#      - suppression des recrut
#      - choix espèce ou sites 
## Etape B : calculs et ajouts des colonnes pour des effets aléatoires
## Etape C : lancement d'une chaine avec un dataj
## Etape D : mise en forme des sorties

rm(list=ls(all=TRUE))
gc() # garbage collector


#### chargement des packages et données #####

Library <- function(Packages) {
  InstallAndLoad <- function(Package) {
    if (!Package %in% installed.packages()[, 1]) {install.packages(Package)}
    require(Package, character.only = TRUE)
  }
  invisible(sapply(Packages, InstallAndLoad))
}

# Ajouter les packages necessaires ici
Library(c("knitr","raster","tidyverse",
          "ggplot2","rstan","bayesplot","rstanarm","dplyr"
          ))
# attention à l'ordre d'installation des packages, conflit de fonctions entre raster et tidyverse


# load("prepa_donnees/Data16esp_s0_to_join.Rdata") # issu de Donnees_guyafor-pre-traitement.Rmd
load("prepa_donnees/Data_guyafor_mch_join.Rdata")
# metadata 
     #nbmin_arbres,            # nombre minimal d'arbres par esp et par site
     #nbmin_mes,               # nombre minimal de mesures par arbres
     #an_min_meteo,            # donnee meteo les plus anciennes
     #tabMesuresSelect,        # table de données
     #tab_arbres_select,       # Id et statu des arbres sélectionnés
     #tab_esp,                 # effectif de d'arbre et de mesure par esp et par traitement 
     #tab_esp_site_logg,       # effectif d'arbres et mesure par esp, par site et par traitement
     #tab_site_logg,           # effectif d'arbres et mesures par site, par traitement et par espèce
     #esp_range_dmax,          # Min et max des dmax95 par espèce
     #Dmax_esp_site,           # dmax par espèce et par site
     #cc_espGFClim,            # correspondance espèce et code-espèce
     #ParacouInvIsh,           # valeur indice stress hydrique par forest et par int. d'inventaire
     #index_trait,              # valeur indice exploitation par site
     #datacr_cl_lo,            # données pour le modèle de croissance
     #datamo_cl_lo,            # données pour le modèlede mortalité
     #file=fich_sauv_circ) #"Data16esp_s0_to_join.Rdata") #"S0" pour seuil à 0 pour le nb minimal d'arbres par esp et par site, le 12/06/2019, avec le % de biomasse perdue sur IndexLogg





rstan_options(auto_write = T) # creer un fichier rds et evite donc l'etape de compilation si le .stan n'a pas etait modifie
options(mc.cores = parallel::detectCores()) # fixe le nombre de noyaux utilises. ici la fonction detectCores de la library parallel, renvoie le nombre total de coeur de l'ordi
#options(mc.cores = 4) # sur un serveur du reseau pour ne pas prendre tous les coeurs


#### A- simplification optionnelle des donnees, centrage et réduction des variables #####

##1 paramètres
#defaut
code_esp_cible<-"multiesp"
espliste<-espGFClim
parametres<-c(FALSE,
              NA,
              FALSE,
              FALSE,
              code_esp_cible) # valeur par défaut
names(parametres)<-c("selparacou",
                     "pcpara",
                     "prelogg",
                     "chablis",
                     "codeEsp")

selparacou<-TRUE # selection de quelques données de Paracou 
pcpara<-30 #pourcentage final de donnée de paracou

prelogg<-TRUE #exclusion de données pre-exploitation
# sellogg<-TRUE # exclusion des dispositifs sans exploitation

chablis<-TRUE # exclusion des chablis primaires et secondaires 

code_sorties<-"paracou30ncl_mch"

datacr_s<-datacr_cl_lo
datamo_s<-datamo_cl_lo

"avant selection"
nbl_cr<-nrow(datacr_s)
nbl_mo<-nrow(datamo_s)
nbl_cr
nbl_mo
# tableau résumé espèces
tabtemp_cr<-datacr_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) 

tabtemp<-datamo_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
  left_join(tabtemp_cr,by="idEsp")


##2 selection des donnees
 
if (prelogg) {
  parametres["prelogg"]<-prelogg
  datacr_s<-datacr_s %>% 
    left_join(dplyr::select(index_trait,Forest,Plot,Post_CY),by=c("Forest","Plot")) %>% 
    filter(is.na(Post_CY)|(Year1>=Post_CY))
  datamo_s<-datamo_s %>% 
    left_join(dplyr::select(index_trait,Forest,Plot,Post_CY),by=c("Forest","Plot")) %>% 
    filter(is.na(Post_CY)|(Year1>=Post_CY))
  code_esp_cible<-code_sorties
}
"perte après prelogg"
nrow(datacr_s)-nbl_cr
nrow(datamo_s)-nbl_mo
nbl_cr<-nrow(datacr_s)
nbl_mo<-nrow(datamo_s)
nbl_mo_mort<-nrow(datamo_s[datamo_s$to_Death==1,])
nbl_cr
nbl_mo
nbl_mo_mort

if (chablis) {
  parametres["chablis"]<-chablis
  datamo_s<-datamo_s %>% 
    filter(!((to_Death==1)&((CodeMeas3==6)|(CodeMeas3==7))))
  code_esp_cible<-code_sorties
}
"perte après chablis"
nrow(datamo_s)-nbl_mo
nbl_mo<-nrow(datamo_s)
nbl_mo_mort<-nrow(datamo_s[datamo_s$to_Death==1,])
nbl_mo
nbl_mo_mort

if (selparacou) { 
  parametres["selparacou"]<-selparacou
  parametres["pcpara"]<-pcpara 
  code_esp_cible<-code_sorties
  #calcul des nb de mesure totaux et par arbres sur les autres sites
  rtreemes_cr<-datacr_s%>% # calcul du nombre d'arbre, nb d'accroissement et d'accroissement par arbre sur paracou
    filter(Forest!="Paracou") %>% 
    group_by(idTree) %>% 
    summarise(nbmes=n()) %>% 
    summarise(nbtree_cr=n(),nb_cr=sum(nbmes)) %>% 
    mutate(rmes_cr=round(nb_cr/nbtree_cr))
  
  rtreemes<-datamo_s%>%
    filter(Forest!="Paracou") %>% 
    group_by(idTree) %>% 
    summarise(nbmes=n()) %>% 
    summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
    mutate(rmes_mo=round(nb_mo/nbtree_mo)) %>% 
    bind_cols(rtreemes_cr) 
  
  if (rtreemes$nb_cr/nrow(datacr_s)<((100-pcpara)/100)){ # test proportion d'accroissement de paracou supérieur à la proportion souhaitee
 
  ##croissance  
  #calcul du nombre d'individus de paracou nécessaire pour avoir la proportion souhaitee
  paratree_cr<-round((rtreemes$nb_cr)/rtreemes$rmes_cr*pcpara/(100-pcpara))
  parames_cr<-paratree_cr*rtreemes$rmes_cr
  
  #selection des arbres de paracou 
  Treeselpara_cr<-datacr_s %>%
    filter(Forest=="Paracou") %>% 
    group_by(idTree) %>% 
    summarise(nb_mes=n()) %>% 
    dplyr::select(idTree) 
 
  Treeselpara_cr<-Treeselpara_cr %>% 
    sample_n(min(paratree_cr,nrow(Treeselpara_cr)),replace=FALSE) # choix par les arbres de Paracou sans remise

  datacr_par<-datacr_s %>% 
    semi_join(Treeselpara_cr,by="idTree") %>% 
    sample_n(parames_cr,replace = FALSE) # choix parmis toutes les mesures, 
  
  # ajout des donnée de Paracou
  datacr_s<-datacr_s %>% 
    filter(Forest!="Paracou") %>% 
    bind_rows(datacr_par) 
   
   
     
  ##mortalité  
  #calcul du nombre d'individus de paracou nécessaire pour avoir la proportion souhaitee
  paratree_mo<-round(rtreemes$nb_mo/rtreemes$rmes_mo*pcpara/(100-pcpara))
  parames_mo<-paratree_mo*rtreemes$rmes_mo
  

  # recensement les arbres qui meurent effectivement sur Paracou
  tree_to_death<-datamo_s %>%
    filter((Forest=="Paracou") & (to_Death==1))

  #selection des arbres de paracou 
  Treeselpara_mo<-datamo_s %>%
    filter(Forest=="Paracou") %>% 
    group_by(idTree) %>% 
    summarise(nb_mes=n()) %>% 
    dplyr::select(idTree) %>%            
    anti_join(tree_to_death,by="idTree") # on enlève les acc avec des arbres qui meurent 
  

  paratree_mo<-paratree_mo-nrow(tree_to_death) # recalcul du nombre d'arbres nécessaire
  parames_mo<-parames_mo-(rtreemes$rmes_mo*nrow(tree_to_death)) # recalcul du nombre d'accroissement nécessaire
  # on ne prends que le nombre moyens d'acroissement par arbres car si on prends tous les acc des arbres qui meurent on peut dépasser le nombre de mesures necessaire

  # choix des arbres de paracou 
  # choix des acc non mortels parmi les arbres qui meurent
  datamo_par_mort<-datamo_s %>% 
    semi_join(tree_to_death,by="idTree") %>% 
    filter(to_Death==0) %>%  # acc des arbres qui meurent sans les acc mortels
    sample_n((rtreemes$rmes_mo-1)*nrow(tree_to_death),replace = FALSE) %>%  #nb choix d'acc vaut rmes_mo en comptant les acc mortel
    bind_rows(tree_to_death) # ajout acc mortels
  
  Treeselpara_mo<-Treeselpara_mo %>% 
    sample_n(paratree_mo, replace = FALSE)  # choix du reste des arbres necessaires

  datamo_par<-datamo_s %>% 
    semi_join(Treeselpara_mo,by="idTree") %>% 
    sample_n(parames_mo,replace=FALSE) %>% #choix des mesures
    bind_rows(datamo_par_mort) # ajout des arbres qui meurent et de leurs mesures
 

   test<-datamo_par %>%
     group_by(idTree) %>%
     summarise(nb_acc=n())
   median(test$nb_acc)
   sum(test$nb_acc)
 
  datamo_s<-datamo_s %>% 
    filter(Forest!="Paracou") %>% 
    bind_rows(datamo_par) 
  }
}
"perte après selparacou"
nrow(datacr_s)-nbl_cr
nrow(datamo_s)-nbl_mo
nrow(datamo_s[datamo_s$to_Death==1,])-nbl_mo_mort

nbl_cr<-nrow(datacr_s)
nbl_mo<-nrow(datamo_s)
nbl_mo_mort<-nrow(datamo_s[datamo_s$to_Death==1,])
nbl_cr
nbl_mo
nbl_mo_mort


##3 Centrage et réduction des variables Logg et climat
CentreReduit <- function(vec)
  return(vec-mean(vec,na.rm=T))/sd(vec,na.rm=T)

datacr_s$IshInv<-CentreReduit(datacr_s$IshInv)
datamo_s$IshInvMo<-CentreReduit(datamo_s$IshInvMo)

datamo_s$TxLogg<-CentreReduit(datamo_s$TxLogg)
datacr_s$TxLogg<-CentreReduit(datacr_s$TxLogg)


##4 graphes et tableau résumé 

# tableau résumé espèces
tabspe_cr<-datacr_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) %>% 
  arrange(desc(nb_acc))

tab_wd<-datacr_s %>% 
  group_by(idEsp) %>% 
  summarise(WDmoy=mean(WD))

tabspe<-datamo_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
  left_join(tabspe_cr,by="idEsp")

ggplot(tabspe, aes(x=reorder(idEsp,nb_acc),y=nb_acc))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.7))

# tableau résumé sites

tabforests_cr<-datacr_s%>%
  group_by(Forest,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(Forest) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) %>% 
  arrange(desc(nb_acc))

tabforests<-datamo_s%>%
  group_by(Forest,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(Forest) %>% 
  summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
  left_join(tabforests_cr,by="Forest")

ggplot(tabforests, aes(x=reorder(Forest,nb_acc),y=nb_acc))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.7))

ggplot(tabforests, aes(x=reorder(Forest,nb_mo),y=nb_mo))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.7))

"ratio données Paracou"
"croissance" 
tabforests$nb_acc[tabforests$Forest=="Paracou"]/sum(tabforests$nb_acc)
"mortalité" 
tabforests$nb_mo[tabforests$Forest=="Paracou"]/sum(tabforests$nb_mo)

# distribution des variables

par(mfrow=c(2,3))
 hist(datamo_cl_lo$IshInvMo)
 hist(datamo_s$IshInvMo)
 hist(datamo_cl_lo$IshInvMo[datamo_cl_lo$Forest=="Paracou"])
 hist(datamo_cl_lo$TxLogg)
 hist(datamo_s$TxLogg)
 hist(datamo_cl_lo$TxLogg[datamo_cl_lo$Forest=="Paracou"])


 par(mfrow=c(2,3))
 hist(datamo_cl_lo$Diam1)
 hist(datamo_s$Diam1)
 hist(datamo_cl_lo$Diam1[datamo_cl_lo$Forest=="Paracou"])
 hist(datacr_cl_lo$Diam1)
 hist(datacr_s$Diam1)
 hist(datacr_cl_lo$Diam1[datacr_cl_lo$Forest=="Paracou"])

 
code_esp_cible

#### B- Calculs et ajouts des colonnes pour des effets aléatoires ####

#1 construction du vecteur pour effet aleatoire
    #sur individus
  larbres<-sort(unique(c(unique(datacr_s$idTree),unique(datamo_s$idTree)))) # liste complete des arbres
  tarbres<-data.frame(idTree=larbres,rank=c(1:length(larbres)))

  datacr <- datacr_s%>%
    left_join(tarbres,by =c("idTree")) %>%  # on ajoute le rang de chaque arbres dans la liste complete des arbres
    rename(TreeRank=rank)
  datamo<- datamo_s%>%
    left_join(tarbres,by =c("idTree")) %>%  # on ajoute le rang de chaque arbres dans la liste complete des arbres
    rename(TreeRank=rank)

    # sur especes  
  tesp<-data.frame(idEsp=as.character(espliste),rank=c(1:length(espliste)))
  
  datacr <- datacr%>%
    left_join(tesp,by =c("idEsp")) %>%   # on ajoute le rang de chaque espece dans la liste complete des arbres
    rename(EspRank=rank)
  
  datamo <- datamo%>%
    left_join(tesp,by =c("idEsp")) %>%   # on ajoute le rang de chaque espece dans la liste complete des arbres
    rename(EspRank=rank)
  
 

  #3 construction de la liste pour Stan
  dataj <- list(
    Ncr = nrow(datacr),           # Sample size pop1 : accroissement observes
    Nmo = nrow(datamo),         # Sample size pop2 evenements de survie et de mort
    Nt = length(larbres),           # nombre d'arbre distincts
    Nesp = length(espliste),           # nombre d'arbre distincts
    Acc_cr = datacr$AGR,           # accroissement entre n-1 et n
    Acc_mo = datamo$AGR,        # accroissement entre n-2 et n-1 suivi d'un evenement de survie/mort
    clim_cr = datacr$IshInv,        # covariable climat n-1 modèle de croissance
    clim_vig = datamo$IshInvVig,        # covariable climat entre n-2 et n-3 calcul pred accroissement avant un evenement de survie/mort
    clim_mo = datamo$IshInvMo,      # covariable climat n-1 calcul logit d'un evenement de survie/mort
    logg_cr = datacr$TxLogg,        # covariable exploitation modèle de croissance
    logg_mo = datamo$TxLogg,        # covariable exploitation modèle de mortalité
    dbh_cr =  datacr$Diam1,        # covariable diametre initial modele croissance
    dbh_mo = datamo$Diam2,        # covariable diametre n-1 calcul proba mort avant un evenement de survie/mort
    dbh_vig = datamo$Diam1,        # covariable diametre n-2 calcul pred accroissement  avant un evenement de survie/mort
    dbhmax_mo = datamo$Dmax95,        # dmax pour modèle de mortalité
    dbhmax_cr = datacr$Dmax95,        # dmax pour modèle de croissance
    WD_mo = datamo$WD,        # Wood density pour modèle de mortalité
    WD_cr = datacr$WD,        # Wood density pour modèle de croissance
    morts = datamo$to_Death,     # vecteur de 1 = evenement de mort et de 0 = evenements de survie
    rgtree_cr = datacr$TreeRank ,        # rang des individu dans la liste  arbres pour croissance
    rgtree_mo = datamo$TreeRank,       # rang des individu dans la liste arbres pour la mortalite
    rgEsp_cr = datacr$EspRank,        # rang des espece dans la liste des espece pour croissance
    rgEsp_mo = datamo$EspRank        # rang des espece dans la liste des espece pour la mortalité
  )

  save(dataj,parametres,tesp,code_esp_cible,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data.Rdata'))
  save(dataj,parametres,tesp,code_esp_cible,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data2.Rdata'),version=2)  # pour export vers Rstudio serveur
  
#### C- Lancement des chaines et étude de convergence ####

# modèle complet ####
pars_save<-c("Ks","Dopt","sigma",
             "cr_clim","cr_Clesp","cr_sigClesp",
             "cr_logg","cr_sigLoesp","cr_Loesp",
             "cr_int","cr_Inesp","cr_sigInesp",
             "oo_Gmax","cr_Gesp","cr_dmax","cr_sigGesp",
             "oo_logit","mo_Ooesp","vig","onto","onto_sq","mo_sigOoesp",
             "mo_clim","mo_Clesp","mo_sigClesp",
             "mo_logg","mo_Loesp","mo_sigLoesp",
             "mo_int","mo_Inesp","mo_sigInesp")

mtype<-"joint"
  
temps_depart <-Sys.time()
temps_depart
fitjm <- stan('join_model5_multiesp.stan', data = dataj,pars=pars_save,include=TRUE,
               chain=2,
               iter=4000,warmup=3000
              )

Sys.time()- temps_depart
save(fitjm,parametres,pars_save,code_esp_cible,tesp,mtype,
     file=paste('stan_sorties/stan_',code_esp_cible,'_jo5_sortie_c2i2.Rdata',sep="")) # voir debut du fichier stan pour les specificite


chain<-fitjm
mtype<-"joint"
code_sim<-paste('stan_',code_esp_cible,'_jo5_c2i2',sep="")
# tesp<-tesp %>% transmute(especes=idEsp,nesp=rank)

# modèle de croissance seul ####

pars_save<-c("Ks","Dopt","sigma",
             # "Ks2","KsC","DoptC",
             "cr_clim","cr_Clesp","cr_sigClesp",
             "cr_logg","cr_sigLoesp","cr_Loesp",
             "oo_Gmax","cr_Gesp","cr_dmax","cr_sigGesp",
             "cr_int","cr_Inesp","cr_sigInesp")

mtype<-"croiss"

temps_depart <-Sys.time()
temps_depart
fitj_cr <- stan('cr_model5_multiesp.stan', data = dataj,pars=pars_save,include=TRUE,
                chain=2,
                iter=4000,warmup=3000
                )
Sys.time()- temps_depart
save(fitj_cr,parametres,pars_save,code_esp_cible,tesp,mtype,
     file=paste('stan_sorties/stan_',code_esp_cible,'_cr5_sortie_i4.Rdata'))

chain<-fitj_cr
mtype<-"croiss"
code_sim<-paste('stan_',code_esp_cible,'_cr5_i4')
tesp<-tesp %>% transmute(especes=idEsp,nesp=rank)

# modèle mortalité seul ####
pars_save<-c("oo_logit","mo_Ooesp","onto","onto_sq","mo_sigOoesp",
             "mo_clim","mo_Clesp","mo_sigClesp",
             "mo_logg","mo_Loesp","mo_sigLoesp",
             "mo_int","mo_Inesp","mo_sigInesp")

mtype<-"morta"

temps_depart <-Sys.time()
temps_depart
fitj_mo <- stan('mo_model5_multiesp.stan', data = dataj,pars=pars_save,include=TRUE,
                chain=4,
                iter=2000,warmup=1000,
               )
Sys.time()- temps_depart
save(fitj_mo,parametres,pars_save,code_esp_cible,tesp,mtype,
     file=paste('stan_sorties/stan_',code_esp_cible,'_mo5_ncl_sortie_i2.Rdata'))

chain<-fitj_mo
mtype<-"morta"
code_sim<-paste('stan_',code_esp_cible,'_mo5_i2')
tesp<-tesp %>% transmute(especes=idEsp,nesp=rank)


## calcul vecteurs de  sous ensembles de paramètres ####
cr_Gesp_multi<-paste("cr_Gesp[",1:nrow(tesp),"]",sep="")
cr_Clesp_multi<-paste("cr_Clesp[",1:nrow(tesp),"]",sep="")
cr_Loesp_multi<-paste("cr_Loesp[",1:nrow(tesp),"]",sep="")
cr_Inesp_multi<-paste("cr_Inesp[",1:nrow(tesp),"]",sep="")

mo_Ooesp_multi<-paste("mo_Ooesp[",1:nrow(tesp),"]",sep="")
mo_Clesp_multi<-paste("mo_Clesp[",1:nrow(tesp),"]",sep="")
mo_Loesp_multi<-paste("mo_Loesp[",1:nrow(tesp),"]",sep="")
mo_Inesp_multi<-paste("mo_Inesp[",1:nrow(tesp),"]",sep="")



pars_all<-chain@model_pars

pars_cr1<-c("oo_Gmax","Ks","Dopt","cr_dmax","cr_clim","cr_logg","cr_int","sigma","cr_sigGesp","cr_sigClesp","cr_sigLoesp")
pars_cr11<-c("oo_Gmax","Ks","Dopt","cr_dmax","cr_clim","cr_logg","sigma")
pars_crG<-c("oo_Gmax","cr_sigGesp",cr_Gesp_multi)
pars_crCl<-c("cr_clim","cr_sigClesp",cr_Clesp_multi)
pars_crLogg<-c("cr_logg","cr_sigLoesp",cr_Loesp_multi)
pars_crInt<-c("cr_int","cr_sigInesp",cr_Inesp_multi)
pars_cr_all<-unique(c(pars_cr1,pars_crG,pars_crCl,pars_crLogg,pars_crInt))



pars_mo1<-c("oo_logit","onto","onto_sq","mo_clim","mo_logg","mo_int","mo_sigOoesp","mo_sigClesp","mo_sigLoesp")
pars_mo11<-c("oo_logit","onto","onto_sq","mo_clim","mo_logg")
pars_moOo<-c("oo_logit","mo_sigOoesp",mo_Ooesp_multi)
pars_moCl<-c("mo_clim","mo_sigClesp",mo_Clesp_multi)
pars_moLogg<-c("mo_logg","mo_sigLoesp",mo_Loesp_multi)
pars_moInt<-c("mo_int","mo_sigInesp",mo_Inesp_multi)
if(mtype=="joint"){
  pars_mo1<-c(pars_mo1,"vig")
  pars_mo11<-c(pars_mo11,"vig")
}
pars_mo_all<-unique(c(pars_mo1,pars_moOo,pars_moCl,pars_moLogg,pars_moInt))

## 1 Rhat ####
print(chain, pars = pars_save) 

## 2 Traces des chaines ####
traceplot(chain, pars=c("lp__")) # plot les chaines de la vraisemblance
nrtg<-ifelse(mtype=="joint",12,6)
traceplot(chain, pars=pars_save, nrow=nrtg) # nrow : nombre de lignes de graphe
traceplot(chain,pars=pars_mo1)
## 3 calculs : création sortie en format df et pour ggplot - choix éventuel de chaines ####

chain_to_df<-function (chain_ini){
nbchain<-chain@sim$chains # extraction du nombre de chaines de sampling après warmup
nbsample<-chain@sim$iter-chain@sim$warmup
chain_out<-as.data.frame(chain) 
chain_out$chaines<-as.factor(rep(1:nbchain,each=nbsample)) # construction d'une colonnes iteration et numero de chaine
chain_out$iterations<-rep(1:nbsample,nbchain)
return(chain_out)
}

chain_df<-chain_to_df(chain) 

chain_sel<-chain_df
  # filter(chaines!=2) %>% 
  # filter (chaines!=4) 

chain_ggpl<-chain_sel %>% 
  pivot_longer(-c(iterations,chaines),names_to ="variables",values_to = "valeurs" ) %>% 
  mutate(test=ifelse(substr(variables,nchar(variables),nchar(variables))=="]",    # extraction numéro espèce du titre variables
                     ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",1,2)
                     ,0)) %>% 
  mutate(nesp=ifelse(test==1,as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     ifelse(test==2,as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)),
                            -1))) %>% 
  select(-test) %>% 
  left_join(tesp,by="nesp") %>% # ajout nom complet espèce 
  mutate(especes_var=especes, idtemp=paste(chaines,iterations,variables))

# Création colonne pour graphiques avec nom de l'expèces si la variables spécifique et celui de la varible si nonspé ou 
 # Remarque  : mutate(especes_var=ifelse(is.na(especes),variables,especes) affecte, si is.na=FALSE, la valeur du level de espèce 
 # au lieu d'affecter la chaine de caractère
temp<-chain_ggpl %>% 
  filter(is.na(especes)) %>% 
  mutate(especes_var=variables)

chain_ggpl<-chain_ggpl %>% 
  filter(!is.na(especes)) %>% 
  bind_rows(temp) %>% 
  arrange(chaines,iterations,nesp,variables)
  
save(chain_sel,chain_ggpl,code_esp_cible,parametres,tesp,pars_save,mtype,
      file=paste('stan_sorties/',code_sim,'_chain.Rdata',sep=""))

## 4 traces chaines choisies éventuellement ####
ggplot(chain_ggpl) + 
  geom_line(aes(x = iterations, y = valeurs ,color= chaines)) +
  xlab("itérations") +
  facet_wrap(~variables,scales="free_y")


## 5 nuage de points pour couple de paramètres donnés ####

mcmc_pairs(as.array(chain), pars = pars1)
mcmc_pairs(as.array(chain), pars = parsG)
mcmc_pairs(as.array(chain), pars = pars_crCl)
mcmc_pairs(as.array(chain), pars = parcLogg)
# GGally::ggpairs() pour faire des pairs plus lisibles

mcmc_pairs(as.array(chain), pars = pars_mo11)
mcmc_pairs(as.array(chain), pars = pars_mo1)
mcmc_pairs(as.array(chain), pars = pars_moOo)
mcmc_pairs(as.array(chain), pars = pars_moCl)
mcmc_pairs(as.array(chain), pars = pars_moLogg1)


## 6 Posteriors ####


poster_gg <-function (vectorpars,legende){
  ggplot(chain_ggpl %>% filter(variables%in%vectorpars)) + 
    geom_freqpoly(aes(x=valeurs,binwidth=0.015,color= chaines)) +
    # geom_freqpoly(aes(x=valeurs,binwidth=0.015),size=1,color="red") +
    xlab(legende) +
    geom_vline(xintercept = 0, linetype = "solid")+
    facet_wrap(~especes_var,scales="free")
}

# test<-chain_ggpl %>% 
#   filter(substr(variables,1,5)=="cr_In")

poster_gg(pars_cr1,"paramètres non spécifiques")
poster_gg(pars_crG,"paramètre de l'ordonnée à l'origine")
poster_gg(pars_crLogg,"paramètre exploitation")
poster_gg(pars_crCl,"paramètres climat")
poster_gg(pars_crInt,"paramètres interactions Climat et exploitation")

poster_gg(pars_mo1,"paramètres non spécifiques")
poster_gg(pars_moOo,"paramètre de l'ordonnée à l'origine")
poster_gg(pars_moLogg,"paramètre exploitation")
poster_gg(pars_moCl,"paramètres climat")
poster_gg(pars_moInt,"paramètres interaction climat et exploitation")

poster_gg(c(pars_cr1,pars_mo1),"paramètres non spécifiques modèle joint")


## 7 calculs Predictions et graphe de sortie ####
temp<-chain_sel %>%   
  select(-chaines,-iterations)
param<-apply(temp,2,median) # on extrait les médianes des paramètres
names(param)<-colnames(temp)
load(file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data.Rdata',sep=""))

###1 calcul des prédictions pour chaque observé ####

if (mtype %in% c("croiss","joint")){
 donnee_cr<-as.data.frame(dataj[][c("Acc_cr","clim_cr","logg_cr","dbh_cr","dbhmax_cr","WD_cr","rgEsp_cr")]) #ensemble des donnees liees au modele de croissance
 donnee_cr<-donnee_cr %>% left_join(transmute(tesp,especes=especes,rgEsp_cr=nesp),by="rgEsp_cr")
 Nesp_cr<-donnee_cr$rgEsp_cr
 
 donnee_cr$pred <- (param["oo_Gmax"]+ param[paste("cr_Gesp[",Nesp_cr,"]",sep="")] +
              param["cr_dmax"]*donnee_cr$dbhmax_cr+
              (param["cr_clim"]+param[paste("cr_Clesp[",Nesp_cr,"]",sep="")])* donnee_cr$clim_cr+
              (param["cr_logg"]+param[paste("cr_Loesp[",Nesp_cr,"]",sep="")])* donnee_cr$logg_cr)*
   exp((-0.5)*(log(donnee_cr$dbh_cr/(donnee_cr$dbhmax*param["Dopt"]))/(param["Ks"]*donnee_cr$WD_cr))*
             (log(donnee_cr$dbh_cr/(donnee_cr$dbhmax*param["Dopt"]))/(param["Ks"]*donnee_cr$WD_cr)))

ggplot(donnee_cr) +
  geom_point(aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1))) +
  geom_line (aes(x=dbh_cr/dbhmax_cr,y=pred),color="red",size=1) +
  facet_wrap(~especes)+
  labs(title="Modèle de croissance",x="diam/dmax",y="Log(Acc+1)")
}

if (mtype=="morta"){
 donnee_mo<-as.data.frame(dataj[][c("Acc_mo","clim_vig","clim_mo","logg_mo","dbh_mo","dbhmax_mo","WD_mo","morts","rgEsp_mo")]) #donnees du modele de mortalité
 donnee_mo<-donnee_mo %>% left_join(transmute(tesp,especes=especes,rgEsp_mo=nesp),by="rgEsp_mo")
 Nesp_mo<-donnee_mo$rgEsp_mo
 
 donnee_mo$pred<-(param["oo_logit"]+ param[paste("mo_Ooesp[",Nesp_mo,"]",sep="")] +
                     param["onto"]*donnee_mo$dbh_mo/donnee_mo$dbhmax_mo+
                     param["onto_sq"]*(donnee_mo$dbh_mo/donnee_mo$dbhmax_mo)*(donnee_mo$dbh_mo/donnee_mo$dbhmax_mo)+
                     (param["mo_clim"]+param[paste("mo_Clesp[",Nesp_mo,"]",sep="")])*donnee_mo$clim_mo+
                     (param["mo_logg"]+param[paste("mo_Loesp[",Nesp_mo,"]",sep="")])*donnee_mo$logg_mo)+
                     (param["mo_int"]+param[paste("mo_Inesp[",Nesp_mo,"]",sep="")])*donnee_mo$logg_mo*donnee_mo$clim_mo

 donnee_mo_logit<-donnee_mo %>%
   mutate(logitm1=exp(pred)/(1+exp(pred))) %>% # on calcule le logit-1(pred) c'est à dire la probabilité de mort prévu pour chaque individus
   arrange(morts,logitm1) %>% 
   mutate(mort_ch=ifelse(morts==1,"arbres morts","arbres vivants")) %>% #légende de futurs graphes
   group_by(morts) %>% 
   mutate(numlig=row_number())

 logitm1_pc<-data.frame(matrix(data=NA,nrow=50,ncol=3))
 colnames(logitm1_pc)<-c("quantile","morts","vivants")
 logitm1_pc$quantile<-seq(0.02,1,by=0.02)

 for(p in 1:nrow(logitm1_pc)){
  logitm1_pc$morts[p]<-quantile(donnee_mo_logit$logitm1[donnee_mo_logit$morts==1],
                    probs=logitm1_pc$quantile[p])
  logitm1_pc$vivants[p]<-quantile(donnee_mo_logit$logitm1[donnee_mo_logit$morts==0],
                    probs=logitm1_pc$quantile[p])
  
 }
  
 ggplot(donnee_mo_logit) +
  geom_line (aes(x=numlig,y=logitm1,color=mort_ch),size=1) +
  facet_wrap(~mort_ch,scales="free_x")+
  labs(title="Modèle de mortalité",x="individu",y="proba mort")
 
 ggplot(logitm1_pc) +
  geom_line (aes(x=quantile,y=morts,color="red"),size=1) +
  geom_line (aes(x=quantile,y=vivants,color="blue"),size=1)+
  labs(title="Modèle de mortalité",x="quantile population",y="proba mort")
 
}

if (mtype=="joint"){
  donnee_jo<-as.data.frame(dataj[][c("Acc_mo","clim_vig","clim_mo","logg_mo","dbh_mo","dbhmax_mo","WD_mo","morts","rgEsp_mo")]) #donnees du modele de mortalité
  donnee_jo<-donnee_jo %>% left_join(transmute(tesp,especes=especes,rgEsp_mo=nesp),by="rgEsp_mo")
  Nesp_mo<-donnee_jo$rgEsp_mo
  
  donnee_jo$logVig<-(param["oo_Gmax"]+ param[paste("cr_Gesp[",Nesp_mo,"]",sep="")] +
                     param["cr_dmax"]*donnee_jo$dbhmax_mo+
                    (param["cr_clim"]+param[paste("cr_Clesp[",Nesp_mo,"]",sep="")])* donnee_jo$clim_vig+
                    (param["cr_logg"]+param[paste("cr_Loesp[",Nesp_mo,"]",sep="")])* donnee_jo$logg_mo+
                    (param["cr_int"]+param[paste("cr_Inesp[",Nesp_mo,"]",sep="")])*
                                                                 donnee_jo$logg_mo*donnee_jo$clim_vig )*
         exp((-0.5)*(log(donnee_jo$dbh_mo/(donnee_jo$dbhmax_mo*param["Dopt"]))/(param["Ks"]*donnee_jo$WD_mo))*
                    (log(donnee_jo$dbh_mo/(donnee_jo$dbhmax_mo*param["Dopt"]))/(param["Ks"]*donnee_jo$WD_mo)))
  
  donnee_jo$pred <- (param["oo_logit"]+ param[paste("mo_Ooesp[",Nesp_mo,"]",sep="")] +
                     param["vig"]*(log(donnee_jo$Acc_mo+1)-donnee_jo$logVig))+   
                     param["onto"]*donnee_jo$dbh_mo/donnee_jo$dbhmax+
                     param["onto_sq"]*(donnee_jo$dbh_mo/donnee_jo$dbhmax)*(donnee_jo$dbh_mo/donnee_jo$dbhmax)+
                     (param["mo_clim"]+param[paste("mo_Clesp[",Nesp_mo,"]",sep="")])*donnee_jo$clim_mo+
                     (param["mo_logg"]+param[paste("mo_Loesp[",Nesp_mo,"]",sep="")])*donnee_jo$logg_mo+
   (param["mo_int"]+param[paste("mo_Inesp[",Nesp_mo,"]",sep="")])*donnee_jo$logg_mo*donnee_jo$clim_mo

 donnee_jo_logit<-donnee_jo %>%
  mutate(logitm1=exp(pred)/(1+exp(pred))) %>% # on calcule le logit-1(pred) c'est à dire la probabilité de mort prévu pour chaque individus
  arrange(morts,logitm1) %>% 
  mutate(mort_ch=ifelse(morts==1,"arbres morts","arbres vivants")) %>% #légende de futurs graphes
  group_by(morts) %>% 
  mutate(numlig=row_number()) # numero de ligne par type mort ou vivant

 #graphe
 logitm1_pc<-data.frame(matrix(data=NA,nrow=50,ncol=3)) # construction des abscisses
 colnames(logitm1_pc)<-c("quantile","morts","vivants")
 logitm1_pc$quantile<-seq(0.02,1,by=0.02)
 
  for(p in 1:nrow(logitm1_pc)){
  logitm1_pc$morts[p]<-quantile(donnee_jo_logit$logitm1[donnee_jo_logit$morts==1],
                                probs=logitm1_pc$quantile[p])
  logitm1_pc$vivants[p]<-quantile(donnee_jo_logit$logitm1[donnee_jo_logit$morts==0],
                                  probs=logitm1_pc$quantile[p])
  
 }

 ggplot(donnee_jo_logit) +
  geom_line (aes(x=numlig,y=logitm1,color=mort_ch),size=1) +
  facet_wrap(~mort_ch,scales="free_x")+
   labs(title="Modèle de mortalité",x="individu",y="proba mort")
 

 ggplot(logitm1_pc) +
  geom_line (aes(x=quantile,y=morts,color="red"),size=1) +
  geom_line (aes(x=quantile,y=vivants,color="blue"),size=1) +
   labs(title="Modèle de mortalité",x="quantile population",y="proba mort")
 
}

###2 calcul de prévision et erreur pour diamètre (ontogénie) donné par espèce et sur un range de variables climatiques ####

## Fonctions de prédictions
# pour accélérer la boucle, stocker 1000 tirages de rnorm dans un vecteur 
# et ensuite aller chercher les valeurs dans ce vecteur

Fpred_cr<-function(i,vparamj,ontos,datamed,sigmas){
  return (as.numeric(vparamj["oo_Gmax"])+as.numeric(sigmas[i,]["cr_sigGesp"])+
            as.numeric(vparamj["cr_dmax"])*datamed["dbhmax_cr"]+
                (as.numeric(vparamj["cr_clim"])+as.numeric(sigmas[i,]["cr_sigClesp"]))*datamed["clim_cr"]+
                (as.numeric(vparamj["cr_logg"])+as.numeric(sigmas[i,]["cr_sigLoesp"]))*datamed["logg_cr"])*
    exp((-0.5)*(log(ontos[i]/as.numeric(vparamj["Dopt"]))/(as.numeric(vparamj["Ks"])*datamed["WD_cr"]))*
               (log(ontos[i]/as.numeric(vparamj["Dopt"]))/(as.numeric(vparamj["Ks"])*datamed["WD_cr"])))+
    sigmas[i,]["sigma"]
}

Fpred_cr_spe<-function(i,vparamj,ontos,datamed,dataspe,sigmas,nesp){
  #vparamj : vecteur de paramètres ( longueur =100*nb chaine =40000)
  #ontos : abscisse du graph : ratio diam/dmax = ontogénie
  #datamed : mediane des varaibles explicatives non dpdte de l'espèce : clim, logg
  #dataspe : varaible explicatie dépendante de l'espèce : dbhmax, WD
  #sigmas : tableau de valeur tiré dans les loi normale des effet aléatoir et de l'erreiur du modèle
  #nesp : numéro d'espèce (rang dans la liste des espèces)
  return (as.numeric(vparamj["oo_Gmax"])+
            vparamj[paste("cr_Gesp[",nesp,"]",sep="")]+
            as.numeric(sigmas[i,]["cr_sigGesp"])+
            as.numeric(vparamj["cr_dmax"])*dataspe[dataspe$rgEsp_cr==nesp,"dbhmaxspe_cr"]+
           # (as.numeric(vparamj["cr_logg"])+
           #       vparamj[paste("cr_Loesp[",nesp,"]",sep="")]+
           #       as.numeric(sigmas[i,]["cr_sigLoesp"]))*datamed["logg_cr"]+
           (as.numeric(vparamj["cr_clim"])+
                 vparamj[paste("cr_Clesp[",nesp,"]",sep="")]+
                 as.numeric(sigmas[i,]["cr_sigClesp"]))*datamed["clim_cr"]
          )*
          exp((-0.5)*(log(ontos[i]/as.numeric(vparamj["Dopt"]))/
                        (as.numeric(vparamj["Ks"])*dataspe[dataspe$rgEsp_cr==nesp,"WDspe_cr"]))*
                     (log(ontos[i]/as.numeric(vparamj["Dopt"]))/
                         (as.numeric(vparamj["Ks"])*dataspe[dataspe$rgEsp_cr==nesp,"WDspe_cr"]))
              )+
          sigmas[i,]["sigma"]
}

Fpred_cr_spe_clim<-function(i,vparamj,ontos,dataclim,datamed,dataspe,sigmas,nesp){
  #vparamj : vecteur de paramètres ( longueur =100*nb chaine =40000)
  #ontos : abscisse du graph : ratio diam/dmax = ontogénie
  #dataclim : variable explicative ciblee : Clim
  #datamed : mediane des variables explicatives non ciblees : logg
  #dataspe : variable explicatie dépendante de l'espèce : dbhmax, WD
  #sigmas : tableau de valeur tiré dans les loi normale des effet aléatoire et de l'erreur du modèle
  #nesp : numéro d'espèce (rang dans la liste des espèces)
  return (as.numeric(vparamj["oo_Gmax"])+
            vparamj[paste("cr_Gesp[",nesp,"]",sep="")]+
            as.numeric(sigmas[i,]["cr_sigGesp"])+
            as.numeric(vparamj["cr_dmax"])*dataspe[dataspe$rgEsp_cr==nesp,"dbhmaxspe_cr"]+
            # (as.numeric(vparamj["cr_logg"])+
            #       vparamj[paste("cr_Loesp[",nesp,"]",sep="")]+
            #       as.numeric(sigmas[i,]["cr_sigLoesp"]))*datamed["logg_cr"]+
            (as.numeric(vparamj["cr_clim"])+
               vparamj[paste("cr_Clesp[",nesp,"]",sep="")]+
               as.numeric(sigmas[i,]["cr_sigClesp"]))*dataclim[i]
  )*
    exp((-0.5)*(log(ontos[i]/as.numeric(vparamj["Dopt"]))/
                  (as.numeric(vparamj["Ks"])*dataspe[dataspe$rgEsp_cr==nesp,"WDspe_cr"]))*
               (log(ontos[i]/as.numeric(vparamj["Dopt"]))/
                  (as.numeric(vparamj["Ks"])*dataspe[dataspe$rgEsp_cr==nesp,"WDspe_cr"]))
    )+
    sigmas[i,]["sigma"]
}


Fpred_mo_spe_clim<-function(i,vparamj,diam,mo_dataClim,mo_datamed,mo_dmaxspe,nesp){
  #vparamj : vecteur de paramètres ( longueur =100*nb chaine =40000)
  #diam : abscisse du graph : ratio diam/dmax = ontogénie
  #mo_dataClim : variable explicative ciblee : Clim
  #mo_datamed : mediane des variables explicatives non ciblees : logg
  #modmaxspe : variable explicatie dépendante de l'espèce : dbhmax
  #nesp : numéro d'espèce (rang dans la liste des espèces)
  return  (as.numeric(vparamj["oo_logit"])+
             vparamj[paste("mo_Ooesp[",nesp,"]",sep="")]+
             as.numeric(vparamj["onto"])*diam[i]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]+
             as.numeric(vparamj["onto_sq"])*diam[i]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]*
             diam[i]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]+
             (as.numeric(vparamj["mo_logg"])+vparamj[paste("mo_Loesp[",nesp,"]",sep="")])*
                    mo_datamed["logg_mo"]+
             (as.numeric(vparamj["mo_clim"])+vparamj[paste("mo_Clesp[",nesp,"]",sep="")])*
                    mo_dataClim[i]+
             (as.numeric(vparamj["mo_int"])+vparamj[paste("mo_Inesp[",nesp,"]",sep="")])*
                    mo_dataClim[i]*mo_datamed["logg_mo"])
}


Fpred_mo_inter_spe<-function(i,param,mo_dataClim,mo_dataLogg,mo_datamed,mo_dmaxspe,nesp){
  return  (param["oo_logit"]+
             param[paste("mo_Ooesp[",nesp,"]",sep="")]+
             param["onto"]*mo_datamed["dbh_mo"]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]+
             param["onto_sq"]*mo_datamed["dbh_mo"]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]*
                              mo_datamed["dbh_mo"]/mo_dmaxspe[mo_dmaxspe$rgEsp_mo==nesp,"dbhmaxspe_mo"]+
             (param["mo_logg"]+param[paste("mo_Loesp[",nesp,"]",sep="")])*
             mo_dataLogg[i]+
             (param["mo_clim"]+param[paste("mo_Clesp[",nesp,"]",sep="")])*
             mo_dataClim[i]+
             (param["mo_int"]+param[paste("mo_Inesp[",nesp,"]",sep="")])*
             mo_dataClim[i]*mo_dataLogg[i])
  
}

Fpred_mo_inter<-function(i,param,mo_dataClim,mo_dataLogg,mo_datamed){
  return  (param["oo_logit"]+
           param["onto"]*mo_datamed["dbh_mo"]/mo_datamed["dbhmax_mo"]+
           param["onto_sq"]*mo_datamed["dbh_mo"]/mo_datamed["dbhmax_mo"]*
                            mo_datamed["dbh_mo"]/mo_datamed["dbhmax_mo"]+
           param["mo_logg"]*mo_dataLogg[i]+
           param["mo_clim"]*mo_dataClim[i]+
           param["mo_int"]*mo_dataClim[i]*mo_dataLogg[i])
  
}


## tableau d'initialisation boucle pred_cr ####
vparam<-chain_sel %>% 
  select(-lp__) %>% 
  # slice(3951:4000)
  slice((nrow(chain_sel)-50):nrow(chain_sel))


datamed<-apply(select(donnee_cr,dbh_cr,dbhmax_cr,clim_cr,logg_cr,WD_cr),2,median) # on extrait les médianes des variables
names(datamed)<-colnames(select(donnee_cr,dbh_cr,dbhmax_cr,clim_cr,logg_cr,WD_cr))

dataspe<-select(donnee_cr,dbhmax_cr,WD_cr,rgEsp_cr) %>% # tableau des variables dpdtes de l'espèce
  group_by(rgEsp_cr) %>% 
  summarise(WDspe_cr=median(WD_cr),dbhmaxspe_cr=median(dbhmax_cr))

onto<-seq(0.08,2,by=0.08) #abscisse pour calcul du modèle"
dataclim<-seq(min(donnee_cr$clim_cr),max(donnee_cr$clim_cr),
              by=(max(donnee_cr$clim_cr)-min(donnee_cr$clim_cr))/20)
datalogg<-seq(min(donnee_cr$logg_cr),max(donnee_cr$logg_cr),
              by=(max(donnee_cr$logg_cr)-min(donnee_cr$logg_cr))/20)

repet<-10 # répétitions par ligne de paramètre et par abscisse (onto)


# Boucles Fpred_cr_spe et graphe 2 variables ####
nn<-length(onto)*repet # nombre de simulations pour un jeu de paramètres 
ontos<-rep(onto,each=repet) # ontogénies pour une valeur de paramètres
# Nbesp<-nrow(tesp)
Nbesp<-2

pred_pparam<-matrix(data=NA,nrow=nn,ncol=0)
# i<-2
# j<-1
# k<-1
# for(k in 1:max(Nesp)){
temps_depart<-Sys.time()
for(j in 1:nrow(vparam)){
  srt<-rep(NA,nn)
  cr_sigGesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigGesp"]))
  cr_sigClesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigClesp"]))
  # cr_sigLoesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigLoesp"]))
  sigma<-rnorm(nn,mean=0,as.numeric(vparam[j,]["sigma"]))
  sig<-data.frame(cr_sigGesp,cr_sigClesp,
                  # cr_sigLoesp,
                  sigma)
    for(k in 1:Nbesp){
      for (i in 1:nn){
        srt[i]<-Fpred_cr_spe(i,vparam[j,],ontos,datamed,dataspe,sig,k)
      }
    pred_pparam<-cbind(pred_pparam,srt)
    colnames(pred_pparam)[k+(j-1)*Nbesp]<-paste("param",j,"_sp",k,sep="")
   }
}  
Sys.time()- temps_depart 

pred_pc<-as.data.frame(pred_pparam)
pred_pc$ontos<-ontos 

pred_pc<-pred_pc%>% 
  pivot_longer(
    cols= starts_with("param"),
    names_to="num_param",
    # names_prefix = "param",
    # names_ptypes=list(param=integer()),# pour tranformer les titre de colonne parame en entier
    values_to = "predictions") %>% 
  mutate(nesp=as.numeric(substr(num_param,nchar(num_param),nchar(num_param)))) %>% 
  left_join(tesp,by="nesp") %>% 
  arrange(nesp,ontos)

pred_pc$predictions<-unlist(pred_pc$predictions) # quantile() plante si les données sont sous forme de liste 

pred_pc<-pred_pc%>% 
  group_by(especes,ontos) %>% 
  # summarise(test=n())
  summarise(pc01=quantile(predictions,probs=0.01),
            pc05=quantile(predictions,probs=0.05),
            mediane=quantile(predictions,probs=0.50),
            pc95=quantile(predictions,probs=0.95),
            pc99=quantile(predictions,probs=0.99))

ggplot(pred_pc) +
  facet_wrap(~especes,scales="free_y")+
  geom_point(data=donnee_cr,aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1))) +
  geom_ribbon(aes(x=ontos,ymin=pc05,ymax=pc95),color="grey")+
  geom_line (aes(x=ontos,y=mediane),color="red",size=1) +
  geom_line (aes(x=ontos,y=pc01),color="grey",size=0.5,linetype = "dashed") +
  geom_line (aes(x=ontos,y=pc99),color="grey",size=0.5,linetype = "dashed")


# Boucles Fpred_cr_spe_clim et grapha 3 variables ####

nnvi<-length(onto)*repet*length(dataclim) # nombre de simulations pour un jeu de paramètres 
ontos<-rep(onto,each=(repet*length(dataclim))) # diamètres pour une valeur de paramètres
clims<-rep(rep(dataclim,each=repet),length(onto)) # indice stress hydrique pour une valeur de paramètres


# Nbesp<-nrow(tesp)
Nbesp<-6
temps_depart<-Sys.time()
pred_pparam<-matrix(data=NA,nrow=nnvi,ncol=0)
# i<-1
# j<-1
# k<-1
# for(k in 1:max(Nesp)){
for(j in 1:nrow(vparam)){
  srt<-rep(NA,nnvi)
  cr_sigGesp<-rnorm(nnvi,mean=0,as.numeric(vparam[j,]["cr_sigGesp"]))
  cr_sigClesp<-rnorm(nnvi,mean=0,as.numeric(vparam[j,]["cr_sigClesp"]))
  # cr_sigLoesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigLoesp"]))
  sigma<-rnorm(nnvi,mean=0,as.numeric(vparam[j,]["sigma"]))
  sig<-data.frame(cr_sigGesp,cr_sigClesp,
                  # cr_sigLoesp,
                  sigma)
  for(k in 1:Nbesp){
    for (i in 1:nnvi){
      srt[i]<-Fpred_cr_spe_clim(i,vparam[j,],ontos,clims,datamed,dataspe,sig,k)
    }
    
    # srt<-lapply(X=1:nn, FUN=Fpred_cr(,vparam[j,],ontos,datamed,sig))
    pred_pparam<-cbind(pred_pparam,srt)
    colnames(pred_pparam)[k+(j-1)*Nbesp]<-paste("param",j,"_sp",k,sep="")
  }
}  
Sys.time()- temps_depart 

pred_pc<-as.data.frame(pred_pparam)
pred_pc$Ontos<-ontos 
pred_pc$Clims<-clims
pred_pc$Clims_ch<-as.character(clims)

pred_pc<-pred_pc%>% 
  pivot_longer(
    cols= starts_with("param"),
    names_to="num_param",
    # names_prefix = "param",
    # names_ptypes=list(param=integer()),# pour tranformer les titre de colonne parame en entier
    values_to = "predictions") %>% 
  mutate(nesp=as.numeric(substr(num_param,nchar(num_param),nchar(num_param)))) %>% 
  left_join(tesp,by="nesp") %>% 
  arrange(nesp,Ontos,Clims)

pred_pc$predictions<-unlist(pred_pc$predictions) # quantile() plante si les données sont sous forme de liste 

pred_pc1<-pred_pc%>% 
  group_by(especes,Ontos,Clims,Clims_ch) %>% 
  # summarise(test=n())
  summarise(pc01=quantile(predictions,probs=0.01),
            pc05=quantile(predictions,probs=0.05),
            mediane=quantile(predictions,probs=0.50),
            pc95=quantile(predictions,probs=0.95),
            pc99=quantile(predictions,probs=0.99))
  


ggplot(pred_pc1) +
  facet_wrap(~especes,scales="free_y")+
  geom_point(data=donnee_cr,aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1))) +
  geom_ribbon(aes(x=Ontos,ymin=pc05,ymax=pc95),fill="grey",alpha=0.1)+
  geom_line (aes(x=Ontos,y=mediane,color=Clims_ch),size=1) 

ggplot(pred_pc1 %>% filter(especes=="Carapa surinamensis")) +
  # facet_wrap(~clims,scales="free_y")+
  geom_point(data=donnee_cr,aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1))) +
  # geom_ribbon(aes(x=Ontos,ymin=pc05,ymax=pc95),fill="grey",alpha=0.1)+
  geom_line (aes(x=Ontos,y=mediane,color=Clims_ch),size=1) 

ggplot(pred_pc1,aes(x=Ontos,y=Clims))+
  geom_raster(aes(fill=mediane))+
  geom_contour(aes(z=))
  facet_wrap(~especes)
  
  geom_line (aes(x=ontos,y=pc01),color="grey",size=0.5,linetype = "dashed") +
  geom_line (aes(x=ontos,y=pc99),color="grey",size=0.5,linetype = "dashed")
  
## tableau d'initialisation boucle pred_mo ####
vparam<-chain_sel %>% 
    select(-lp__) %>% 
    slice((nrow(chain_sel)-250):nrow(chain_sel))
  
  
mo_datamed<-apply(select(donnee_mo,dbh_mo,dbhmax_mo,clim_mo,logg_mo),2,median) # on extrait les médianes des variables
names(mo_datamed)<-colnames(select(donnee_mo,dbh_mo,dbhmax_mo,clim_mo,logg_mo))
  
mo_dmaxspe<-select(donnee_mo,dbhmax_mo,rgEsp_mo) %>% # tableau des variables dpdtes de l'espèce
    group_by(rgEsp_mo) %>% 
    summarise(dbhmaxspe_mo=median(dbhmax_mo))

diam<-seq(10,130,by=5) # exploration range des diamètres
mo_dataClim<-seq(min(donnee_mo$clim_mo),max(donnee_mo$clim_mo),
              by=(max(donnee_mo$clim_mo)-min(donnee_mo$clim_mo))/50) # exploration range des stress hydrique
mo_dataLogg<-seq(min(donnee_mo$logg_mo),max(donnee_mo$logg_mo),
              by=(max(donnee_mo$logg_mo)-min(donnee_mo$logg_mo))/50) # exploration range des indexLogg



# Boucles Fpred_mo_spe_clim et graphe a 3 variables ####

nnvi<-length(diam)*length(mo_dataClim) # nombre de simulations pour un jeu de paramètres 
diams<-rep(diam,each=(repet*length(mo_dataClim))) # diamètres pour une valeur de paramètres
clims<-rep(rep(mo_dataClim,each=repet),length(diam)) # indice stress hydrique pour une valeur de paramètres
# Nbesp<-nrow(tesp)
Nbesp<-16
temps_depart<-Sys.time()
temps_depart
pred_pparam<-matrix(data=NA,nrow=nnvi,ncol=0)
# i<-1
# j<-1
# k<-1
# for(k in 1:max(Nesp)){
pb <- txtProgressBar(min = 0, max = nrow(vparam), style = 3)
db=0

for(j in 1:nrow(vparam)){
  db=db+1
  srt<-rep(NA,nnvi)
  for(k in 1:Nbesp){
    for (i in 1:nnvi){
      srt[i]<-Fpred_mo_spe_clim(i,vparam[j,],diams,clims,mo_datamed,mo_dmaxspe,k)
    }
    
    pred_pparam<-cbind(pred_pparam,srt)
    colnames(pred_pparam)[k+(j-1)*Nbesp]<-paste("param",j,"_sp",k,sep="")
  }
  setTxtProgressBar(pb, db)
}  
Sys.time()- temps_depart 

pred_pc<-as.data.frame(pred_pparam)
pred_pc$Diams<-diams 
pred_pc$Clims<-clims
pred_pc$Clims_ch<-as.character(clims)

# boucle  Fpred_mo_inter_spe####
nnvi<-length(mo_dataLogg)*length(mo_dataClim) # nombre de simulations pour un jeu de paramètres 
loggs<-rep(mo_dataLogg,each=length(mo_dataClim)) # Logg pour une valeur de paramètres
clims<-rep(mo_dataClim,length(mo_dataLogg)) # indice stress hydrique pour une valeur de paramètres
# Nbesp<-nrow(tesp)
Nbesp<-16

temps_depart<-Sys.time()
pred_pparam<-matrix(data=NA,nrow=nnvi,ncol=0)
i<-1
j<-1
k<-1
# for(k in 1:max(Nesp)){
 for(k in 1:Nbesp){
  srt<-rep(NA,nnvi)
   for(i in 1:nnvi){
     srt[i]<-Fpred_mo_inter_spe(i,param,clims,loggs,mo_datamed,mo_dmaxspe,k)
    }
    pred_pparam<-cbind(pred_pparam,srt)
    colnames(pred_pparam)[k]<-paste("param",k,sep="")
}  
Sys.time()- temps_depart 

pred_pc<-as.data.frame(pred_pparam)
pred_pc$Loggs<-loggs 
pred_pc$Clims<-clims

save(pred_pc,file="stan_paracou30ncl_mo5_i4_predinter_spe.Rdata")

#boucle  Fpred_mo_inter####
nnvi<-length(mo_dataLogg)*length(mo_dataClim) # nombre de simulations pour un jeu de paramètres 
loggs<-rep(mo_dataLogg,each=length(mo_dataClim)) # Logg pour une valeur de paramètres
clims<-rep(mo_dataClim,length(mo_dataLogg)) # indice stress hydrique pour une valeur de paramètres
# Nbesp<-nrow(tesp)

temps_depart<-Sys.time()
pred_pparam<-matrix(data=NA,nrow=nnvi,ncol=0)
i<-1
srt<-rep(NA,nnvi)
for(i in 1:nnvi){
    srt[i]<-Fpred_mo_inter(i,param,clims,loggs,mo_datamed)
  }
pred_pparam<-cbind(pred_pparam,srt)
colnames(pred_pparam)[1]<-paste("valeurs")
Sys.time()- temps_depart 

pred_pc_nonspe<-as.data.frame(pred_pparam)
pred_pc_nonspe$Loggs<-loggs 
pred_pc_nonspe$Clims<-clims
pred_pc_nonspe$predictions<-exp(pred_pc_nonspe$valeurs)/(1+exp(pred_pc_nonspe$valeurs))

save(pred_pc_nonspe,file="stan_paracou30ncl_mo5_i4_predinter.Rdata")


## graphe prédictions ####

pred_gg<-pred_pc%>% 
  pivot_longer(
    cols= starts_with("param"),
    names_to="num_param",
    # names_prefix = "param",
    # names_ptypes=list(param=integer()),# pour tranformer les titre de colonne parame en entier
    values_to = "predictions") %>% 
  mutate(nesp=as.numeric(substr(num_param,nchar("param")+1,nchar(num_param)))) %>% 
  left_join(tesp,by="nesp")  
  # arrange(nesp,Diams,Clims)  

pred_gg$predictions<-unlist(pred_gg$predictions) # quantile() plante si les données sont sous forme de liste 
pred_gg<-pred_gg %>% 
  mutate(logitm1=exp(predictions)/(1+exp(predictions))) # on calcule le logit-1(pred) c'est à dire la probabilité de mort prévu pour chaque individus


pred_pc1<-pred_pc%>% 
  group_by(especes,Diams,Clims,Clims_ch) %>% 
  # summarise(test=n())
  summarise(pc01=quantile(predictions,probs=0.01),
            pc2_5=quantile(predictions,probs=0.025),
            pc05=quantile(predictions,probs=0.05),
            mediane=quantile(predictions,probs=0.50),
            pc95=quantile(predictions,probs=0.95),
            pc975=quantile(predictions,probs=0.975),
            pc99=quantile(predictions,probs=0.99))



ggplot(pred_pc1) +
  facet_wrap(~especes,scales="free_y")+
  geom_ribbon(aes(x=Diams,ymin=pc05,ymax=pc95),fill="grey",alpha=0.1)+
  geom_line (aes(x=Diams,y=mediane,color=Clims_ch),size=1) 

ggplot(pred_pc1 %>% filter(especes=="Carapa surinamensis")) +
  # facet_wrap(~clims,scales="free_y")+
  # geom_ribbon(aes(x=Ontos,ymin=pc05,ymax=pc95),fill="grey",alpha=0.1)+
  geom_line (aes(x=Diams,y=mediane,color=Clims_ch),size=1) 

ggplot(pred_pc1,aes(x=Diams,y=Clims))+
  geom_raster(aes(fill=mediane))+
  geom_contour(aes(z=))
facet_wrap(~especes)

geom_line (aes(x=ontos,y=pc01),color="grey",size=0.5,linetype = "dashed") +
  geom_line (aes(x=ontos,y=pc99),color="grey",size=0.5,linetype = "dashed")


ggplot(pred_gg,aes(x=Loggs,y=Clims,z=logitm1))+
  facet_wrap(~especes)+
  # geom_raster(aes(fill=predictions))+
  stat_contour(geom="polygon",aes(fill=..level..))+
  geom_tile (aes(fill=logitm1))+
  stat_contour(bins = 15) +
  xlab("Exploitation Forestière") +
  ylab("Stress Hydrique") +
  guides(fill = guide_colorbar(title = "Proba mort"))

# par espèces
ggplot(filter(pred_gg,especes=="Goupia glabra"),aes(x=Loggs,y=Clims,z=logitm1))+
  # geom_raster(aes(fill=predictions))+
  stat_contour(geom="polygon",aes(fill=..level..))+
  geom_tile (aes(fill=logitm1))+
  stat_contour(bins = 15) +
  xlab("Exploitation Forestière") +
  ylab("Stress Hydrique") +
  guides(fill = guide_colorbar(title = "Proba mort"))+
  ggtitle("Goupia glabra")


# graphe nonspe
ggplot(pred_pc_nonspe,aes(x=Loggs,y=Clims,z=predictions))+
  # geom_raster(aes(fill=predictions))+
  geom_contour(aes(colour="black"))+
  geom_tile (aes(fill=predictions))+
  stat_contour(bins = 15) +
  xlab("Exploitation Forestière") +
  ylab("Stress Hydrique") +
  guides(fill = guide_colorbar(title = "Proba mort"))

               # lineend = "butt", linejoin = "round",
               # linemitre = 10, na.rm = FALSE, show.legend = NA,
               # inherit.aes = TRUE)
               # 


## 8 tableau paramètres ####
poster_abs<-function(chaine,legende){
  ggplot(chaine) + 
    geom_freqpoly(aes(x=Abs,binwidth=0.015,color="red")) +
    # geom_freqpoly(aes(x=valeurs,binwidth=0.015),size=1,color="red") +
    geom_vline(xintercept = 0, linetype = "solid")+
    xlab(legende) +
    facet_wrap(~especes,scales="free")
  
}

tab_spe<-function(chaine){
  tab<-chaine %>% 
    group_by(especes) %>% 
    summarise(pc1_abs=quantile(Abs,probs=0.01),
              pc2_5_abs=quantile(Abs,probs=0.025),
              pc5_abs=quantile(Abs,probs=0.05),
              mediane_abs=median(Abs),
              pc95_abs=quantile(Abs,probs=0.95),
              pc975_abs=quantile(Abs,probs=0.975),
              pc99_abs=quantile(Abs,probs=0.99))
  return(tab)
  
}

tab_nonspe<-function(chaine,parametre){
  tab<-chaine %>% 
    filter(variables%in%parametre) %>% 
    group_by(variables) %>% 
    summarise(pc1=quantile(valeurs,probs=0.01),
              pc2_5=quantile(valeurs,probs=0.025),
              pc5=quantile(valeurs,probs=0.05),
              mediane=median(valeurs),
              pc95=quantile(valeurs,probs=0.95),
              pc975=quantile(valeurs,probs=0.975),
              pc99=quantile(valeurs,probs=0.99))
  return(tab)
}
  
if(mtype%in%c("croiss","joint")){
## Climat
  chain_clim<-chain_sel %>% 
     select(one_of(pars_crCl),-cr_sigClesp) %>% #on laisse le vecteur cr_clim qui sera dupliquer lors du pivot_longer
     pivot_longer(
        cols= starts_with("cr_Clesp"),
        names_to="variables",
        values_to = "valeurs") %>% 
     mutate(Abs=cr_clim+valeurs) %>% 
     mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                            as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                            as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
     left_join(tesp,by="nesp")  

  poster_abs(chain_clim,"Paramètres Clim croissance")
  cr_clim_pc_abs<-tab_spe(chain_clim)
  
  cr_pc_nonspe<-tab_nonspe(chain_ggpl,pars_cr1)

  save(cr_clim_pc_abs,file=paste("PC_",code_sim,"_cr_climAbs.Rdata",sep=""))
  write.csv2(cr_clim_pc_abs,file=paste(code_sim,"_cr_clim_pc_abs.csv",sep=""),dec=".")
  save(cr_pc_nonspe,file=paste("PC_",code_sim,"_nonspe.Rdata",sep=""))
  write.csv2(cr_pc_nonspe,file=paste(code_sim,"_nonspe.csv",sep=""),dec=".")
  
# logg
  chain_logg<-chain_sel %>% 
   select(one_of(pars_crLogg),-cr_sigLoesp) %>% 
   pivot_longer(
    cols= starts_with("cr_Loesp"),
    names_to="variables",
    values_to = "valeurs") %>% 
   mutate(Abs=cr_logg+valeurs) %>% 
   mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                     as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
   left_join(tesp,by="nesp")  

  poster_abs(chain_logg,"Paramètres Logg croissance")
  cr_logg_pc_abs<-tab_spe(chain_logg)
  
  save(cr_logg_pc_abs,file=paste("PC_",code_sim,"_cr_LoggAbs.Rdata",sep=""))
  write.csv2(cr_logg_pc_abs,file=paste(code_sim,"_cr_logg_pc_abs.csv",sep=""),dec=".")

  # interaction clima_Logg
  chain_int<-chain_sel %>% 
   select(one_of(pars_crInt),-cr_sigInesp) %>% 
   pivot_longer(
    cols= starts_with("cr_Inesp"),
    names_to="variables",
    values_to = "valeurs") %>% 
   mutate(Abs=cr_int+valeurs) %>% 
   mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                     as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
   left_join(tesp,by="nesp")  

  poster_abs(chain_int,"Paramètres int croissance")
  cr_int_pc_abs<-tab_spe(chain_int)
  
 ggplot(chain_int) + 
  geom_freqpoly(aes(x=Abs,binwidth=0.015,color="red")) +
  # geom_freqpoly(aes(x=valeurs,binwidth=0.015),size=1,color="red") +
  geom_vline(xintercept = 0, linetype = "solid")+
  xlab("Paramètres Logg croissance") +
  facet_wrap(~especes,scales="free")

  save(cr_int_pc_abs,file=paste("PC_",code_sim,"_cr_IntAbs.Rdata",sep=""))
  write.csv2(cr_int_pc_abs,file=paste(code_sim,"_cr_int_pc_abs.csv",sep=""),dec=".")
}
if(mtype%in%c("morta","joint")){
## Climat
 chain_clim<-chain_sel %>% 
   select(one_of(pars_moCl),-mo_sigClesp) %>% 
   pivot_longer(
    cols= starts_with("mo_Clesp"),
    names_to="variables",
    values_to = "valeurs") %>% 
   mutate(Abs=mo_clim+valeurs) %>% 
   mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                     as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
   left_join(tesp,by="nesp")  

  poster_abs(chain_clim,"Paramètres Clim mortalité")
  mo_clim_pc_abs<-tab_spe(chain_clim)

  mo_pc_nonspe<-tab_nonspe(chain_ggpl,pars_mo1)

  save(mo_clim_pc_abs,file=paste("PC_",code_sim,"_mo_ClimAbs.Rdata",sep=""))
  write.csv2(mo_clim_pc_abs,file=paste(code_sim,"_mo_clim_pc_abs.csv",sep=""),dec=".")
  save(mo_pc_nonspe,file=paste("PC_",code_sim,"_nonspe.Rdata",sep=""))
  write.csv2(mo_pc_nonspe,file=paste(code_sim,"_nonspe.csv",sep=""),dec=".")
  
# logg
  chain_logg<-chain_sel %>% 
   select(one_of(pars_moLogg),-mo_sigLoesp) %>% 
   pivot_longer(
    cols= starts_with("mo_Loesp"),
    names_to="variables",
    values_to = "valeurs") %>% 
   mutate(Abs=mo_logg+valeurs) %>% 
   mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                     as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
   left_join(tesp,by="nesp")  

  poster_abs(chain_logg,"Paramètres Logg mortalité")
  mo_logg_pc_abs<-tab_spe(chain_logg)

  save(mo_logg_pc_abs,file=paste("PC_",code_sim,"_mo_LoggAbs.Rdata",sep=""))
  write.csv2(mo_logg_pc_abs,file=paste(code_sim,"_mo_logg_pc_abs.csv",sep=""),dec=".")

# interaction climat Logg

  chain_int<-chain_sel %>% 
   select(one_of(pars_moInt),-mo_sigInesp) %>% 
   pivot_longer(
    cols= starts_with("mo_Inesp"),
    names_to="variables",
    values_to = "valeurs") %>% 
   mutate(Abs=mo_int+valeurs) %>% 
   mutate(nesp=ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",
                     as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)))) %>% 
   left_join(tesp,by="nesp")  

  poster_abs(chain_int,"Paramètres Interaction mortalité")
  mo_int_pc_abs<-tab_spe(chain_int)

  save(mo_int_pc_abs,file=paste("PC_",code_sim,"_mo_IntAbs.Rdata",sep=""))
  write.csv2(mo_int_pc_abs,file=paste(code_sim,"_mo_int_pc_abs.csv",sep=""),dec=".")
}

## 9 graphe vulnérabilité ####
# load("PC_paracou30ncl_cr5_i4_pc_IntAbs.Rdata") # cr_int_pc_abs


graphe_vuln<-function(tab_mo,tab_cr_,titre){
  Vuln4graph<-as.data.frame(matrix(data=NA,nrow=16,ncol=0))
  Vuln4graph$VulnGmed<-tab_cr_$mediane_abs
  Vuln4graph$VulnGmedlow<-tab_cr_$pc2_5_abs
  Vuln4graph$VulnGmedhigh<-tab_cr_$pc975_abs
  Vuln4graph$VulnMmed<-tab_mo$mediane_abs
  Vuln4graph$VulnMmedlow<-tab_mo$pc2_5_abs
  Vuln4graph$VulnMmedhigh<-tab_mo$pc975_abs
  Vuln4graph$Species<-tab_mo$especes
 
  ggplot(data=Vuln4graph, aes(x=VulnMmed, y=VulnGmed, col=Species, label=Species,alpha=1),alpha=1) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin=VulnGmedlow, ymax=VulnGmedhigh, col=Species),
                    size=0.5,alpha=0.5, linetype = "solid") +
    geom_errorbarh(aes(xmin= VulnMmedlow, xmax=VulnMmedhigh, col=Species),
                   height = 0, size=0.5,alpha=0.5, linetype = "solid") +
    xlab("Impact sur la mortalité") + ylab("Impact sur la croissance") + 
    theme_bw() + theme(text = element_text(size=15))+
    ggtitle(titre)+ theme(plot.title = element_text(size=15,hjust=0.55))
}

#climat
load(paste("PC_",code_sim,"_mo_climAbs.Rdata",sep="")) # mo_clim_pc_abs
load(paste("PC_",code_sim,"_cr_ClimAbs.Rdata",sep="")) # cr_clim_pc_abs
graphe_vuln(mo_clim_pc_abs,cr_clim_pc_abs,"Effet du Stress Hydrique")


#logg
load(paste("PC_",code_sim,"_mo_LoggAbs.Rdata",sep="")) # mo_logg_pc_abs
load(paste("PC_",code_sim,"_cr_LoggAbs.Rdata",sep="")) # cr_logg_pc_abs
graphe_vuln(mo_logg_pc_abs,cr_logg_pc_abs,"Effet de l'exploitation")

#int
load(paste("PC_",code_sim,"_mo_IntAbs.Rdata",sep="")) # mo_logg_pc_abs
load(paste("PC_",code_sim,"_cr_ClimAbs.Rdata",sep="")) # cr_logg_pc_abs
graphe_vuln(mo_int_pc_abs,cr_int_pc_abs,"Effet de l'interaction")

# plot
ggplot(data=Vuln4graph, aes(x=VulnMmed, y=VulnGmed, col=Species, label=Species,alpha=1),alpha=1) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin=VulnGmedlow, ymax=VulnGmedhigh, col=Species),
                  size=0.5,alpha=0.5, linetype = "solid") +
  geom_errorbarh(aes(xmin= VulnMmedlow, xmax=VulnMmedhigh, col=Species),
                 height = 0, size=0.5,alpha=0.5, linetype = "solid") +
  xlab("Impact sur la mortalité") + ylab("Impact sur la croissance") + 
  theme_bw() + theme(text = element_text(size=15))+
  ggtitle("Effet du stress hydrique")+ theme(plot.title = element_text(size=15,hjust=0.55))
    
# geom_label(data=dta_label,aes(x=Cx,y=Cy,label=Clabel))

# plot
# ggplot(data=Vuln4graph, aes(x=VulnMmed, y=VulnGmed, col=Species, label=SpCode)) +
# ggplot(data=Vuln4graph, aes(x=VulnMmed, y=VulnGmed, col=Species, label=Species)) +
#   geom_point(size=3) +
#   geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 1, linetype = "dashed") +
#   xlab("Impact on mortality") + ylab("Impact on growth") +
#   # geom_text_repel(hjust=-0.1, vjust=-0.5, size=8) +
#   theme_bw() + theme(legend.position = "none") + theme(text = element_text(size=25))


## 10 graphe indices REW journaliers ####
clim_site<-REW_sites_gfclim %>%  # REW_sites_gfclim est issu du ficher Donnees_guyafor_pre-traitement.Rmd
  mutate(REW04=ifelse(REW>=0.4,NA,REW)) %>% 
  filter(Forest=="Organabo") %>% # choix du site
  filter(year>=2010&year<2019) # choix des années

ggplot(clim_site)+
  geom_line(aes(x=day_julian, y=REW), color="blue")+
  facet_wrap(~year,scale="free",nrow=2)+ # fixe le nombre de ligne de graphes
  geom_hline(yintercept = 0.4, linetype = "dashed")+
  geom_ribbon(aes(x=day_julian,ymin=REW04, ymax=0.4,fill="grey"),color="blue")+
  xlab("site Organabo (Mana)") + # toitre du graphe sur l'axe des X
  theme(legend.position = "none") # pour ne pas faire apparaitre la légende



# tests ####


datacr_gg<-as.data.frame(cbind(datacr,AGR1))
ggplot(datacr_gg, aes(diam_1, log(AGR1+1))) +
  geom_point() +
  geom_line(aes(y = predict(datacr[,"diam_1"],var1[,1])), col = "red")


-	Pour déclarer où est C++14, il te faut un fichier Makevars.win, et le remplir :
  https://community.rstudio.com/t/error-in-shlib-internal-args-c-14-standard-requested-but-cxx14-is-not-defined/16819/2,
cf la contribution de Mara
-	Pour savoir où est ton Makevars, utiliser les commandes tools::makevars_user() et tools::makevars_site() 
https://stackoverflow.com/questions/43597632/understanding-the-contents-of-the-makevars-file-in-r-macros-variables-r-ma. 
Sur mon ordinateur, il n’y en a pas mais je n’utilise pas RStan. C’est dans un de ces fichiers que l’emplacement de C++ est enregistré, dans la « variable » CXX14.

pb <- txtProgressBar(min = 0, max = total, style = 3)
i=0
for (k in DB_dup$idTree) {
  i=i+1
  
  setTxtProgressBar(pb, i)
}


