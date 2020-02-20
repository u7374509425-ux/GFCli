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


#### initialisation et chargement des données #####

Library <- function(Packages) {
  InstallAndLoad <- function(Package) {
    if (!Package %in% installed.packages()[, 1]) {install.packages(Package)}
    require(Package, character.only = TRUE)
  }
  invisible(sapply(Packages, InstallAndLoad))
}

# Ajouter les packages necessaires ici
Library(c("knitr","raster","tidyverse","dplyr",
          "ggplot2","rstan","bayesplot","rstanarm"
          ))
# attention à l'ordre d'installation des packages, conflit de fonctions entre raster et tidyverse

# library(knitr)
# library(raster)
# # library(leaflet)
# library(tidyverse)
# library(rstan)
# library(bayesplot)
# # library(shinystan)
# library(rstanarm)
# library(brms)
# library(dplyr)
# library(ggplot2)

# load("prepa_donnees/Data16esp_s0_to_join.Rdata") # issu de Donnees_guyafor-pre-traitement.Rmd
load("prepa_donnees/Data_guyafor_16esp_Cl_join.Rdata")
# metadata 
     #nbmin_arbres,            # nombre minimal d'arbres par esp et par site
     #nbmin_mes,               # nombre minimal de mesures par arbres
     #an_min_meteo,            # donnee meteo les plus anciennes
     #tabMesuresSelect,        # table de données
     #tab_arbres_select,       # Id et statu des arbres sélectionnés
     #tab_esp,                 # effectif de d'arbre et de mesure par esp et par traitement 
     #tab_esp_site_logg,       # effectif d'arbres et mesure par esp, par site et par traitement
     #"tab_site_logg,           # effectif d'arbres et mesures par site, par traitement et par espèce
     #"esp_range_dmax,          # Min et max des dmax95 par espèce
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
#fich_sauv<-"Data_16esp_allyear.Rdata"
#defaut
code_esp_cible<-"multiesp"
espliste<-espGFClim
parametres<-c(FALSE, NA,FALSE,FALSE,NA,FALSE,paste(espliste,collapse = ", "),code_esp_cible)
names(parametres)<-c("selann","annee_u","delrecru","selsite","site","selEsp","espliste","codeEsp")
annee_u <-sort(unique(c(datacr_cl_lo$Year1,datacr_cl_lo$Year2,datamo_cl_lo$Year1,datacr_cl_lo$Year3)))

#parametres
selann<-FALSE            # si TRUE selection des annees via un vecteur ci dessous
annee_sel <-seq(1991,2012)
# delrecru<-TRUE         # si TRUE suppression des recruts = arbres absent du premier inventaire, necessite la tablede mesure source "tabMesuresSelect"

selsite<-TRUE
# site<-c("Paracou")
site<-c("Acarouany","BAFOG","Laussat","Montagne Plomb","Montagne Tortue","Nouragues",
"Régina St Georges","Tibourou","Trésor") # tous sauf Paracou  soit 7% des données

selEsp<-FALSE
#esp_sel<-c("Sextonia rubra","Manilkara bidentata","Dicorynia guianensis","Goupia glabra")
esp_sel<-c("Bocoa prouacensis", "Carapa surinamensis", "Chrysophyllum sanguinolentum", "Dicorynia guianensis",
           "Goupia glabra", "Qualea rosea", "Sextonia rubra", "Symphonia globulifera", "Virola michelii")


code_sorties<-"s_paracou"

datacr_s<-datacr_cl_lo
datamo_s<-datamo_cl_lo

"avant selection"
nrow(datacr_s)
nrow(datamo_s)

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
if (selann) {
  parametres[1]<-selann
  parametres[2]<-paste(as.character(min(annee_sel)),"-",as.character(max(annee_sel)))
  annees_u<-annee_sel
  datacr_s<-filter(datacr_s,(Year1 %in% annee_u) & (Year2 %in% annee_u))
  datamo_s<-filter(datamo_s,(Year1 %in% annee_u) & (Year2 %in% annee_u) & (Year3 %in% annee_u))
}

"après selann"
nrow(datacr_s)
nrow(datamo_s)


# 
# if (delrecru) {
#   parametres[3]<-delrecru
#   
#  unrecruts<-tabMesuresSelect %>% 
#   filter(CensusYear==min(annee_sel)) %>% 
#   filter(!is.na(CircCorr)) %>% 
#   select(idTree,Forest,Plot,SubPlot,codeEsp)
#  unrecruts<-unique(unrecruts)
#  datacr_s<-semi_join(datacr_s,unrecruts,by="idTree")
#  datamo_s<-semi_join(datamo_s,unrecruts,by="idTree")
# }    
# 
# "après delrecru"
# nrow(datacr_s)
# nrow(datamo_s)

 
if (selsite) {
  parametres[4]<-selsite
  parametres[5]<-paste(site,collapse = ", ") # Colle tous les element de site sur une seule chaine de caracteres
  datacr_s<-filter(datacr_s,Forest %in% site)
  datamo_s<-filter(datamo_s,Forest %in% site)
  code_esp_cible<-code_sorties
}  

if (selEsp) {
  espliste<-esp_sel
  parametres[6]<-selEsp
  parametres[7]<-paste(espliste,collapse = ", ")
  code_esp_cible<-code_sorties
  datacr_s<-filter(datacr_s,idEsp %in% espliste)
  datamo_s<-filter(datamo_s,idEsp %in% espliste)
}

"après selsite et selEsp"
nrow(datacr_s)
nrow(datamo_s)

##3 Centrage et réduction des variables Logg et climat
CentreReduit <- function(vec)
  return(vec-mean(vec,na.rm=T))/sd(vec,na.rm=T)

datacr_s$IshInv<-CentreReduit(datacr_s$IshInv)
datamo_s$IshInvMo<-CentreReduit(datamo_s$IshInvMo)

datamo_s$TxLogg<-CentreReduit(datamo_s$TxLogg)
datacr_s$TxLogg<-CentreReduit(datacr_s$TxLogg)


##4 graphes et tableau résumé 
#tableau arbre en ligne et annee de mesure en colonne
ArbresAnnee<-datacr_s%>%
  select(idTree,idEsp,idLogg,Diam1,Forest,Year1,Dmax95)%>%
  spread(key=Year1,value=Diam1)

# tableau résumé espèces
tabtemp_cr<-datacr_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) %>% 
  arrange(desc(nb_acc))

tabtemp<-datamo_s%>%
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
  left_join(tabtemp_cr,by="idEsp")

ggplot(tabtemp, aes(x=reorder(idEsp,nb_acc),y=nb_acc))+
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


# distribution des variables
ggplot(datacr_s, aes(x=reorder(Forest,nb_acc),y=nb_acc))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.7))

 hist(datacr_s$IshInv)
 hist(datacr_cl_lo$TxLogg)
# hist(datacr_s$AGR)

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
  
  # apurement éventuel
#visudatacr <- datacr %>%
 #   filter(is.na(IshInv))  
           
#   datacr <- datacr%>%
 #  filter() # arbre des pb d'accroissement (gain de 25 cm sur le diamètre, )
 
#temp <- datamo%>%
#filter(is.na(IshInvMo))
  

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
    WD_cr = datacr$WD,        # Wood desnsity pour modèle de croissance
    morts = datamo$to_Death,     # vecteur de 1 = evenement de mort et de 0 = evenements de survie
    rgtree_cr = datacr$TreeRank ,        # rang des individu dans la liste  arbres pour croissance
    rgtree_mo = datamo$TreeRank,       # rang des individu dans la liste arbres pour la mortalite
    rgEsp_cr = datacr$EspRank,        # rang des espece dans la liste des espece pour croissance
    rgEsp_mo = datamo$EspRank        # rang des espece dans la liste des espece pour la mortalité
  )

  save(dataj,parametres,tesp,code_esp_cible,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data.Rdata'))
  save(dataj,parametres,tesp,code_esp_cible,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data2.Rdata'),version=2)  # pour export vers Rstudio serveur
  
#### C- Lancement des chaines et étude de convergence ####
  
##1 graphe de vérif ####
  
  graphAcc<-function (dataf,n){
    Dgraph<-dataf[["dbh_cr"]][dataf[["rgEsp_cr"]]==n]/dataf[["dbhmax_cr"]][dataf[["rgEsp_cr"]]==n]
    Accgraph<-dataf[["Acc_cr"]][dataf[["rgEsp_cr"]]==n]
    plot(x=Dgraph,y=Accgraph)
  }
  
  # plot(x=dataj[["dbh_cr"]],y=dataj[["Acc_cr"]]) # observation accroissement : il faut borner éventuellement le paramètre Gmax et Dopt en fonction de ce nuage

  # graphe pour n especes
par(mfrow=c(1,dataj[["Nesp"]]))
for(i in 1:dataj[["Nesp"]])  graphAcc(dataj,i)

 # observation accroissement : il faut borner éventuellement le paramètre Gmax et Dopt en fonction de ce nuage

par(mfrow=c(1,1))
graphAcc(dataj,1) 

Dplot<-dataj[["dbh_cr"]][dataj[["rgEsp_cr"]]==4]/dataj[["dbhmax_cr"]][dataj[["rgEsp_cr"]]==4]
Accplot<-dataj[["Acc_cr"]][dataj[["rgEsp_cr"]]==4]
points(x=Dplot,y=Accplot,col="green") 

##2 lancement des chaines ####
# modèle complet ####
temps_depart <-Sys.time()
fitj_c <- stan('join_model4_multiesp.stan', data = dataj,chain=4)
Sys.time()- temps_depart
save(fitj_c,parametres,file=paste('stan_sorties/stan_',code_esp_cible,'_join_vcr_sortie.Rdata')) # voi debut du fichier stan pour les specificite

# selection e parametres de sortiespour exclure des paramètre lies au effet aleatoire : un par observation ("logacc_mu_cr","logit_mo", "logacc_mu_mo", "vig_mo")
# non testé
#pars_save<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","oo_logit","vig","onto","onto_sq","mo_clim","mo_logg","sigma","cr_sigGesp","cr_sigClesp",
#             "cr_sigLoesp","cr_Clesp","cr_Clesp","cr_Loesp","mo_sigesp","mo_sigClesp",
#             "mo_sigLoesp","mo_esp","mo_Clesp","mo_Loesp")
#temps_depart <-Sys.time()
#fitj_c <- stan('join_model4_multiesp.stan', data = dataj,chain=4,pars=pars_save,include=TRUE)
#Sys.time()- temps_depart
#save(fitj_c,file=paste('stan_sorties/stan_',code_esp_cible,'_join_multiesp_sortie.Rdata')) # voi debut du fichier stan pour les specificite



# modèle de croissance seul ####
pars_save<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","sigma","cr_sigGesp","cr_sigClesp",
             "cr_sigLoesp","cr_Clesp","cr_Clesp","cr_Loesp")

temps_depart <-Sys.time()
fitj_cr <- stan('cr_model4_multiesp_camila.stan', data = dataj,pars=pars_save,include=TRUE,
                chain=4,
                iter=4000,warmup=3000,
                # control = list(adapt_delta = 0.99,max_treedepth = 15)
                )
Sys.time()- temps_depart
save(fitj_cr,parametres,file=paste('stan_sorties/stan_',code_esp_cible,'_cr_camila_vcr_3alea_sortie_i4.Rdata'))
#3alea = 3 effet aleatoire espèce sur le modèle de croissance
#  temps_depart <-Sys.time()
#  fitjo <- stan('join_model4_multiesp.stan', data = dataj)
#  Sys.time()- temps_depart]
#  save(fitjo,file=paste('stan_',code_esp_cible,'_joint_sortie.Rdata'))

##1 vecteurs de paramètres  et calcul de rhat  ####

chain<-fitj_cr
# pars<-chain@model_pars

# paste(chain@model_pars,collapse = ",") # pour sélectionner facilement les paramètre dans la console

# liste des noms des colonnes de parametres pour graphes de sortie. A limiter aux parametres  : exclure transformed parameters et generated quantities 
espliste<-parametres["espliste"]

cr_Gesp_multi<-paste("cr_Gesp[",1:length(espliste),"]",sep="")
cr_Clesp_multi<-paste("cr_Clesp[",1:length(espliste),"]",sep="")
cr_Loesp_multi<-paste("cr_Loesp[",1:length(espliste),"]",sep="")

pars1<-c("oo_Gmax","Ks","Dopt","cr_dmax","sigma","cr_sigGesp")
parsG<-c("oo_Gmax","cr_sigGesp",cr_Gesp_multi)
parsCl<-c("cr_clim","cr_sigClesp",
          cr_Clesp_multi)
parsLogg<-c("cr_logg","cr_sigLoesp",
            cr_Loesp_multi)


# modèle joint sans interactions Logg climat
# pars<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","oo_logit","vig","onto","onto_sq","mo_clim","mo_logg","sigma","cr_sigGesp","cr_sigClesp",
# "cr_sigLoesp","cr_Clesp","cr_Clesp","cr_Loesp","mo_sigesp","mo_sigClesp",
# "mo_sigLoesp","mo_esp","mo_Clesp","mo_Loesp")

#modèle de croissance sans interactions logg climat
pars_cr<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","sigma","cr_sigGesp","cr_sigClesp",
        "cr_sigLoesp","cr_Clesp","cr_Clesp","cr_Loesp","lp__")
#parametres modèle de croissance sans parametres par esp
pars_non_spe<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","sigma","cr_sigGesp","cr_sigClesp",
        "cr_sigLoesp","lp__")

# modèle de mortalité sans interactions logg climat
# pars <- c("Gmax", "Dopt", "Ks","cr_clim","vig","onto","onto_sq","mo_clim","sigma", "sigGt") 

pars<-pars_cr


# calcul rhat
print(chain, pars = pars_cr) 
print(chain, pars = cr_Clesp_multi) 

##2 Traces des chaines ####
traceplot(chain, pars=c("lp__")) # plot les chaines dela vraisemblance

traceplot(chain, pars=pars, nrow=4) # nrow : nombre de lignes de graphe
# mcmc_trace(as.array(chain), pars = c(pars),facet_args = list(labeller = label_parsed))
# mcmc_trace(as.array(chain), pars = c(pars),facet_args = list(labeller = label_parsed),window =c(1,700))



## tracer des chaines choisie : ex suppression d'une chaine non convergente et choix de paramètres
chain_df<-as.data.frame(chain) 
chain_df$chaines<-as.factor(c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000)))
chain_df$iterations<-rep(1:1000,4)
chain_sel<-chain_df %>% 
  filter(chaines!=4) %>% 
  filter (chaines !=1) 
  # select(-starts_with("logacc_mu_cr")) %>% 
  # select(-starts_with("logacc_mu_mo")) %>% 
  # select(-starts_with("logit_mo")) %>% 
  # select(-starts_with("vig_mo")) 

espmodel<-as.data.frame(strsplit(parametres["espliste"],", "))
espmodel$espliste<-as.character(espmodel$espliste)
espmodel$nesp<-1:length(espmodel$espliste)

chain_ggpl<-chain_sel %>% 
  pivot_longer(-c(iterations,chaines),names_to ="variables",values_to = "valeurs" ) %>% 
  mutate(test=ifelse(substr(variables,nchar(variables),nchar(variables))=="]",
                     ifelse(substr(variables,nchar(variables)-2,nchar(variables)-2)=="[",1,2)
                     ,0)) %>% 
  mutate(nesp=ifelse(test==1,as.numeric(substr(variables,nchar(variables)-1,nchar(variables)-1)),
                     ifelse(test==2,as.numeric(substr(variables,nchar(variables)-2,nchar(variables)-1)),
                            -1))) %>% 
  left_join(espmodel,by="nesp") %>% 
  mutate(especes=ifelse(is.na(espliste),variables,espliste)) %>% 
  select(-test,-espliste)


# chain_ggpl<-chain_sel %>% 
#   pivot_longer(-c(iterations,chaines),names_to ="variables",values_to = "valeurs" ) %>% 
#   filter(variables %in% pars)

ggplot(chain_ggpl) + 
  geom_line(aes(x = iterations, y = valeurs ,color= chaines)) +
  xlab("itérations") +
  facet_wrap(~variables,scales="free_y")


##3 nuage de points pour paramètres donnés ####

mcmc_pairs(as.array(chain), pars = pars1)
mcmc_pairs(as.array(chain), pars = parsG)
mcmc_pairs(as.array(chain), pars = parsCl)
mcmc_pairs(as.array(chain), pars = parcLogg)
# GGally::ggpairs() pour faire des pairs plus lisibles



##4 Posteriors ####

#for(i in 1:length(pars))
#  mcmc_areas(as.array(chain), prob = 0.8,pars = pars[i])

#mcmc_areas(as.array(chain), prob = 0.8,pars = pars)
# mcmc_areas(as.array(chain), prob = 0.8,pars = pars1)
mcmc_areas(as.array(chain), prob = 0.8,pars = parsG)
mcmc_areas(as.array(chain), prob = 0.8,pars = parsCl)
mcmc_areas(as.array(chain), prob = 0.8,pars = parsLogg[1:3])
mcmc_areas(as.array(chain), prob = 0.8,pars = parsLogg[4:6])

mcmc_areas(as.array(chain), prob = 0.8,pars = c("oo_Gmax","Ks","Dopt","cr_dmax","sigma","cr_sigGesp"))
mcmc_areas(as.array(chain), prob = 0.8,pars = c("oo_Gmax","cr_sigGesp","Ks"))
mcmc_areas(as.array(chain), prob = 0.8,pars = c("Dopt"))
mcmc_areas(as.array(chain), prob = 0.8,pars = c("cr_dmax"))
mcmc_areas(as.array(chain), prob = 0.8,pars = c("sigma"))

ggplot(chain_ggpl %>% filter(variables%in%parsG)) + 
  geom_freqpoly(aes(x=valeurs,binwidth=0.015,color= chaines)) +
  xlab("Parametres OO de G") +
  geom_vline(xintercept = 0, linetype = "solid")+
  facet_wrap(~especes)

ggplot(chain_ggpl %>% filter(variables%in%parsCl)) + 
  geom_freqpoly(aes(x=valeurs,binwidth=0.015,color= chaines)) +
  xlab("Parametres Climat") +
  geom_vline(xintercept = 0, linetype = "solid")+
  facet_wrap(~especes)

ggplot(filter(chain_ggpl,variables%in%parsLogg)) + 
  geom_freqpoly(aes(x=valeurs,binwidth=0.015,color= chaines)) +
  xlab("Parametres Logg") +
  geom_vline(xintercept = 0, linetype = "solid")+
  facet_wrap(~especes)
#launch_shinystan(chain)

#### D- Predictions et graphe de sortie ####
temp<-chain_sel %>%   # on reprends le tableau avec seulement les "bonnes " chaines
  select(-chaines,-iterations)

###1 calcul des prédiction pour chaque observé
param<-apply(temp,2,median) # on extrait les médianes des paramètres
names(param)<-colnames(temp)


donnee_cr<-as.data.frame(dataj[c(5,7,10,12,16,18,20,22)]) #ensemble des donnees liees au modele de croissance
donnee_mo<-as.data.frame(dataj[-c(1:5,7,10,12,16,18,20,22)]) #donnees du modele de mortalité

nesp<-donnee_cr$rgEsp_cr

# on utilise les valeurs des effet aléatoire de chaque espèce
donnee_cr$pred <- (param["oo_Gmax"]+ param[paste("cr_Gesp[",nesp,"]",sep="")] +
              param["cr_dmax"]*donnee_cr$dbhmax+
              (param["cr_clim"]+param[paste("cr_Clesp[",nesp,"]",sep="")])* donnee_cr$clim_cr+
              (param["cr_logg"]+param[paste("cr_Loesp[",nesp,"]",sep="")])* donnee_cr$logg_cr)*
  exp((-0.5)*(log(donnee_cr$dbh_cr/(donnee_cr$dbhmax*param["Dopt"]))/(param["Ks"]*donnee_cr$WD_cr))*
             (log(donnee_cr$dbh_cr/(donnee_cr$dbhmax*param["Dopt"]))/(param["Ks"]*donnee_cr$WD_cr)))

donnee_cr<-donnee_cr %>% left_join(transmute(tesp,idEsp=idEsp,rgEsp_cr=rank),by="rgEsp_cr")

# temps_depart<-Sys.time()
ggplot(donnee_cr) +
  geom_point(aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1))) +
  geom_line (aes(x=dbh_cr/dbhmax_cr,y=pred),color="red",size=1) +
  facet_wrap(~idEsp)
# Sys.time()- temps_depart 


# pour accélérer la boucle, stocker 1000 tirages de rnorm dans un vecteur 
# et ensuite aller chercher les valeurs dans ce vecteur

# i<-2
# vparamj<-vparam[j,]
# ontos<-Abscisses
# sigmas<-sig

###2 calcul de prévision et erreur pour diamètre (ontogénie) donné par espèce et sur un range de variables climatiques

Fpred_cr<-function(i,vparamj,ontos,vmed,sigmas){
  return (as.numeric(vparamj["oo_Gmax"])+as.numeric(sigmas[i,]["cr_sigGesp"])+
            as.numeric(vparamj["cr_dmax"])*vmed["dbhmax_cr"]+
                (as.numeric(vparamj["cr_clim"])+as.numeric(sigmas[i,]["cr_sigClesp"]))*vmed["clim_cr"]+
                (as.numeric(vparamj["cr_logg"])+as.numeric(sigmas[i,]["cr_sigLoesp"]))*vmed["logg_cr"])*
    exp((-0.5)*(log(onto[i]/as.numeric(vparamj["Dopt"]))/(as.numeric(vparamj["Ks"])*vmed["WD_cr"]))*
               (log(onto[i]/as.numeric(vparamj["Dopt"]))/(as.numeric(vparamj["Ks"])*vmed["WD_cr"])))+
    sigmas[i,]["sigma"]
}


vparam<-chain_sel %>% 
  select(oo_Gmax,Ks,Dopt,cr_clim,cr_logg,cr_dmax,sigma,cr_sigGesp,cr_sigClesp,cr_sigLoesp) #%>% 
  # slice(50:100)

vmed<-apply(select(donnee_cr,dbhmax_cr,clim_cr,logg_cr,WD_cr),2,median) # on extrait les médianes des variables
names(vmed)<-colnames(select(donnee_cr,dbhmax_cr,clim_cr,logg_cr,WD_cr))

onto<-seq(0.08,2,by=0.08) #abscisse pour calcul du modèle"

repet<-100 # répétitions par ligne de paramètre et par abscisse (onto)
nn<-length(onto)*repet # nombre de simulations pour un jeu de paramètres 
abscisses<-rep(onto,each=repet) # abscisses pour une valeur de paramètres
pred_pparam<-matrix(data=NA,nrow=nn,ncol=0)

# j<-2

temps_depart<-Sys.time()
for(j in 1:nrow(vparam)){
  srt<-rep(NA,nn)
  cr_sigGesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigGesp"]))
  cr_sigClesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigClesp"]))
  cr_sigLoesp<-rnorm(nn,mean=0,as.numeric(vparam[j,]["cr_sigLoesp"]))
  sigma<-rnorm(nn,mean=0,as.numeric(vparam[j,]["sigma"]))
  sig<-data.frame(cr_sigGesp,cr_sigClesp,cr_sigLoesp,sigma)
  
  for (i in 1:nn){
    srt[i]<-Fpred_cr(i,vparam[j,],ontos,vmed,sig)
    }
    
  # srt<-lapply(X=1:nn, FUN=Fpred_cr(,vparam[j,],abscisses,vmed,sig))
  pred_pparam<-cbind(pred_pparam,srt)
  colnames(pred_pparam)[j]<-paste("param",j,sep="")
}
Sys.time()- temps_depart 


# colnames(pred_pparam)<-rep("")
#calcul des medianes et pc par abscisses
# pred_ontos$Abscisse<-abscisses
pred_pc<-as.data.frame(pred_pparam)
pred_pc$abscisses<-abscisses 

pred_pc<-pred_pc%>% 
  pivot_longer(
    cols= starts_with("param"),
    names_to="num_param",
    # names_prefix = "param",
    # names_ptypes=list(param=integer()),# pour tranformer les titre de colonne parame en entier
    values_to = "predictions") %>% 
  group_by(abscisses) %>% 
  summarise(pc01=quantile(predictions,probs=0.01),
            pc05=quantile(predictions,probs=0.05),
            mediane=quantile(predictions,probs=0.50),
            pc95=quantile(predictions,probs=0.95),
            pc99=quantile(predictions,probs=0.99))

  
  

ggplot(pred_pc) +
  geom_point(data=donnee_cr,aes(x=dbh_cr/dbhmax_cr,y=log(Acc_cr+1)))+ 
  geom_ribbon (aes(x=abscisses,ymin=pc05,ymax=pc95),color="grey")+
  geom_line (aes(x=abscisses,y=mediane),color="red",size=1) +
  geom_line (aes(x=abscisses,y=pc01),color="grey",size=0.5,linetype = "dashed") +
  geom_line (aes(x=abscisses,y=pc99),color="grey",size=0.5,linetype = "dashed") 



par(mfrow=c(1,1))
Dplot<-dataj[["dbh_cr"]]
Accplot<-dataj[["Acc_cr"]]
plot(x=Dplot,y=Accplot) # observation accroissement : il faut borner éventuellement le paramètre Gmax et Dopt en fonction de ce nuage
points(x=Dplot,y=pred_cr,col="blue")


par(mfrow=c(1,2))
Dplot<-dataj[["dbh_cr"]][dataj[["rgEsp_cr"]]==1]/dataj[["dbhmax_cr"]][dataj[["rgEsp_cr"]]==1]
Accplot<-dataj[["Acc_cr"]][dataj[["rgEsp_cr"]]==1]
plot(x=Dplot,y=Accplot) # observation accroissement : il faut borner éventuellement le paramètre Gmax et Dopt en fonction de ce nuage
points(x=Dplot,y=pred_cr[dataj[["rgEsp_cr"]]==1],col="blue")

Dplot<-dataj[["dbh_cr"]][dataj[["rgEsp_cr"]]==2]/dataj[["dbhmax_cr"]][dataj[["rgEsp_cr"]]==2]
Accplot<-dataj[["Acc_cr"]][dataj[["rgEsp_cr"]]==2]
points(x=Dplot,y=Accplot,col="red") # observation accroissement : il faut borner éventuellement le paramètre Gmax et Dopt en fonction de ce nuage
points(x=Dplot,y=pred_cr[dataj[["rgEsp_cr"]]==2],col="blue")







#chain_pred<- chain_mat %>%  # si effet aléatoire autant de colonne que d'observations * 4000 lignes
#  select(starts_with("predacc1"))

pars_median<-apply(as.matrix(chain_pars),2,median)

diam_ini<-dataj[["dbh1"]]
obs1<-dataj[["Acc1"]]

#logacc_mu_cr[n1] = (oo_Gmax+cr_Clesp[rgEsp_cr[n1]]+cr_clim*clim_cr[n1]+cr_logg*logg_cr[n1])*exp ((-0.5)*pow(log(dbh_cr[n1]/(dbhmax_cr[n1]*(Dopt+cr_Desp[rgEsp_cr[n1]])))/Ks,2)); // obligation de faire un boucle sinon erreur de la fonction pow()

mod<-pars_median["Gmax"]*exp(-0.5*(log(diam_ini/pars_median["Dopt"])/pars_median["Ks"])^2)
mod<-exp(mod)-1



plot(x = diam_ini,y=obs1)
#points(x = diam_ini,y=pred1,col="red")
points(x = diam_ini,y=mod,col="blue")



chain_mat<-as.matrix(chain) # lignes : les 1000 derniere valeur des 4 chaines, colonnes : les parametres + lp__ ( vraisemblance)+ les Nt effet aleatoire de Gt, puis les Nt *nb generated quantities soit pres de 12000 colonnes 
chain_pars <- chain_mat[, pars] # selection des colonnes parametres
pars_mat<-dimnames(chain_mat)[[2]]
crmu1<-pars_mat[which("cr_mu1" %in% pars_mat[])]
nb_crmu1<-grep("cr_mu1",pars_mat)

predict <- function(dbh,var1) 
  (pars_opt["Gmax"]+pars_opt["cr_clim"]*var1)*exp(-0.5*(log(dbh/pars_opt["Dopt"])/pars_opt["Ks"])^2)

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


