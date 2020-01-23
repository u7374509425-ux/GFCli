## Mise en forma donnee pour modele Stan
## Modele joint de croissance et de mortalite

###Résumé
## Etape A : simplification optionnelle des donnée
#      - suppression des recrut
#      - choix espèce ou sites 
## Etape B : calculs et ajouts des colonnes pour des effets aléatoires
## Etape C : lancement d'une chaine avec un dataj
## Etape D : mise en forme des sorties



#### initialisation et chargement des données #####
rm(list=ls(all=TRUE))
gc() # garbage collector


library(knitr)
library(raster)
#library(leaflet)
library(tidyverse)
library(rstan)
library(bayesplot)
#library(shinystan)
library(rstanarm)
library(brms)
library(dplyr)

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
parametres<-c(FALSE, NA,FALSE,FALSE,NA,FALSE,NA,code_esp_cible)
names(parametres)<-c("selann","annee_u","delrecru","selsite","site","selEsp","espliste","codeEsp")
annee_u <-sort(unique(c(datacr_cl_lo$Year1,datacr_cl_lo$Year2,datamo_cl_lo$Year1,datacr_cl_lo$Year3)))
espliste<-espGFClim

#parametres
selann<-FALSE            # si TRUE selection des annees via un vecteur ci dessous
annee_sel <-seq(1991,2012)
delrecru<-TRUE         # si TRUE suppression des recruts = arbres absent du premier inventaire, necessite la tablede mesure source "tabMesuresSelect"

selsite<-FALSE
site<-c("Paracou")

selEsp<-FALSE
#esp_sel<-c("Sextonia rubra","Manilkara bidentata","Dicorynia guianensis","Goupia glabra")
esp_sel<-c("Sextonia rubra","Manilkara bidentata","Ruizterania albiflora","Symphonia globulifera")
code_sorties<-"Sr_Mb_Ra_Sg_CLmulti"

datacr_s<-datacr_cl_lo
datamo_s<-datamo_cl_lo

"avant selection"
nrow(datacr_s)
nrow(datamo_s)

# tableau résumé espèces
tabtemp_cr<-datacr_s%>%
  filter(Forest%in%site)%>% 
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) 

tabtemp<-datamo_s%>%
  filter(Forest%in%site)%>% 
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



if (delrecru) {
  parametres[3]<-delrecru
  
 unrecruts<-tabMesuresSelect %>% 
  filter(CensusYear==min(annee_sel)) %>% 
  filter(!is.na(CircCorr)) %>% 
  select(idTree,Forest,Plot,SubPlot,codeEsp)
 unrecruts<-unique(unrecruts)
 datacr_s<-semi_join(datacr_s,unrecruts,by="idTree")
 datamo_s<-semi_join(datamo_s,unrecruts,by="idTree")
}    

"après delrecru"
nrow(datacr_s)
nrow(datamo_s)

 
if (selsite) {
  parametres[4]<-selsite
  parametres[5]<-paste(site,collapse = ", ") # Colle tous les element de site sur une seule chaine de caracteres
  datacr_s<-filter(datacr_s,Forest %in% site)
  datamo_s<-filter(datamo_s,Forest %in% site)
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
datacr_s$IshInv<-(datacr_s$IshInv-mean(datacr_s$IshInv,na.rm=T))/sd(datacr_s$IshInv,na.rm=T)
datamo_s$IshInvMo<-(datamo_s$IshInvMo-mean(datamo_s$IshInvMo,na.rm=T))/sd(datamo_s$IshInvMo,na.rm=T)
datamo_s$IshInvVig<-(datamo_s$IshInvVig-mean(datamo_s$IshInvVig,na.rm=T))/sd(datamo_s$IshInvVig,na.rm=T)


##4 Pour visu construction d'un tableau arbre en ligne et annee de mesure en colonne
ArbresAnnee<-datacr_s%>%
  select(idTree,idEsp,idLogg,Diam1,Forest,Year1,Dmax95)%>%
  spread(key=Year1,value=Diam1)
# tableau résumé espèces
tabtemp_cr<-datacr_s%>%
  filter(Forest%in%site)%>% 
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_cr=n(),nb_acc=sum(nbmes)) 

tabtemp<-datamo_s%>%
  filter(Forest%in%site)%>% 
  #  semi_join(unrecruts,by="idTree") %>%
  group_by(idEsp,idTree) %>% 
  summarise(nbmes=n()) %>% 
  group_by(idEsp) %>% 
  summarise(nbtree_mo=n(),nb_mo=sum(nbmes)) %>% 
  left_join(tabtemp_cr,by="idEsp")



#### B- Calculs et ajouts des colonnes pour des effet aléatoire ####

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

  save(dataj,parametres,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data.Rdata'))
  save(dataj,parametres,file=paste('stan_sorties/stan_',code_esp_cible,'_multisite_vcr_data2.Rdata'),version=2)  # pour export vers Rstudio serveur
  
#### C- Lancement des chaines ####
  
##1 graphe de vérif ##
  
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

  #2 lancement des chaines
# modèle complet
temps_depart <-Sys.time()
fitj_c <- stan('join_model4_multiesp.stan', data = dataj,chain=4)
Sys.time()- temps_depart
save(fitj_c,file=paste('stan_sorties/stan_',code_esp_cible,'_join_multiesp_sortie.Rdata')) # voi debut du fichier stan pour les specificite

# selection e parametres de sortiespour exclure des paramètre lies au effet aleatoire : un par observation ("logacc_mu_cr","logit_mo", "logacc_mu_mo", "vig_mo")
# non testé
#pars_save<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","oo_logit","vig","onto","onto_sq","mo_clim","mo_logg","sigma","cr_sigGesp","cr_sigClesp",
#             "cr_sigLoesp","cr_Gesp","cr_Clesp","cr_Loesp","mo_sigesp","mo_sigClesp",
#             "mo_sigLoesp","mo_esp","mo_Clesp","mo_Loesp")
#temps_depart <-Sys.time()
#fitj_c <- stan('join_model4_multiesp.stan', data = dataj,chain=4,pars=pars_save,include=TRUE)
#Sys.time()- temps_depart
#save(fitj_c,file=paste('stan_sorties/stan_',code_esp_cible,'_join_multiesp_sortie.Rdata')) # voi debut du fichier stan pour les specificite

#modèle de croissance seul
pars_save<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","sigma","cr_sigGesp","cr_sigClesp",
             "cr_sigLoesp","cr_Gesp","cr_Clesp","cr_Loesp")

temps_depart <-Sys.time()
fitj_cr <- stan('cr_model4_multiesp_camila.stan', data = dataj,pars=pars_save,include=TRUE,
                chain=4,
                iter=2000,warmup=1000,
                # control = list(adapt_delta = 0.99,max_treedepth = 15)
)
Sys.time()- temps_depart
save(fitj_cr,file=paste('stan_sorties/stan_',code_esp_cible,'_cr_camila_multiesp_3alea_sortier.Rdata'))
#3alea = 3 effet aleatoire espèce sur le modèle de croissance
#  temps_depart <-Sys.time()
#  fitjo <- stan('join_model4_multiesp.stan', data = dataj)
#  Sys.time()- temps_depart]
#  save(fitjo,file=paste('stan_',code_esp_cible,'_joint_sortie.Rdata'))


#### D- Conctruction des graphes de sortie ####
## 
  
chain<-fitj_cr
# pars<-chain@model_pars

# paste(chain@model_pars,collapse = ",") # pour sélectionner facilement les paramètre dans la console

# liste des noms des colonnes de parametres pour graphes de sortie. A limiter aux parametres  : exclure transformed parameters et generated quantities 

# modèle joint sans interactions Logg climat
# pars<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","oo_logit","vig","onto","onto_sq","mo_clim","mo_logg","sigma","cr_sigGesp","cr_sigClesp",
# "cr_sigLoesp","cr_Gesp","cr_Clesp","cr_Loesp","mo_sigesp","mo_sigClesp",
# "mo_sigLoesp","mo_esp","mo_Clesp","mo_Loesp")

#modèle de croissance sans interactions logg climat
pars<-c("oo_Gmax","Ks","Dopt","cr_clim","cr_logg","cr_dmax","sigma","cr_sigGesp","cr_sigClesp",
        "cr_sigLoesp","cr_Gesp","cr_Clesp","cr_Loesp","lp__")

# modèle de mortalité sans interactions logg climat
# pars <- c("Gmax", "Dopt", "Ks","cr_clim","vig","onto","onto_sq","mo_clim","sigma", "sigGt") 

## calcul de rhat par paramètres ####
print(chain, pars = c(pars)) 

## Traces des chaines ####
traceplot(chain, pars=c("lp__")) # plot les chaines dela vraisemblance

traceplot(chain, pars=pars, nrow=4) # nrow : nombre de lignes de graphe
# mcmc_trace(as.array(chain), pars = c(pars),facet_args = list(labeller = label_parsed))
# mcmc_trace(as.array(chain), pars = c(pars),facet_args = list(labeller = label_parsed),window =c(1,700))



## tracer des chaines choisie : ex suppression d'une chaine non convergente et choix de paramètres
chain_df<-as.data.frame(chain) 
chain_df$chaines<-as.factor(c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000)))
chain_df$iterations<-rep(1:1000,4)
chain_sel<-chain_df %>% 
  filter(chaines!=2)
  # select(-starts_with("logacc_mu_cr")) %>% 
  # select(-starts_with("logacc_mu_mo")) %>% 
  # select(-starts_with("logit_mo")) %>% 
  # select(-starts_with("vig_mo")) 

chain_ggpl<-chain_sel %>% 
  pivot_longer(-c(iterations,chaines),names_to ="variables",values_to = "valeurs" )

ggplot(chain_ggpl) + 
  geom_line(aes(x = iterations, y = valeurs ,color= chaines)) +
  xlab("itérations") +
  facet_wrap(~variables,scales="free_y")


## nuage de points pour paramètre donnés
pars1<-c("oo_Gmax","Ks","Dopt","cr_dmax","sigma","cr_sigGesp")
parsG<-c("oo_Gmax","cr_Gesp[1]", "cr_Gesp[2]", "cr_Gesp[3]", "cr_Gesp[4]")
parsCl<-c("cr_clim","cr_sigClesp",
         "cr_Clesp[1]", "cr_Clesp[2]", "cr_Clesp[3]", "cr_Clesp[4]")
parsLogg<-c("cr_logg","cr_sigLoesp",
            "cr_Loesp[1]","cr_Loesp[2]","cr_Loesp[3]","cr_Loesp[4]")


pars_pairs<-unique(chain_ggpl %>% select(variables))
pars_pairs<-pars_pairs$variables
mcmc_pairs(as.array(chain), pars = pars1)
mcmc_pairs(as.array(chain), pars = parsG)
mcmc_pairs(as.array(chain), pars = parsCl)
mcmc_pairs(as.array(chain), pars = parcLogg)
# GGally::ggpairs() pour faire des pairs plus lisibles


ggplot(chain_ggpl) +
  geom_point()

## Posteriors
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

#launch_shinystan(chain)

#### E- Predictions ####
nesp<-1
prediction <- function(donnee,param,nesp) {
cr_Gesp<-param[paste("cr_Gesp[",nesp,"]",sep="")]
cr_Clesp<-param[paste("cr_Clesp[",nesp,"]",sep="")]
cr_Loesp<-param[paste("cr_Loesp[",nesp,"]",sep="")]

logacc <- (param["oo_Gmax"]+ cr_Gesp +
              param["cr_dmax"]*donnee["dbhmax"]+
              (param["cr_clim"]+cr_Clesp)* donnee["clim_cr"]+
              (param["cr_logg"]+cr_Loesp)* donnee["logg_cr"])*
  exp((-0.5)*pow(log(donnee["dbh_cr"]/(donnee["dbhmax"]*param["Dopt"]))/(param["Ks"]*donnee["WD_cr"]),2)) 
}
temp<-chain_sel %>% 
  select(-chaines,-iterations)
  
param<-apply(temp,2,median) 
names(param)<-colnames(temp)
chain_sel %>% select(-chaines,-iterations))
  
#pars <- c("Gmax", "Dopt", "Ks","cr_clim","vig","onto","onto_sq","mo_clim","sigma","sigGt") 
#chain_pars <- chain_mat %>%
#  select(one_of(pars)) # selection des colonnes paramètres


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

#logacc_mu_cr[n1] = (oo_Gmax+cr_Gesp[rgEsp_cr[n1]]+cr_clim*clim_cr[n1]+cr_logg*logg_cr[n1])*exp ((-0.5)*pow(log(dbh_cr[n1]/(dbhmax_cr[n1]*(Dopt+cr_Desp[rgEsp_cr[n1]])))/Ks,2)); // obligation de faire un boucle sinon erreur de la fonction pow()

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


