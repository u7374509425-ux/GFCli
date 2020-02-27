// modele join croissance et mortalite
// reduit au modele de mortalité
// pas d'effet vigueur
// effet logg effet climat 
// effet alea especes sur logit, cima et logg
// Interaction Logg*clima avec effet aléatoire
//  logit_mo[n2]= oo_logit+
//                onto*dbh_mo[n2]/dbhmax_mo[n2]+
//                onto_sq*pow(dbh_mo[n2]/dbhmax_mo[n2],2)+
//                mo_clim*clim_mo[n2]+
//                mo_logg*logg_mo[n2]+
//                mo_int*clim_mo[n2]*logg_mo[n2]; // logit(mort reussie) 


data {                           // Data block
//  int<lower=0> Ncr;               // Sample size pop1 : accroissement observes
  int<lower=0> Nmo;              // Sample size pop evenements de survie ou de mort
//  int <lower=0> Nt;              // longueur liste des arbres
  int <lower=0> Nesp;            // longueur liste des especes
//  vector<lower=0> [Ncr] Acc_cr;     // accroissement entre n-1 et n
  vector<lower=0> [Nmo] Acc_mo;  // accroissement entre n-2 et n-2 suivi d'un evenement de survie ou de mort
//  vector [Ncr] clim_cr;             // covariable climat n-1 modele de croissance
  vector [Nmo] clim_vig;          // covariable climat n-2 calcul pred accroissement avant un evenement de survie
  vector [Nmo] clim_mo;            // covariable climat n-1 calcul logit d'un evenement de survie
//  vector [Ncr] logg_cr;             // covariable exploitation modele de croissance
  vector [Nmo] logg_mo;          // covariable exploitation pred accroissement avant un evenement de survie
//  vector<lower=0> [Ncr] dbh_cr;     // covariable diametre initial modele croissance
  vector<lower=0> [Nmo] dbh_mo;   // covariable diametre n-1 calcul proba mort avant un evenement de mort/survie
  vector<lower=0> [Nmo] dbh_vig; // covariable diametre n-2 calcul pred accroissement avant un evenement de survie
  vector<lower=0> [Nmo] dbhmax_mo;   // diametre max de l'espèces par foret (une donnée par observation)
//  vector<lower=0> [Ncr] dbhmax_cr;   // diametre max de l'espèces par foret (une donnée par observation)
  vector<lower=0> [Nmo] WD_mo;   // densité de bois de l'espece, pour mortalite
//  vector<lower=0> [Ncr] WD_cr;   // densité de bois de l'espece, pour croissance
  int<lower=0,upper=1> morts[Nmo];// vecteur de 0 = evenement de survie et de 1 = evenement de mort
//  int rgEsp_cr [Ncr];              // rang d'espece dans la liste des especes pour modele accroissement
  int rgEsp_mo [Nmo];           // rang d'espece dans la liste des especes pour modele de mortalite
}
//  int rgtree_cr [N1];              // numero d'arbre dans la liste des arbres pour modele accroissement
//  int rgtree_mo[Nmo];            // numero d'arbre dans la liste des arbres pour modele mortalite

 

parameters {
//  real <lower=0.01> oo_Gmax;  // Parameter croissance : ordonnee a l'origine de Gmax
//  real <lower=0> Ks;       // Parameter croissance
//  real <lower=0> Dopt;    // Parameter croissance
//  real  cr_clim;           // Parameter croissance
//  real  cr_logg;           // Parameter croissance
//  real  cr_dmax;           // Parameter croissance
  real  oo_logit;               // Parameter mortalite : ordonnee a l'origine du logit
//  real  vig;               // Parameter mortalite
//  real <upper=0> onto;      // Parameter mortalite
  real  onto;      // Parameter mortalite
  real <lower=0> onto_sq;   // Parameter mortalite
  real  mo_clim;             // Parameter mortalite
  real  mo_logg;           // Parameter mortalite
  real  mo_int;           // Parameter mortalite
//  real <lower=0> sigma;     // variance modèle de croissance
  real <lower=0> mo_sigOoesp;     // variance effet alea espece sur Gmax
  real <lower=0> mo_sigClesp;     // variance effet alea espece sur Gmax
  real <lower=0> mo_sigLoesp;     // variance effet alea espece sur Gmax
  real <lower=0> mo_sigInesp;     // variance effet alea espece sur Gmax
  vector [Nesp] mo_Ooesp;     //effet aleatoire especes sur oo du logit : a declare ici et non dans la partie modele car il entre dans le calcul de transformed parameters
  vector [Nesp] mo_Clesp;     //effet aleatoire especes sur Climat du logit : a declare ici et non dans la partie modele car il entre dans le calcul de transformed parameters
  vector [Nesp] mo_Loesp;     //effet aleatoire especes sur Logf du logit : a declare ici et non dans la partie modele car il entre dans le calcul de transformed parameters
  vector [Nesp] mo_Inesp;     //effet aleatoire especes sur Logf du logit : a declare ici et non dans la partie modele car il entre dans le calcul de transformed parameters
}


transformed parameters {    // Transformed parameters block
// modele croissance n-1 n  
//  vector [Ncr] logacc_mu_cr;            
  vector [Nmo] logit_mo;       // nombre reel=logit(proba mort)    
//  vector [Nmo] logacc_mu_mo;          // prediction log(accroisssement) (donnees mortalite) = logacc_pred
//  vector [Nmo] vig_mo;           // vigueur = log(acc_obs/acc_pred)

// modele d'acroissement
//  for(n1 in 1:Ncr){ 
//    logacc_mu_cr[n1] = (oo_Gmax+cr_Gesp[rgEsp_cr[n1]]+
//                         cr_dmax*dbhmax_cr[n1]+
//                         cr_clim*clim_cr[n1]+
//                         cr_logg*logg_cr[n1])*
//                      exp((-0.5)*
//                      pow(log(dbh_cr[n1]/(dbhmax_cr[n1]*Dopt))/(Ks*WD_cr[n1]),2)); // obligation de faire un boucle sinon erreur de la fonction pow()
//  } // avec effet alea sur oo
//  for(n1 in 1:Ncr){ 
//    logacc_mu_cr[n1] = (oo_Gmax+cr_dmax*dbhmax_cr[n1]+cr_clim*clim_cr[n1]+cr_logg*logg_cr[n1])*exp((-0.5)*pow(log(dbh_cr[n1]/(dbhmax_cr[n1]*Dopt))/(Ks*WD_cr[n1]),2)); // obligation de faire un boucle sinon erreur de la fonction pow()
//  } // sans effet alea

// modele mortalite : data_survie
  for(n2 in 1:Nmo) {
//    logacc_mu_mo[n2] = (oo_Gmax+cr_Gesp[rgEsp_mo[n2]]+cr_clim*clim_vig[n2]+cr_logg*logg_mo[n2])*exp ((-0.5)*pow(log(dbh_vig[n2]/dbhmax_mo[n2]/Dopt)/(Ks*WD_mo[n2]),2)); 
//    vig_mo[n2]=log(Acc_mo[n2]+1)-logacc_mu_mo[n2];                      
//    logit_mo[n2]= oo_logit+mo_vig*vig_mo[n2]+onto*dbh_mo[n2]/dbhmax_mo[n2]+onto_sq*pow(dbh_mo[n2]/dbhmax_mo[n2],2)+mo_clim*clim_mo[n2]+mo_logg*logg_mo[n2]+mo_esp[rgEsp_mo[n2]]; // logit(mort reussie) 
    logit_mo[n2]= oo_logit+
                  mo_Ooesp[rgEsp_mo[n2]]+
                  onto*dbh_mo[n2]/dbhmax_mo[n2]+
                  onto_sq*pow(dbh_mo[n2]/dbhmax_mo[n2],2)+
                  (mo_clim+mo_Clesp[rgEsp_mo[n2]])*clim_mo[n2]+
                  (mo_logg+mo_Loesp[rgEsp_mo[n2]])*logg_mo[n2]+
                  (mo_int+mo_Inesp[rgEsp_mo[n2]])*clim_mo[n2]*logg_mo[n2]; // logit(mort reussie) 

  }
}



model {                            // Model block
//  oo_Gmax ~ normal(0,100);       // priors
//  Ks ~ normal(0,100);         // priors
//  Dopt ~ normal(0,100);       // priors
//  cr_clim ~ normal(0,100);    // priors
//  cr_logg ~ normal(0,100);    // priors
//  sigma ~ gamma(10^2,10^2);       // priors
  mo_sigOoesp ~ normal(0,5);;       // priors toncature min sur 0
  mo_sigClesp ~ normal(0,5);;       // priors toncature min sur 0
  mo_sigLoesp ~ normal(0,5);;       // priors toncature min sur 0
  mo_sigInesp ~ normal(0,5);;       // priors toncature min sur 0
  oo_logit ~ normal(0,100);    // priors
//  vig ~ normal(0,100);         // priors
  onto ~ normal(0,100);        // priors
  onto_sq ~ normal(0,100);    // priors
  mo_clim ~ normal(0,100);   // priors
  mo_logg ~ normal(0,100);   // priors
  mo_int ~ normal(0,100);   // priors

  for(nesp in 1:Nesp){                // priors
    mo_Ooesp[nesp]~normal(0,mo_sigOoesp);
    mo_Clesp[nesp]~normal(0,mo_sigClesp);
    mo_Loesp[nesp]~normal(0,mo_sigLoesp);
    mo_Inesp[nesp]~normal(0,mo_sigInesp);
  }
//  for(n1 in 1:Ncr){
//    log(Acc_cr[n1]+1)~normal(logacc_mu_cr[n1],sigma);       // likelihood
//  }
  for(n2 in 1:Nmo)
    morts[n2]~bernoulli_logit(logit_mo[n2]); // proba de realisation de l'evenement de mort ou de survie

}
