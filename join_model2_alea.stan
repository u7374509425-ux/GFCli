data {                        // Data block
  int<lower=0> N1;             // Sample size pop1 : accroissement observes
  int<lower=0> Nmo;             // Sample size pop evenements de survie ou de mort
  vector<lower=0> [N1] Acc1;          // accroissement entre n-1 et n
  vector<lower=0> [Nmo] Acc_mo;          // accroissement entre n-2 et n-2 suivi d'un evenement de survie ou de mort
  vector [N1] clim1;            // covariable climat n-1
  vector [Nmo] clim_mo;         // covariable climat n-2 calcul pred accroissement avant un evenement de survie
  vector [Nmo] climP;         // covariable climat n-1 calcul logit d'un evenement de survie
  vector<lower=0> [N1] dbh1;          // covariable diametre initial modele croissance
  vector<lower=0> [Nmo] dbh_mo;          // covariable diametre n-2 survie calcul proba mort avant un evenement de survie
  int<lower=0,upper=1> morts[Nmo];      // vecteur de 0 = evenement de survie et de 1 evenement de mort
  int numt1 [N1];      // numero d'arbre dans la liste des arbres pour modele accroissement
  int numt_mo[Nmo];      // numero d'arbre dans la liste des arbres pour modele mortalite
  int <lower=0> Nt;      // longueur liste des arbres
}
parameters {
  real <lower=0.01> Gmax;  // Parameter croissance
  real <lower=0> Ks;       // Parameter croissance
  real <lower=10> Dopt;    // Parameter croissance
  real  cr_clim;           // Parameter croissance
  real  vig;               // Parameter mortalite
  real <upper=0> onto;      // Parameter mortalite
  real <lower=0> onto_sq;   // Parameter mortalite
  real mo_clim;             // Parameter mortalite
  real <lower=0> sigma;     // variance
  real <lower=0> sigGt;     // variance effet alea sur Gmax
  vector [Nt] Gt;            //effet aleatoire arbre sur Gmax : a declare ici et non dans la partie modele car il entre 
}
transformed parameters {    // Transformed parameters block
// modele croissance n-1 n  
  vector [N1] lcr_mu1;            
  vector [Nmo] logit_mo;       // nombre reel=logit(proba mort)    
  vector [Nmo] lcr_mo_mu;          // prediction accroisssement  
  vector [Nmo] vig_mo;           // vigueur = log(acc_obs/acc_pred)

// modele d'acroissement
  for(n1 in 1:N1) {
    lcr_mu1[n1] = (Gmax+Gt[numt1[n1]]+cr_clim*clim1[n1])*exp ((-0.5)*pow(log(dbh1[n1]/Dopt)/Ks,2)); // obligation de faire un boucle sinon erreur de le fonction pow()
  }

// modele mortalite : data_survie
  for(n2 in 1:Nmo) {
    lcr_mo_mu[n2] = (Gmax+Gt[numt_mo[n2]]+cr_clim*clim_mo[n2])*exp ((-0.5)*pow(log(dbh_mo[n2]/Dopt)/Ks,2)); 
    vig_mo[n2]=log(Acc_mo[n2]+1)-lcr_mo_mu[n2];                      
    logit_mo[n2]= vig*vig_mo[n2]+onto*dbh_mo[n2]+onto_sq*pow(dbh_mo[n2],2)+mo_clim*climP[n2]; // logit(mort reussie) 
  }
} 

model {                            // Model block
  Gmax ~ normal(0,100);       // priors
  Ks ~ normal(0,100);         // priors
  Dopt ~ normal(0,100);       // priors
  cr_clim ~ normal(0,100);    // priors
  vig ~ normal(0,100);         // priors
  onto ~ normal(0,100);        // priors
  onto_sq ~ normal(0,100);    // priors
  mo_clim ~ normal(0,100);;   // priors
  sigma ~ gamma(10^2,10^2);       // priors
  sigGt ~ gamma(10^2,10^2);       // priors

  for(nt in 1:Nt){                // priors
    Gt[nt]~normal(0,sigGt);
  }

 for(n1 in 1:N1)
   log(Acc1[n1]+1)~normal(lcr_mu1[n1],sigma);       // likelihood
// log(Acc1 +1)~normal(cr_mu1,sigma) ; /// prends 2.5 fois plus de temps qu'avec la boucle !!
  
  for(n2 in 1:Nmo)
    morts[n2]~bernoulli_logit(logit_mo[n2]); // proba de realisation de l'evenement de mort ou de survie
    //  surv2~bernoulli_logit(logit_mo2);
  

}
generated quantities {        // Generated quantities block. 
  vector [N1] predacc1 ;      // prediction acroissement
  vector [N1] resacc1 ;      // prediction acroissement
  vector [Nmo] predmort;      // prediction survie

  for(n1 in 1:N1) {
   predacc1[n1]=exp(lcr_mu1[n1])-1;       
   resacc1[n1]=Acc1[n1] - predacc1[n1];
  }
  resacc1 = Acc1 - predacc1;
   
  for(n2 in 1:Nmo){
    if (morts[n2]==0) predmort[n2]=1/(1+exp((logit_mo[n2])));  // proba survie 
      else predmort[n2]=1/(1+exp(-(logit_mo[n2]))); // proba mort 
  }
}
