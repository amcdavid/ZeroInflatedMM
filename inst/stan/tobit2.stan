/*
Models:
1) Marginal model
2) Reduced rank VGLM
3) factor analytic GLM


*/

data {
  int<lower=0> N; //number of obs
  int<lower=0> Tf; //number of fixed covariates
  int<lower=0> G; //number of genes
  int<lower=0> NposG; //number of positive obs (length of indexing vector)

  matrix[N, Tf] x; //fixed covariate matrix, fixed across genes
  
  

  //matrix[N, G] y;  //zero-inflated data
  int<lower=0, upper=1> v[G,N]; //indicators
  /*
   For each gene, we will need to index into the linear predictor
   to pull out the positive obs.
   Ipos contains the indices (concatenated over genes) of the positive obs
   IposGI contains the starting index for each gene in Ipos
   We only store the y that are positive (and rely again on the indexing from Ipos)
  */	
  int<lower=1, upper=N> Ipos[NposG]; 
  int<lower=1, upper=NposG+1> IposGI[G+1]; 
  vector[NposG] y;
  
}

/*
transformed data{
  vector<lower=-1, upper=1>[G] rho_G;
  for(g in 1:G){
  	rho_G[g] = 0;
  }

}
*/

parameters {
  vector[Tf] betaD_Tf_G[G]; //fixed effect regression coefficients
  vector[Tf] betaM_Tf_G[G]; //marginal effect
  
  vector<lower=0>[G] sigmaC_G; //dispersion
  real<lower=0> gamaC; //hyper dispersion
  real<lower=0> gambC;

  vector[Tf] muC_Tf; //fixed effect anchors
  vector[Tf] muD_Tf; 
  vector<lower=0>[Tf] sigmaD_Tf; //fixed effect variance
  vector<lower=0>[Tf] sigmaC_Tf;

  vector<lower=-1, upper=1>[G] rho_G;
  //vector<lower=0,upper=1>[G] rho_GUnit; //covariance
  //real<lower=0> betaA; //covariance hyper params
  //real<lower=0> betaB;
}


/* transformed parameters {
  vector<lower=-1,upper=1>[G] rho_G; //covariance
  rho_G = rho_GUnit*2-1;
} */



model {
 //print("X dim", dims(x));	
 //print("betaD dim", dims(betaD));
 
  //rho_GUnit ~ beta(betaA, betaB);
  muC_Tf ~ normal(0, 4);
  muD_Tf ~ normal(0, 2);
  sigmaD_Tf ~ student_t(6, 1, 2);
  sigmaC_Tf ~ student_t(6, 1, 2);
  gamaC ~ student_t(6, 1, 2);
  gambC ~ student_t(6, 1, 2);
  // sigmaD_Tf ~ units are between gene variation of effect, g prior might make sense
  // sigmaD_Tr ~ same 
  for (g in 1:G){
      vector[N] etaD;
      int iposStart;
      int iposEnd;
      vector[IposGI[g+1]- IposGI[g]] etaC;
      //int ipos[N];
      int lenIpos;

      // Because local variable must precede all statements, these are over-allocated
      iposStart = IposGI[g];
      iposEnd = IposGI[g+1]-1;
      lenIpos = iposEnd-iposStart;
      //ipos = Ipos[iposStart:iposEnd];
      //print("G=", g, "; iposStart=", iposStart, "; iposEnd=", iposEnd, "; ipos=", Ipos[iposStart:iposEnd]);

      
      //should turn this into some sort of multi normal
      betaD_Tf_G[g] ~ normal(muD_Tf, sigmaD_Tf);
      betaM_Tf_G[g] ~ normal(muC_Tf, sigmaC_Tf);
      rho_G[g] ~ normal(0, .3);

      etaD = -x*betaD_Tf_G[g];
      etaC = (x[Ipos[iposStart:iposEnd],:]*betaM_Tf_G[g]);
      sigmaC_G[g] ~ gamma(gamaC, gambC);
      //log normal likelihood
      y[iposStart:iposEnd] ~ normal(etaC, sigmaC_G[g]);
      etaD[Ipos[iposStart:iposEnd]] = ( -etaD[Ipos[iposStart:iposEnd]]  + (rho_G[g]/sigmaC_G[g])*(y[iposStart:iposEnd]-etaC))/sqrt(1-rho_G[g]^2);

      for (n in 1:N){
      	  target += log(Phi_approx(etaD[n]));
      }
      //print("len(y)= ", num_elements(y[iposStart:iposEnd]), "; len(etaC)=", num_elements(etaC[Ipos[iposStart:iposEnd]]), "; len(sigma)=", num_elements(sigmaC_G[g]));

      }

}


