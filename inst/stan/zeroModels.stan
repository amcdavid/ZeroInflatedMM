/*
Models:
1) Marginal model
2) Reduced rank VGLM
3) factor analytic GLM
*/

data {
  // FIXED effects
  int<lower=0> N; //number of obs
  int<lower=0> Tf; //number of fixed covariates
  int<lower=0> G; //number of genes
  matrix[N, Tf] x; //fixed covariate matrix, fixed across genes
  
  // random effect shizzle
  int<lower=0> Tr; //number of random covariates
  int<lower=0> Nr; //number of random groups
  int<lower=0, upper=N> rr[Nr+1]; //starting index of random group (hence must be sorted..)
  matrix[N, Tr] xr; //random effects covariate matrix, fixed across genes


  // zero-inflated shizzle
  int<lower=0, upper=1> marginal; //use the marginal model, or not.
  int<lower=0> NposG; //number of positive obs (length of indexing vector)
  /*
   For each gene, we will need to index into the linear predictor
   to pull out the positive obs.
   Ipos contains the indices (concatenated over genes) of the positive obs
   IposGI contains the starting index for each gene in Ipos
   We only store the y that are positive (and rely again on the indexing from Ipos)
 */
 
  int<lower=1, upper=N> Ipos[NposG]; 
  int<lower=1, upper=NposG+1> IposGI[G+1];

  // DATA
  int<lower=0, upper=1> v[G,N]; //indicators
  vector[NposG] y;
  
}

transformed data {
  // Index into concatenated parameter vectors
  int id_Tr = 1:Tr;
  int ic_Tr = (Tr+1):(2*Tr);
  int id_Tf = 1:Tf;
  int ic_Tf = (Tf+1):(2*Tf);
}

parameters {
  vector[2*Tf] beta_Tf_G[G]; //fixed effect regression coefficients
  
  vector<lower=0>[G] sigmaC_G; //dispersion
  /* real<lower=0> gamaC; //hyper dispersion
  real<lower=0> gambC; */

  vector[2*Tf] mu_Tf; //fixed effect anchors
  vector<lower=0>[2*Tf] sigma_Tf; //fixed effect variance

  //Random effects
  vector[2*Tr] beta_Tr_G_Nr[G,Nr]
  vector[2*Tr] mu_Tr; //random effect anchors

  cholesky_factor_cov[2*Tr] Sigma_Tr_G[G] //variance components
 
}


model {
  real sigma_Tf = 1;
  real gamaC = 1.5;
  real gambC = .5;
 //print("X dim", dims(x));	
 //print("betaD dim", dims(betaD));
  mu_Tf ~ normal(0, 4);
  mu_Tr ~ normal(0, 4);
  /* gamaC ~ student_t(6, 1, 2);
  gambC ~ student_t(6, 1, 2); */
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
      beta_Tf_G[g] ~ normal(mu_Tf, sigma_Tf);
      etaD = x*beta_Tf_G[g][id_Tf];
      etaC = (x[Ipos[iposStart:iposEnd],:]*beta_Tf_G[g][ic_Tf]);
      for(r in 1:Nr){
      	    int rposStart;
	    int rposEnd;
      	    rposStart = rr[r];
	    rposEnd = rr[r+1]-1;
      	    beta_Tr_G_nr[r, g] ~ multi_normal_cholesky(0, Sigma_Tr_G[g]);
      	    etaD[rposStart:rposEnd] = etaD[rposStart:rposEnd] + x[rposStart:rposEnd,:]*beta_Tr_G_nr[r,g,id_Tr]
      }
      if(marginal==1){
	etaC = etaC .* (1+exp(-etaD[Ipos[iposStart:iposEnd]]));
      }
      v[g] ~ bernoulli_logit(etaD);
      sigmaC_G[g] ~ gamma(gamaC, gambC);
      //print("len(y)= ", num_elements(y[iposStart:iposEnd]), "; len(etaC)=", num_elements(etaC[Ipos[iposStart:iposEnd]]), "; len(sigma)=", num_elements(sigmaC_G[g]));
      //target += normal_lpdf( (y[iposStart:iposEnd]) | (etaC[Ipos[iposStart:iposEnd]]) , sigmaC_G[g]);
      y[iposStart:iposEnd] ~ normal(etaC, sigmaC_G[g]);
      }
}

