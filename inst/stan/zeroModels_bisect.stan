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
  int<lower=1, upper=N+1> rr[Nr+1]; //starting index of random group (hence must be sorted..)
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

  //Group-start for each block in Ipos
  //And y
  int<lower=1, upper=NposG+1> RIpos[G*(Nr+1)];

  // DATA
  int<lower=0, upper=1> v[G,N]; //indicators
  vector[NposG] y;

  // DEBUG
  int<lower=0, upper=3> debug;

  // Random effects
  int<lower=0, upper=2> ranefs;
}


transformed data {
  // Index into concatenated parameter vectors
  int c_Tr = Tr+1;
  int c_Tf = (Tf+1);
  real<lower=0> gamaC = 1.5;
  real<lower=0> gambC = .1;
  real<lower=0> sigma_Tf = 4;
  vector[2*Tr] zero_Tr = rep_vector(0.0, Tr*2);
}

parameters {
  matrix[2*Tf,G] beta_Tf_G; //fixed effect regression coefficients
  
  vector<lower=0>[G] sigmaC_G; //dispersion
  // real<lower=0> gamaC; //hyper dispersion
  // real<lower=0> gambC;

  vector[2*Tf] mu_Tf; //fixed effect anchors
  //vector<lower=0.01>[2*Tf] sigma_Tf; //fixed effect variance

  //Random effects
  // vector[2*Tr] mu_Tr; //random effect anchors
  //Whitened prior
  matrix[2*Tr,Nr*G] z_Tr_GNr;
  cholesky_factor_corr[2*Tr] Sigma_Tr_G;//[G]; //variance components
  matrix<lower=0>[2*Tr,G] tau_Tr_G;
}

transformed parameters {
    vector[2*Tr] beta_Tr_G_Nr[Nr,G];
    for(g in 1:G){
      for(r in 1:Nr){
	beta_Tr_G_Nr[r,g,:] = (diag_pre_multiply(tau_Tr_G[:,g],Sigma_Tr_G) * z_Tr_GNr[:,r+(g-1)*Nr]);//'
	}
      }
  //Cauchy 0, 2.5 for scale
  /* matrix<lower=0>[2*Tr,G] tau_Tr_G;
  for(g in 1:G){
    for (k in 1:(2*Tr)){
      tau_Tr_G[k,g] = 2.5 * tan(tau_unif_Tr_G[k,g]);
    }
    } */
}


model {
  vector[N] etaD;
  int iposStart;
  int iposEnd;
  int rposStart;
  int rposEnd;
  int ripi;  //index into RIpos, points at current left-open end point
  //rii+1 points at current right-closed end point
  int ipos_ri_start; //index into Ipos/Y (RIpos[ripi])
  int ipos_ri_end; //index into Ipos/Y (RIpos[ripi+1])	   

  ripi = 0;
  if(debug>2){
    print("X dim", dims(x));	
    print("beta dim", dims(beta_Tf_G));
  }
  mu_Tf ~ normal(0, 4);
  Sigma_Tr_G ~ lkj_corr_cholesky(2);
  to_vector(tau_Tr_G) ~ normal(0, 5);
  to_vector(z_Tr_GNr) ~ normal(0,1);
  for(i in 1:Tf){
    beta_Tf_G[i,] ~ normal(mu_Tf[i], sigma_Tf);
  }
  sigmaC_G ~ gamma(gamaC, gambC);
  //mu_Tr ~ normal(0, 4);
  /* gamaC ~ student_t(6, 1, 2);
     gambC ~ student_t(6, 1, 2); */
  // sigmaD_Tf ~ units are between gene variation of effect, g prior might make sense
  // sigmaD_Tr ~ same 
  for (g in 1:G){
    vector[IposGI[g+1]- IposGI[g]] etaC;
    //index into X/V
    int ipos[IposGI[g+1]- IposGI[g]];
    //index for random effects

    // Because local variable must precede all statements, these are over-allocated
    iposStart = IposGI[g];
    iposEnd = IposGI[g+1]-1;
    ipos = Ipos[iposStart:iposEnd];
    if(debug>1){
      print("G=", g, "; iposStart=", iposStart, "; iposEnd=", iposEnd, "; ipos=", Ipos[iposStart:iposEnd]);
      print("dim(beta_Tr_G_Nr)=", dims(beta_Tr_G_Nr));
    }

    etaD = x*beta_Tf_G[1:Tf,g];
    etaC = (x[ipos,:]*beta_Tf_G[c_Tf:(2*Tf),g]);

    //Init random effect
    if(ranefs>0){
      ripi = ripi +1;
      for(r in 1:Nr){
	// to avoid an "aliasing" message, but probably doesn't save any time
	vector[rr[r+1] -rr[r]] eta_tmp;
	rposStart = rr[r];
	rposEnd = rr[r+1]-1;
	ipos_ri_start = RIpos[ripi];
	ipos_ri_end = RIpos[ripi+1]-1;
	if(debug>0){
	  print("ipos_ri_start=", ipos_ri_start, " ; ipos_ri_end=", ipos_ri_end);
	}
	    
	eta_tmp = etaD[rposStart:rposEnd];
	    
	//beta_Tr_G_Nr[r, g] ~ multi_normal_cholesky(zero_Tr, Sigma_Tr_G); //[g]);
	//beta_Tr_G_Nr[r, g] ~ multi_normal_cholesky(zero_Tr, quad_form_diag(Sigma_Tr_G, tau_Tr_G[:,g])); //[g]);
	etaD[rposStart:rposEnd] = eta_tmp + xr[rposStart:rposEnd,:]*beta_Tr_G_Nr[r,g,1:Tr];
	    
	if((ipos_ri_start <= ipos_ri_end) && (ranefs>1)){ //non-empty group
	  // ipos_ri_start-iposStart+1: start point in etaC for this group.
	  eta_tmp[1:(ipos_ri_end-ipos_ri_start+1)] = etaC[(ipos_ri_start-iposStart+1):(ipos_ri_end-iposStart+1)];
	  etaC[(ipos_ri_start-iposStart+1):(ipos_ri_end-iposStart+1)] =  xr[Ipos[ipos_ri_start:ipos_ri_end],]*beta_Tr_G_Nr[r,g,c_Tr:(2*Tr)];
	}
	ripi = ripi + 1;
      }
    }
    if(marginal==1){
      etaC = etaC .* (1+exp(-etaD[ipos]));
    }
    v[g] ~ bernoulli_logit(etaD);
    // print("len(y)= ", num_elements(y[iposStart:iposEnd]), "; len(etaC)=", num_elements(etaC[Ipos[iposStart:iposEnd]]), "; len(sigma)=", num_elements(sigmaC_G[g]));
    //target += normal_lpdf( (y[iposStart:iposEnd]) | (etaC[Ipos[iposStart:iposEnd]]) , sigmaC_G[g]);
    y[iposStart:iposEnd] ~ normal(etaC, sigmaC_G[g]);
  }
}
