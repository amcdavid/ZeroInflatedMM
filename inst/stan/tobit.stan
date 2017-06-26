data {
  int<lower=0> N; //number of obs
  int<lower=0> Tf; //number of fixed covariates
  int<lower=0> Tr; //number of random covariates
  int<lower=0> Nr; //number of random groups
  int<lower=1, upper=Nr> rr[N]; //observation -> group mapping
  int<lower=0> G; //number of genes

  matrix[N, Tr] xr; //random effects covariate matrix, fixed across genes
  vector[N] y[G]; //zero-inflated data
  int<lower=0, upper=1> v[G,N]; //indicators
}

parameters {
  vector[Tf] betaD_Tf_G[G]; //fixed effect regression coefficients
  vector[Tr] betaD_Tr_G_Nr[G,Nr]; //random effect coefficients
  
  vector<lower=0>[Tr] kappa_Tr_G[G]; //random effect variances
  vector[Tr] muD_Tr_G[G]; //random means

  vector[Tf] muD_Tf; //fixed effect anchors
  vector[Tr] muD_Tr; //random effect anchors
  vector<lower=0>[Tf] sigmaD_Tf; //fixed effect variance
  vector<lower=0>[Tr] sigmaD_Tr; //random effect variance of mean
  vector<lower=0>[Tr] gamaD_Tr; //random effect 
  vector<lower=0>[Tr] gambD_Tr; //random
}

model {
      /*
      rho0 ~ 2*(beta(2, 2)-.5);
      for(g in 1:G){
  	sigma[G][1,2] <- rho0;
	sigma[G][2,1] <- rho0;
	//Ystar[G] ~ mv_norm(eta[
	
	//y[N,G] <- if_else(Ystar[G][1]>0, Ystar[G][2], 0)
  }
  
  mu_alpha ~ normal(0, 100);
  mu_beta ~ normal(0, 100);
  sigmasq_y ~ inv_gamma(0.001, 0.001);
  sigmasq_alpha ~ inv_gamma(0.001, 0.001);
  sigmasq_beta ~ inv_gamma(0.001, 0.001);
  alpha ~ normal(mu_alpha, sigma_alpha); // vectorized
  beta ~ normal(mu_beta, sigma_beta);  // vectorized
*/
 //print("X dim", dims(x));	
 //print("betaD dim", dims(betaD));

  muD_Tr ~ normal(0, 4);
  muD_Tf ~ normal(0, 4);
  sigmaD_Tf ~ student_t(6, 1, 2); //lognormal(1, 2);
  sigmaD_Tr ~ student_t(6, 1, 2);
  gamaD_Tr ~ student_t(6, 1, 2);
  gambD_Tr ~ student_t(6, 1, 2);
  // sigmaD_Tf ~ units are between gene variation of effect, g prior might make sense
  // sigmaD_Tr ~ same 
  for (g in 1:G){
      vector[N] etaD;	
      muD_Tr_G[g] ~ normal(muD_Tr, sigmaD_Tr);
      kappa_Tr_G[g] ~ student_t(6, gamaD_Tr, gambD_Tr);
      betaD_Tf_G[g] ~ normal(muD_Tf, sigmaD_Tf);
      for (r in 1:Nr){
      	  betaD_Tr_G_Nr[g,r] ~ normal(muD_Tr_G[g], kappa_Tr_G[g]);
	}
      etaD <- x*betaD_Tf_G[g]; 
      for (i in 1:N){
      	   etaD[i] <- etaD[i] + xr[i]*betaD_Tr_G_Nr[g, rr[i]];
      }
      //betaD_Tr_G_Nr ~ normal(muD_Tr_G[g], kappa_Tr_G[g]);
      v[g] ~ bernoulli_logit(etaD);
      //y ~ normal(etaC[g], sigma_y);      	  
      }

  }

/*
generated quantities {
  real alpha0;
  alpha0 <- mu_alpha - xbar * mu_beta;
}
*/


/*
y ~ f(mu[l,g,n])
mu[l,g,n] <- betaf[g]*x[l,g,n]+betar[l,g]*xr[l,g,n]
betaf[g] ~ N(beta0, sigma_0) iid over k,g //could be student t
betar[l,g] ~ N(gamma[g], sigma_g) iid over ...?
gamma[g] ~ N(gamma0, tau_0) iid over k,g
sigma_g ~ Gamma(a, b) iid

*/