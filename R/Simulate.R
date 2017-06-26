# SimRanef
# 1. Can we partition into nuisance (intercepts, residual variances) vs target parameters (treatment contrasts)?  In some cases treatment contrasts comprise multiple random effects..?
# 2. Can we assert that there is only a single level -> Prior for mean vector is exchangable along particular orbits?  Certain conditional indepedences hold?
# 3. Can we assume scalar random effect?
# 4. How to parameterize intergene effects?

#' Simulate from a single feature/gene
#'
#' @param X fixed effect design matrix
#' @param beta fixed effect coefficients
#' @param random_dm a RandomDesignMatrix object
#' @param ranef_likelihood random effect generator function, object of class RanefLikelihood
#' @param phi_fun dispersion generator function
#' @param obs_likelihood the observation generator, object of class ObsLikelihood
#' @details 
#' @section Generator Functions
#' The generator functions have the following API:
#'
#' @return vector of observations Y, linear predictor eta, random effects b
#' @export
#'
#' @examples
simulateOne = function(X, beta, random_dm, ranef_likelihood, phi_fun = phi_constant(), obs_likelihood=ObsLikelihood(normal_likelihood)){
    cluster_idx = random_dm$cluster_idx
    Z_list = random_dm$Z_list
    #assert_that(are_equal(unlist(cluster_idx), sort(unlist(cluster_idx))))
    assert_that(sum(sapply(Z_list, nrow))==nrow(X))
    assert_that(ncol(X)==length(beta))
    eta = X %*% beta
    b = ranef_likelihood(length(Z_list))
    for(i in seq_along(cluster_idx)){
        ci = cluster_idx[[i]]
        eta[ci] = eta[ci] + Z_list[[i]] %*% b[i,]
    }
    sigma = phi_fun(eta, X)
    Y = obs_likelihood(eta, sigma)
    llist(Y, eta, b)
}

simulateMany = function(Ng=1, X, beta_fun, random_dm, sigma_b_fun, ranef_family, phi_fun, obs_likelihood){
    beta_list = sigma_list = list()
    Y = matrix(NA, nrow=nrow(X), ncol=Ng)
    for(i in seq_len(Ng)){
        beta_list[[i]] = beta = beta_fun(i)
        sigma_b = sigma_b_fun(i, beta)
        sigma_list[[i]] = sigma_b$Sigma
        ranef_likelihood = RanefLikelihood(ranef_family, sigma_b$Sigma, sigma_b$args)
        ##phi_fun = sigma_phi_fun(i, beta)
        Y[,i] = simulateOne(X, beta, random_dm, ranef_likelihood, phi_fun, obs_likelihood)$Y
    }
    llist(Y, beta_list, sigma_list)
}


## Here's what params can look like:
## Sigma0_Sigma sigma_phi Nn Ng likelihood design px pi0 Nc
## 0.1 0.5 5 10 nb between 0.5 1 20
## 0.1 0.5 20 100 zin between 0.1 1 20
makeY = function(params, i){
    likelihood = as.character(params$likelihood[i])
    design = params$design[i]
    pnum = suppressWarnings(lapply(params, function(x) as.numeric(as.character(x[i]))))
    Ntot = pnum$Nc*pnum$Nn
    clusters = gl(pnum$Nn, pnum$Nc)
    if(design == 'within'){
        treat = rbinom(Ntot, size = 1, p = pnum$px)
    } else{
        ntreat = floor(pnum$Nn*pnum$px)
        treat = c(rep(0, ntreat*pnum$Nc), rep(1, (pnum$Nn - ntreat)*pnum$Nc))
    }
    X = cbind(1, treat)
    r_dm = RandomDesignMatrix(X, f=clusters)
    if (likelihood == 'nb'){
        Likelihood = ObsLikelihood(nb_likelihood, nb_likelihood_control(link = 'log'))
        beta00 = log(1)
        scale = .2
        transfo <- function(x) log2(x+1)
    } else if (likelihood == 'zin'){
        Likelihood = ObsLikelihood(zin_likelihood, zin_control(pi = .2, a = .2, b = .5))
        beta00 = 2
        scale = 1
        transfo <- identity
    }
    gen_beta = GenerateBeta(spike_slab_beta, spike_slab_beta_control(pi0 = pnum$pi0, beta0 = c(beta00, 1*scale), df = 15, Sigma_beta = diag(c(1.5, .5))*scale))
    gen_sigma = GenerateSigma(ln_sigma, ln_sigma_control(Sigma0 = c(.8, .8)*scale, Sigma0_Sigma = rep(pnum$Sigma0_Sigma, 1)*scale))
    Y_list = simulateMany(pnum$Ng,  X,
                      gen_beta,
                      r_dm,
                      gen_sigma,
                      normal_ranef,
                      phi_constant(pnum$sigma_phi), 
                      Likelihood)
    sca = FromMatrix(exprs = t(Y_list$Y) %>% transfo, cData = as.data.frame(cbind(X = X[,-1, drop=FALSE], g = clusters)))
    c(Y_list[-1], list(sca=sca, is_logged = TRUE))
}
