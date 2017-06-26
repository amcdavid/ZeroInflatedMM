#' @include ObsLikelihood.R

#Fixed effects
GenerateBeta = function(fun, args=list()){
    assert_that(has_args(fun, c('i', 'args')))
    args = do.call(fun$control, args)
    fargs = function(i){
        fun(i, args)
    }
    class(fargs)  <- c('function', 'argumentative')
    return(fargs)
}

# constant beta
constant_beta = function(i, args = constant_beta_control()){
 force(args)
 return(args$beta0)
}

constant_beta_control = function(beta0=1){
    llist(beta0)   
}

constant_beta = FunPair(constant_beta, constant_beta_control)

# spike-slab
spike_slab_beta = function(i, args = spike_slab_beta_control()){
    V = runif(length(args$nullcoords)) < args$pi0 # V==1 if null
    beta0 = as.vector(args$beta0 + mvtnorm::rmvt(1, sigma = args$Sigma_beta, df = args$df))
    beta0[args$nullcoords][V == 1] = 0
    beta0
}

spike_slab_beta_control = function(pi0=.9, nullcoords = 2, Sigma_beta=diag(nrow = 2), beta0=rep(1, 2), df=Inf){
    stopifnot(ifelse(is.matrix(Sigma_beta), nrow(Sigma_beta) == length(beta0), length(Sigma_beta) == length(beta0)))
    llist(pi0, nullcoords, Sigma_beta, beta0, df)
}

spike_slab_beta = FunPair(spike_slab_beta, spike_slab_beta_control)

# Random effects
GenerateSigma = function(fun, args=list()){
    assert_that(has_args(fun, c('i', 'beta', 'args')))
    args = do.call(fun$control, args)
    fargs = function(i, beta){
        ff = fun(i, beta, args)
        if (is.null(ff$args)) ff$args = list()
        ff
    }
    class(fargs)  <- c('function', 'argumentative')
    return(fargs)
}

# constant sigma
constant_sigma = function(i, beta, args = constant_sigma_control()){
    return(list(Sigma = args$Sigma0))
}

constant_sigma_control = function(Sigma0=1){
    list(Sigma0 = Sigma0)   
}

constant_sigma = FunPair(constant_sigma, constant_sigma_control)

# log-normal sigma
ln_sigma =  function(i, beta, args = ln_sigma_control()){
    list(Sigma = exp(rnorm(length(args$Sigma0_log), mean = args$Sigma0_log, sd = args$Sigma0_Sigma_log)))
}

ln_sigma_control = function(Sigma0=1, Sigma0_Sigma=1, Sigma0_log, Sigma0_Sigma_log){
    if (missing(Sigma0_log)) {
        cv = Sigma0_Sigma/Sigma0
        Sigma0_log = log(Sigma0/sqrt(1 + cv^2))
        Sigma0_Sigma_log = sqrt(log(1 + cv^2))
    }
    llist(Sigma0_log, Sigma0_Sigma_log)
}

ln_sigma = FunPair(ln_sigma, ln_sigma_control)
