llist = function(..., labels = FALSE) {
    dotlist <- list(...)
    lname <- names(dotlist)
    name <- vname <- as.character(sys.call())[-1]
    for (i in 1:length(dotlist)) {
        vname[i] <- if (length(lname) && lname[i] != "") 
            lname[i]
        else name[i]
        lab <- vname[i]
        if (labels) {
            lab <- attr(dotlist[[i]], "label", exact = TRUE)
            if (length(lab) == 0) 
                lab <- vname[i]
            Hmisc::label(dotlist[[i]]) <- lab
        }
    }
    names(dotlist) <- vname[1:length(dotlist)]
    dotlist
}


ObsLikelihood = function(fun, args=list()){
    assert_that(has_args(fun, c('eta', 'sigma', 'args')))
    args = do.call(fun$control, args)
    fargs = function(eta, sigma){
        fun(eta, sigma, args)
    }
    class(fargs) = c('function', 'argumentative')
    return(fargs)
}

'$.argumentative' = function(x, name){
    environment(x)[[name]]
}

# Bundle a function with its default arguments
# This allows selective over-riding of arguments in the function wrapper
null_fun = function(){
    return(NULL)   
}

FunPair = function(fun, control = null_fun, ...){
    s = structure(fun, control = control, ...)
    class(s) = c('function', 'FunPair')
    s
}

'$.FunPair' = function(x, name){
    assert_that(name %in% names(attributes(x)))
    attr(x, name, exact = TRUE)
}

normal_likelihood = function(eta, sigma, args){
    eta + rnorm(length(eta))*sigma
}

normal_likelihood = FunPair(normal_likelihood)

logit = function(x) log(x/(1 - x))
expit = function(x) exp(x)/(1 + exp(x))


zin_control = function(marginal=TRUE, b=.2, pi=.5, a=1){
    llist(marginal, b, pi, a)
}

##' @importFrom stats runif rnorm rnbinom
zin_likelihood = function(eta, sigma, args = zin_control()){
    eta_marginal = logit(args$pi)
    etaD = eta * args$b
    etaD = etaD - mean(etaD)*args$a + eta_marginal
    piD = expit(etaD)
    if(args$marginal){
        etaC = eta/piD #inflate continuous mean by prob of zero
    } else{
        etaC = eta   
    }
    U <- (runif(length(eta)) < piD)*1
    U[U > 0] = rnorm(sum(U))*sigma + etaC[U > 0]
    pmax(U, 0)
}

zin_likelihood = FunPair(zin_likelihood, zin_control)

nb_likelihood = function(eta, sigma, args){
    if (args$link == 'log') mu = exp(eta) else mu = pmax(0, eta)
    stats::rnbinom(length(eta), size = 15/sigma, mu = mu)
}

nb_likelihood_control = function(link = 'identity'){
    link = match.arg(link, c('identity', 'log'))
    llist(link)
}

nb_likelihood = FunPair(nb_likelihood, nb_likelihood_control)

## Dispersion function
phi_constant = function(sigma=1){
    fun = function(mu, X) sigma
    return(fun)
}