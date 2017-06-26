#    # Slots #
# 
# ## fitted models ##
# -point estimate for pop mean, SD
# -standard error for pop mean, or CI?
# -point estimate for SD of pop mean, or CI?
# -decision function at a level (FDR?)
# -core-time and wall-time
# -which observations were used
# 
# ## Truth objects ##
# -true pop means
# -true sd of pop means
# -families/mispecification/etc?

#' Construct a Fitted random effect object for a scalar parameter of interest
#'
#' @param fixef  vector of fixed effects, eg one per gene
#' @param fixef_se standard errors for fixed effects
#' @param sd vector of population variance components
#' @param fdr_q FDR q-values
#' @param method 
#' @param walltime 
#' @param coretime 
#' @param obs_used 
#'
#' @return list
#' @export
#' @import assertthat
FittedRanefScalar  = function(fixef, fixef_se, sd, fdr_q, method, walltime, coretime, obs_used=NA){
    assert_that(is.vector(fixef), is.numeric(fixef))
    assert_that(length(fixef_se) == length(fixef), is.numeric(fixef_se))
    assert_that(length(sd) == length(fixef), is.numeric(sd))
    assert_that(length(fdr_q) == length(fixef), all(fdr_q <= 1), all(fdr_q >= 0))
    
    obj = llist(fixef, fixef_se, sd, fdr_q, method, walltime, coretime, obs_used)
    class(obj) = c('FittedRanefScalar', 'FittedRanef')
    obj
}

rbind.FittedRanef = function(...){
    alist = list(...)
    if(length(alist)<1) return(NULL)
    out = lapply(seq_along(alist[[1]]), function(i){
           drop(do.call(rbind, lapply(alist, '[[', i)))
    })
    class(out) = class(alist[[1]])
}

TruthRanefScalar = function(fixef, sd, obs_used=NA, other_params=NA){
    obj = llist(fixef, sd, is_sig = abs(fixef) > 0, obs_used, other_params)
    class(obj) = c('TruthRanefScalar', 'TruthRanef')
    obj
}


print.FittedRanefScalar = function(x, ...){
    print(sprintf('%s on %d components, significant on %d', class(x)[1], length(x$fixef), sum(x$fdr_q < .1)), ...)
}

print.TruthRanefScalar = function(x, ...){
    x$fdr_q = 1*(!x$is_sig)
    print.FittedRanefScalar(x, ...)
}

# Should be based on predictive interval
check_calibration = function(ranef0, ranef1){
    calZ = (ranef0$fixef-ranef1$fixef)/ranef0$sd
}

check_fdr_power = function(ranef, truef, expected_fdr = seq(.01, .51, by = .05)){
    assert_that(inherits(ranef, 'FittedRanef'))
    assert_that(inherits(truef, 'TruthRanef'))
    is_sig = outer(ranef$fdr_q, expected_fdr, `<`)
    is_true = apply(is_sig, 2, `&`, truef$is_sig)
    obs_fdr = 1-colSums(is_true)/ #true
        colSums(is_sig) #discoveries
    obs_power = colSums(is_true)/sum(truef$is_sig)
    llist(obs_fdr, obs_power, expected_fdr)
}

check_bias = function(ranef, truef){
    err_fixef = ranef$fixef - truef$fixef
    err_ranef = ranef$sd - truef$sd
    llist(err_fixef, err_ranef)
}

# #Methods#
# -fit
# -time
# -check
