#'@include ObsLikelihood.R

#' Generator for random effects
#'
#' @param fun random distribution, function of sample size n alone
#' @param Sigma variance-covariance matrix or similar
#' @param args additional arguments passed to \code{fun}
#'
#' @return A function that samples from fun with Sigma fixed
#' @export
#'
#' @examples
#' nr = RanefLikelihood(normal_ranef, Sigma=1)
#' nr(2)
RanefLikelihood = function(fun, Sigma, args=list()){
    assert_that(has_args(fun, c('args')))
    args = do.call(fun$control, args)
    if (length(Sigma) > 1 && !is.matrix(Sigma)) Sigma = diag(Sigma)
    Sigma = as.matrix(Sigma)
    fargs = function(n){
        fun(n, Sigma, args)
    }
    class(fargs)  <- c('function', 'argumentative')
    return(fargs)
}

#' @export
#' @describeIn RanefLikelihood normally distributed random effects
normal_ranef = function(n, Sigma, args=normal_ranef_control()){
    force(args)
    
    mvtnorm::rmvt(n = n, sigma = Sigma, df=args$df)
}

normal_ranef_control = function(df=Inf){
    llist(df)
}

normal_ranef = FunPair(normal_ranef, normal_ranef_control)


RandomDesignMatrix = function(Z, Z_list, f, cluster_idx){
    assert_that(missing(Z) & !missing(Z_list) | !missing(Z) & missing(Z_list))
    assert_that(missing(f) & !missing(cluster_idx) | !missing(f) & missing(cluster_idx))
    
    if(!missing(f)) cluster_idx = split(seq_along(f), f)
    if(missing(Z_list)){
        Z_list = split.data.frame(Z, f)
    }
    retval = llist(cluster_idx, Z_list)
    class(retval) = 'RandomDM'    
    retval
}

fixed_design_and_re <- function(mf, model){
    mf = as.data.frame(mf)
    lf <- lme4::lFormula(model, data=mf, control=lme4::lmerControl(check.formula.LHS = "ignore"))
    design <- lf$X
    assert_that(length(lf$reTrms$flist)==1)
    # we could/should use reTerms for the random design matrix; this would require using sparse matrices in stan
    block <- as.factor(lf$reTrms$flist[[1]])
    bars <- lme4::findbars(model)
    r_form <- as.formula(paste0('~', deparse(lme4::subbars(bars[[1]][[2]]))))
    r_design <- model.matrix(r_form, mf)
    list(design=design, block=block, r_design=r_design)
}
