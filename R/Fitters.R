fxc_model_control = function(debug = 0, marginal = 1, ranefs = 2, prior_precision = 1, prior_precision_binomial = 1, center = TRUE){
    list(pass_as_data = llist(debug, marginal, ranefs, prior_precision, prior_precision_binomial), other = llist(center))    
}

fxc_stan_control = function(include = FALSE, pars = c('z_Tr_GNr', 'z_tau_Tr_G', 'L_Sigma_Tr_G', 
                                                      'theta_Tf_G'), ...){
    llist(include, pars, ...)   
}

stan_modelfile = system.file('stan', 'zeroModels_bisect.stan', package='ZeroInflatedMM')


#' Fit a zero-inflated hierarchical mixed model
#' 
#' @param obj SingleCellAssay
#' @param model formula
#' @param contrasts contrasts to test
#' @param stan_control arguments to pass to stan
#' @param model_control arguments to pass to model
#' @param method 
#' @param returnfit return the fit (object of class FxCHM, extending stanfit-class)
#'
#' @importFrom SummarizedExperiment colData rowData assay assays
#' @importClassesFrom rstan stanfit
fit_FxCHM <- function(obj, model, contrasts, stan_control = fxc_stan_control(), model_control= fxc_model_control(), method = c('hmc', 'vb'), returnfit = FALSE){
    ## Extract control args
    control = model_control$other
    ## Need to fix zero-expression columns/groups
    if(any(MAST::freq(obj)==0)){
        warning('Null genes present; dropping')
        obj <- obj[MAST::freq(obj)>0,]
    }
    design_block <-  fixed_design_and_re(colData(obj), model)
    design <- design_block$design
    block <- design_block$block
    t1 <- proc.time()
    v <- 1*(assay(obj)>0)
    G <- nrow(obj)
    IposArr <- apply(v==1, 1, which)
    Ipos <- unlist(IposArr)
    IposGI <- c(1, cumsum(sapply(IposArr, length))+1)
    yraw <- assay(obj)
    if(control$center){
        yraw <- apply(yraw, 1, function(x){
            yi <- which(abs(x)>0)
            if(length(yi)>0){
            x[yi] <- scale(x[yi], scale = length(yi)>1)
            }
            x
        })
    }
    y <- t(yraw)[t(v)==1]
    assert_that(!any(is.na(y)))
    xr <- design_block$r_design
    
    ## Left closed, right open intervals
    rr <- c(findInterval(seq_along(levels(block)), as.numeric(block), left.open = TRUE), length(block))+1
    RIpos <- list()
    for (i in seq_len(G)){
        slice <- IposArr[[i]]
        RIpos[[i]] <- findInterval(rr, slice, left.open=TRUE)+1
        if(i>1) RIpos[[i]] <- RIpos[[i]]+ RIpos[[i-1]][length(rr)]-1
    }
    #Left closed, right open endpoints of y and Ipos for each group
    RIpos <- unlist(RIpos)
    assert_that(all(sort(RIpos) == RIpos), #non decreasing
                length(RIpos) == (nlevels(block) + 1)*G, #length of G*(B+1)
                RIpos[length(RIpos)] == (length(Ipos) + 1)) #last endpoint is right open
    standat <- list(N = ncol(obj), Tf=ncol(design), G=G,
                    NposG=length(Ipos), x=design, v=v,
                    Ipos=Ipos, IposGI=IposGI, RIpos=RIpos,
                    y=y, xr = xr, Nr = nlevels(block), rr = rr,
                    Tr = ncol(xr))
    standat <- c(standat, model_control$pass_as_data)
    if(method=='hmc'){
        stan_control <- c(stan_control, list(file=stan_modelfile, data=standat))
        fit <- do.call(rstan::stan, stan_control)
    } else if(method=='vb'){
        object <- stan_model(stan_modelfile)
        stan_control <- c(stan_control, list(object=object, data=standat))
        fit <- do.call(rstan::vb, stan_control)
    } else if(method=='optimizing'){
        object <- stan_model(stan_modelfile)
        stan_control = stan_control[setdiff(names(stan_control), c('include', 'pars'))]
        stan_control <- c(stan_control, list(object=object, data=standat))
        fit <- do.call(rstan::optimizing, stan_control)
        return(fit)
    }
    
    fixef <- makeParTable(fit, 'beta_Tf_G')
    fixef_interp = expand.grid(jname = colnames(design), comp = c('D', 'C'), stringsAsFactors = FALSE) %>% as.data.table()
    fixef_interp[,j:=.I]
    fixef = dplyr::left_join(fixef, fixef_interp, by = 'j')
    ranef <- makeParTable(fit, 'tau_Tr_G')
    ranef_interp = expand.grid(iname = colnames(xr), comp = c('D', 'C'), stringsAsFactors = FALSE) %>% as.data.table()
    ranef_interp[,i:=.I]
    ranef = left_join(ranef, ranef_interp, by = 'i')
    t2 <- proc.time()
    st <- t2-t1
    # Figure out how to get a bayesian FDR

    fixef_contr = fixef[jname %like% contrasts & comp=='C',]
    if(any(contrasts %in% colnames(xr))) {
      ranef_contr = ranef[iname %like% contrasts & comp == 'D',]
    } else{
      ranef_contr = ranef[iname %like% '(Intercept)']   
    }
    bound = fixef_contr[,.(pval=pnorm(abs(Zrob), lower.tail=FALSE)*2)]
    bound[,fdr:=p.adjust(pval, method='fdr')]
    #bound <- bound[,.(minfdr=min(fdr)),keyby=i]
    ret = FittedRanefScalar(fixef = fixef_contr[,mean], fixef_se=fixef_contr[,sd],
                      fdr =bound[,fdr],
                      sd = ranef_contr[, mean],
                      walltime = st['elapsed'],
                      coretime = st['user.self'],
                      method = method,
                      obs_used = seq_len(ncol(obj)))
    if(returnfit) return(structure(ret, fit = fit)) else ret
}

setOldClass('family')
setClass('FxCHM', contains = 'stanfit', 
         slots = list(design = 'data.frame', 
                      ismarginal = 'logical', 
                      family = 'family',
                      family_zero = 'family',
                      formula = 'formula'))

##' @importFrom stats gaussian binomial
FxCHM = function(stanfit, design, ismarginal, family = gaussian(), family_zero = binomial(), formula){
    new('FxCHM', stanfit, design = design, ismarginal = ismarginal, family = family, family_zero = family_zero)
}

fit = function(x){
    y = attr(x, 'fit')   
    assert_that(!is.null(y))
    return(y)
}

#' @import data.table
makeParTable <- function(fit, pname, rename){
    pdim <- fit@par_dims[[pname]]
    if(missing(rename)) rename <- pname
    if(length(pdim)>0){
        namesMaybe <- c('i', 'j', 'k', 'l')
        if(length(pdim)<5){
            names(pdim) <- namesMaybe[seq_along(pdim)]
        } else{
            stop('Parameter dimension  rank must be < 5')
        }
        pgrid <- do.call(expand.grid, lapply(rev(pdim), seq_len))
        
    }
    nchains <- length(fit@inits)
    ss <- rstan::summary(fit, pars=pname)$summary
    dt <- data.table(pgrid, ss[,c('mean', "sd", "2.5%", "97.5%")], rename)
    dt[,Zrob:=mean/abs((`97.5%`-`2.5%`)*.5)]
    dt[,Z:=mean/sd]
}

nonParEstimate <- function(fit, contrast0, contrast1){
    # get posterior samples of means for each group
    list_of_draws <- extract(fit)
    print(names(list_of_draws))
}

compare_stanfits = function(fitx, fity){
    sx = summary(fitx)$summary %>% melt()
    sy = summary(fity)$summary %>% melt()
    inner_join(sx, sy, by = c('Var1', 'Var2'))
}
