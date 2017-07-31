fxc_model_control = function(debug = 0, marginal = 1, ranefs = 2, center = TRUE){
    list(pass_as_data = llist(debug, marginal, ranefs), other = llist(center))    
}

fxc_stan_control = function(iter = 1e3, chains = 4, cores = 4, ...){
    llist(iter, chains, cores, ...)   
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
    obj <- obj[MAST::freq(obj)>0,]
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
    rr <- c(findInterval(levels(block), block, left.open = TRUE), length(block))+1
    RIpos <- list()
    for (i in seq_len(G)){
        slice <- IposArr[[i]]
        RIpos[[i]] <- findInterval(rr, slice, left.open=TRUE)+1
        if(i>1) RIpos[[i]] <- RIpos[[i]]+ RIpos[[i-1]][length(rr)]-1
    }
    #Left closed, right open endpoints of y and Ipos for each group
    RIpos <- unlist(RIpos)
    assert_that(all(sort(RIpos) == RIpos), #non decreasing
                length(RIpos)==(nlevels(block)+1)*G, #length of G*(B+1)
                RIpos[length(RIpos)]==(length(Ipos)+1)) #last endpoint is right open
    standat <- list(N=ncol(obj), Tf=ncol(design), G=G,
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
    }
    
    contr_fixef <- which(colnames(design) %like% contrasts)
    contr_ranef <- which(colnames(xr) %like% contrasts)
    fixef <- makeParTable(fit, 'beta_Tf_G')
    ranef <- makeParTable(fit, 'tau_Tr_G')
    t2 <- proc.time()
    st <- t2-t1
    browser()
    # Figure out how to get a bayesian FDR
    bound <- fixef[j %in% c(2, 4),.(pval=pnorm(abs(Zrob), lower.tail=FALSE)*2), keyby=list(i,j)]
    bound[,fdr:=p.adjust(pval, method='fdr')]
    bound <- bound[,.(minfdr=min(fdr)),keyby=i]
    FittedRanefScalar(fixef = fixef[j==2,mean], fixef_se=fixef[j==2,sd],
                      fdr =bound[,minfdr],
                      sd = ranef[i==1, mean],
                      walltime = st['elapsed'],
                      coretime = st['user.self'],
                      method = method,
                      obs_used = seq_len(ncol(obj)))
}

setClass('FxCHM', contains = 'stanfit', 
         slots = list(design = 'data.frame', 
                      ismarginal = 'logical', 
                      family = 'family',
                      family_zero = 'family',
                      formula = 'formula'))

##' @importFrom stats gaussian binomial
FxCHM = function(stanfit, design, ismarginal, family = gaussian(), family_zero = binomial(), formula){
    new('FxCHM', stanfit, design=design, ismarginal = ismarginal, family = family, family_zero = family_zero)
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
