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
#' @importFrom rstan summary
fit_FxCHM <- function(obj, model, contrasts, stan_control = fxc_stan_control(), model_control= fxc_model_control(), method = c('hmc', 'vb'), returnfit = FALSE, data0, data1){
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
    rr = design_block$rr
    
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
    fit = FxCHM(fit, design_block, model_control, standat, model = model)
    t2 <- proc.time()
    st <- t2-t1
    if(missing(data0)) return(fit)
    if(missing(data1)){
        retval = predict(fit, newdata = data0, type = 'marginal')
    } else {
        retval = get_marginal(fit, data0, data1, ni = min(100, fit@dims$Ns))
    }

    # Figure out how to get a bayesian FDR

    ret = FittedRanefScalar(fixef = retval$fixef, fixef_se=retval$fixef_se,
                      fdr = rep(1, nrow(retval)),
                      sd = retval$sd,
                      walltime = st['elapsed'],
                      coretime = st['user.self'],
                      method = method,
                      obs_used = seq_len(ncol(obj)))
    if(returnfit) return(structure(ret, fit = fit)) else ret
}

setOldClass('family')
setClass('FxCHM', contains = 'stanfit', 
         slots = list(design_re = 'list', 
                      ismarginal = 'logical', 
                      family = 'family',
                      family_zero = 'family',
                      model = 'formula',
                      converged = 'logical',
                      fixef_s_G_Tf_comp = 'array',
                      ranef_s_Nr_G_comp_Tr = 'array',
                      dims = 'list'
                      ))


##' @importFrom stats gaussian binomial
FxCHM = function(stanfit, design_re, model_control, standat, family = gaussian(), family_zero = binomial(), model){
    ismarginal = model_control$pass_as_data$marginal==1
    summ = rstan::summary(stanfit)$summary
    if(is_sampling(stanfit)){
        params = do.call(rbind, get_sampler_params(stanfit, inc_warmup = FALSE))
        rhat = summ[,'Rhat']
        #maybe we want to do this by gene?
        converged = all(params[, 'divergent__']==0) & all(rhat[names(rhat) %like% 'beta_Tf_G'] < 1.1)
    } else{
        converged = NA
    }
    
    dims = standat[c('G', 'Nr', 'Tf', 'Tr')]
    # unpack the coefficients back into an array
    fixef_s_G_Tf_comp = as.matrix(stanfit, par = 'beta_Tf_G') %>% array(., dim = c(nrow(.), dims$G, dims$Tf, 2), dimnames = list(draw = seq_len(nrow(.)), G = seq_len(dims$G), Tf = seq_len(dims$Tf), comp = c('D', 'C')))
    ranef_s_Nr_G_comp_Tr = as.matrix(stanfit, 'beta_Tr_G_Nr') %>% array(., dim = c(nrow(.), dims$Nr, dims$G, 2, dims$Tr), dimnames = list(draw = seq_len(nrow(.)), Nr = seq_len(dims$Nr), G = seq_len(dims$G), comp = c('D', 'C'), Tr = seq_len(dims$Tr)))
    dims$Ns = dim(fixef_s_G_Tf_comp)[1]
    new('FxCHM', stanfit, design_re = design_re, ismarginal = ismarginal, family = family, family_zero = family_zero, model = model, converged = converged, fixef_s_G_Tf_comp = fixef_s_G_Tf_comp,  ranef_s_Nr_G_comp_Tr= ranef_s_Nr_G_comp_Tr, dims = dims)
}

is_sampling = function(fit){
    fit@stan_args[[1]]$method=='sampling'   
}

getDims = function(fit){
    fit@Dims
}

fit = function(x){
    if(inherits(x, 'FxCHM')) return(x)
    y = attr(x, 'fit')
    assert_that(!is.null(y), msg = 'No fit attribute')
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

Drop = function (x, d) {
    dim(x) <- dim(x)[-d]
    x
}

marg_predict = function(object, newdata, re.form, ...){
    pc = predict(object, newdata = newdata, component = 'C', type = 'response', re.form = re.form, ...)
    pd = predict(object, newdata = newdata, component = 'D', type = 'response', re.form = re.form, ...)
    return(pc *pd)
}



predict_FxCHM = function(object, newdata = NULL, type = c('link', 'response', 'marginal'), component = c('C', 'D'), re.form = NULL, index = 0, ...){
    if(!is.null(newdata)){
        design_block =  fixed_design_and_re(newdata, object@model, match_blocks = object@design_re$block)
    } else{
        design_block = object@design_re
    }
    type = match.arg(type, c('link', 'response', 'marginal'))
    if(type == 'marginal') return(marg_predict(object, newdata, re.form, index = index, ...))
    component = match.arg(component, choices = c("C", "D"))
    # get the means in the desired form
    Ns = object@dims$Ns
    if(index < 1 || index > Ns) index = seq_len(Ns)
    
    beta_G_Tf = colMeans(Drop(object@fixef_s_G_Tf_comp[index,,,component,drop=FALSE], 4), dims = 1)
    G = object@dims$G
    N = nrow(design_block$design)
    eta = matrix(NA, nrow = N, ncol = G)
    
    if(is.null(re.form)){
        beta_Nr_G_Tr = colMeans(Drop(object@ranef_s_Nr_G_comp_Tr[index,,,component,,drop=FALSE], 4), dims = 1)
        block = design_block$block
        rr = design_block$rr
        xr = design_block$r_design
    }
    
    for(g in seq_len(G)){
        eta[,g] = design_block$design %*% beta_G_Tf[g,]
        if(is.null(re.form)){
            for(ri in seq_along(levels(design_block$block))){
                if(rr[ri]< (rr[ri+1]-1)){
                    idx = seq(rr[ri], rr[ri+1]-1)
                    eta[idx,g] = eta[idx,g] + xr[idx,,drop=FALSE] %*%beta_Nr_G_Tr[ri,g,]
                }
            }
        }
    }
    if(type == 'response' & component == 'D'){
        return(object@family_zero$linkinv(eta))
    }
    return(eta)
}

setMethod('predict', 'FxCHM', predict_FxCHM)

compare_stanfits = function(fitx, fity){
    sx = summary(fitx)$summary %>% melt()
    sy = summary(fity)$summary %>% melt()
    inner_join(sx, sy, by = c('Var1', 'Var2'))
}

get_marginal = function(obj, data0, data1, ni = 1){
    pred = array(NA, dim = c(nrow(data1), obj@dims$G, ni))
    for(i in seq_len(ni)){
    if(ni==1) index = 0  else index = i
    pred[,,i] = predict(obj, newdata =data1, type = 'marginal', index = index) -  predict(obj, newdata =data0, type = 'marginal', index = index)
    }
    #averaged change
    meaneffect_sim = apply(pred, 2:3, mean)
    meaneffect = rowMeans(meaneffect_sim)
    
    se_mean = apply(pred, 2, sd)
    cluster_sd = apply((aperm(pred, c(2, 1, 3)) - meaneffect)^2
                        ,1, mean)^.5
    #should be equal to zero if the array broadcasting is working
    assert_that(all(abs(apply((aperm(pred, c(2, 1, 3)) - meaneffect) ,1, mean))<1e-4))
    
    data.table(fixef = meaneffect, fixef_se = se_mean, sd = cluster_sd)
}
