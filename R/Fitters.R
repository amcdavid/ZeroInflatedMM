#' @importFrom SummarizedExperiment colData rowData assay assays
stan_frs <- function(obj, model, contrasts, stanargs = list(iter = 1e3, chains = 4, cores = 4), modelargs= list(debug = 0, marginal = 1, ranefs = 2), method = 'hmc', returnfit = FALSE, center = FALSE){
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
    if(center){
        yraw <- apply(yraw, 2, function(x){
            yi <- which(abs(x)>0)
            x[yi] <- scale(x[yi])
            x
        })
    }
    y <- t(yraw)[t(v)==1]
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
    standat <- c(standat, modelargs)
    mfile <- system.file('stan', 'zeroModels_bisect.stan', package='ZeroInflatedMM')

    if(method=='hmc'){
        stanargs <- c(stanargs, list(file=mfile, data=standat))
        fit <- do.call(rstan::stan, stanargs)
    } else if(method=='vb'){
        object <- stan_model(mfile)
        stanargs <- c(stanargs, list(object=object, data=standat))
        fit <- do.call(rstan::vb, stanargs)
    }
    contr_fixef <- which(colnames(design) %like% contrasts)
    contr_ranef <- which(colnames(xr) %like% contrasts)
    fixef <- makeParTable(fit, 'beta_Tf_G')
    ranef <- makeParTable(fit, 'tau_Tr_G')
    t2 <- proc.time()
    st <- t2-t1
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
