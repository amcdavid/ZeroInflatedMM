#' @importFrom SummarizedExperiment colData rowData assay assays
stan_frs <- function(obj, model, contrasts, stanargs = list(iter=1e3, chains=4, cores=4), modelargs= list(debug = 0, marginal = 0, ranefs = 2), method='hmc'){
    ## Need to fix zero-expression columns/groups
    obj <- obj[MAST::freq(obj)>0,]
    design_block <-  fixed_design_and_re(colData(obj), model)
    design <- design_block$design
    block <- design_block$block
    v <- 1*(assay(obj)>0)
    G <- nrow(obj)
    IposArr <- apply(v==1, 1, which)
    Ipos <- unlist(IposArr)
    IposGI <- c(1, cumsum(sapply(IposArr, length))+1)
    y <- t(assay(obj))[t(v)==1]
    xr <- design_block$r_design
    
    ## Left closed, right open intervals
    rr <- c(findInterval(levels(block), block, left.open = TRUE), length(block))+1
    ## Iposdt <- data.table(which(v==1, arr.ind=TRUE))
    ## setnames(Iposdt, c('i', 'j'))
    ## block_prod <- CJ(i=seq_len(G), j=rr[-length(rr)])
    ## Iposdt[block_prod,,on=c('i', 'j')]
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
    stanargs <- c(stanargs, list(file=mfile, data=standat))
    if(method=='hmc'){
        do.call(rstan::stan, stanargs)
    } else if(method=='vb'){
        object <- stan_model(mfile)
        rstan::vb(object, data=standat)
    }
    #X# rstan::stan(mfile, data=standat, chains=1, cores=1, iter=2e2)
}
