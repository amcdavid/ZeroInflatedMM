vcHook <- function(fit){
    rfc <- ranef(fit@fitC)
    rfd <- ranef(fit@fitD)
    class(rfd) <- class(rfc) <- 'list'
    ranefs <- suppressMessages(reshape2::melt(list(C=rfc, D=rfd)))
    dfC <- as.data.frame(VarCorr(fit@fitC))
    dfD <- as.data.frame(VarCorr(fit@fitD))
    varcor <- suppressMessages(reshape2::melt(list(C=dfC, D=dfD)))
    list(ranefs=ranefs, varcor=varcor, converge=fit@optimMsg)
}

blmer_hurdle_frs <- function(obj, model, contrasts, contrasts0, parallel = FALSE){
    st = system.time(zz <- zlm(model, obj, method='blmer', ebayes=FALSE, hook=vcHook, 
                               parallel=parallel, fitArgsD=list(fixef.prior='t'), strictConvergence=FALSE))
    #deg = summary(zz, doLRT=list(Hypothesis('treatmentCIT'), Hypothesis('treatmentRF')))
    converge = setNames(plyr::ldply(zz@hookOut, '[[', 'converge'), c('primerid', 'C', 'D'))
    setDT(converge)
    converge[,aprob := !is.na(D) | !is.na(C)]
    #wt = list(treat = lrTest(zz, CoefficientHypothesis('treat')))
    #wt = as.data.table(melt(wt))
    #wt = wt[test.type=='hurdle' & metric %like% 'Pr',]
    #wt[,p.adj:=p.adjust(value, method='bonferroni'), key=L1]
    s = MAST::summary(zz, doLRT=TRUE)$datatable
    degConverge = merge(s[component %in% c('H', 'logFC'), .(component, primerid, coef, `Pr(>Chisq)`, z, sd_na=NA_real_)], converge, by='primerid')
    fixef <- degConverge[component=='logFC', coef]
    FittedRanefScalar(fixef=fixef,
                      fixef_se = degConverge[component=='logFC', coef/z],
                      sd = degConverge[component=='logFC', sd_na],
                      fdr = p.adjust(degConverge[component=='H', `Pr(>Chisq)`], method='fdr'),
                      walltime = st['elapsed'],
                      coretime = st['user.self'],
                      method = 'blmer_hurdle',
                      obs_used = seq_len(ncol(obj)))
}
