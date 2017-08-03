##'@import biobroom
limma_frs <- function(obj, model, contrasts){
    design_block <-  fixed_design_and_re(colData(obj), model)
    design <- design_block$design
    block <- design_block$block
    dge <- edgeR::DGEList(counts=(2^assay(obj)-1), samples=as.data.frame(colData(obj)), genes=as.data.frame(mcols(obj)))
    tic = proc.time()
    dge <- edgeR::calcNormFactors(dge)
    v <- limma::voom(dge, design, plot=FALSE, lib.size=dge$samples$lib.size)
    dupcor <- limma::duplicateCorrelation(v,design,block=block)
    y <- limma::voom(dge,design,plot=FALSE,block=block,correlation=dupcor$consensus, lib.size=v$samples$lib.size)
    fit <- limma::lmFit(y,design,block=block,correlation=dupcor$consensus.correlation)
    contrasts <- limma::makeContrasts(contrasts=contrasts, levels=design)
    fit2 <- limma::contrasts.fit(fit, contrasts)
    fit2 <- limma::eBayes(fit2, trend=TRUE)
    sd_group <- sqrt(dupcor$consensus.correlation * fit2$sigma^2)
    toc = proc.time()
    st = toc - tic
    fit_tab <- tidy(fit2) %>% mutate(fixef_se = estimate/statistic, sd = sd_group)
    with(fit_tab, FittedRanefScalar(fixef = estimate,
                                    fixef_se = fixef_se,
                                    fdr = p.adjust(p.value, method='fdr'),
                                    sd = sd,
                                    walltime = st['elapsed'],
                                    coretime = st['user.self'],
                                    method = 'limma',
                                    obs_used = seq_len(ncol(obj))))
}