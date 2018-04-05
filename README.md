# ZeroInflatedMM
Multi-level models and simulations for zero-inflated data

## (1) A multi-level model that shares information about a variance component across genes

At the moment, this *seems* to work for the traditional hurdle model $E(U|V = 1); \logit E(V = 1)$, but has MCMC convergence issues with the marginal parametrization
$E(U); \logit E(V = 1)$.

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/inst/stan/zeroModels.stan.

## (2) A simulation framework for data sampled from zero-inflated mixed models.

This allows simulation of zero-inflated log-normal and negative binomial data with random effects.  The "population" treatment effect for each gene may follow a spike-and-slab model.  Todo would be to use actual data sets to estimate the treatment parameters. Maybe we should just use Splatter here, but it would be good to simulate from the assumed model for (1).

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Simulate.R and 
https://github.com/amcdavid/ZeroInflatedMM/blob/master/tests/testthat/test-Simulate.R



## (3) An evaluation framework for methods fit to SingleCellExperiment objects

At the moment this wraps three methods (Limma + its repeated measure framework, MAST + blmer, and the multi-level model proposed at item 1).  Ultimately aims to check the calibration/coverage of estimated FDR and confidence intervals and evalaute timings and convergence information.  In the case of simulations, sensitivities and MSE can be calculated.  In real data sets, the level of tests can be evaluated.

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Fitter-limma.R and 
https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Fitter-blmer-hurdle.R

