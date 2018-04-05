# ZeroInflatedMM
Multi-level models and simulations for zero-inflated data

## A multi-level model that shares information about a variance component across genes

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/inst/stan/zeroModels.stan.
At the moment, this seems to work for the traditional hurdle model $E(U|V = 1); \logit E(V = 1)$, but has MCMC convergence issues with the marginal parametrization
$E(U); \logit E(V = 1)$.

## A simulation framework for data sampled from zero-inflated mixed models.

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Simulate.R and 
https://github.com/amcdavid/ZeroInflatedMM/blob/master/tests/testthat/test-Simulate.R

## An evaluation framework for methods fit to SingleCellExperiment objects

See https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Fitter-limma.R and 
https://github.com/amcdavid/ZeroInflatedMM/blob/master/R/Fitter-blmer-hurdle.R
