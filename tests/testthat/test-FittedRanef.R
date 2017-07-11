library(testthat)
context('Test metrics')

ngene = 5
neff = 1
fixef = -2:2
pid = seq_len(ngene*neff)

frf = FittedRanefScalar(fixef, fixef_se = rep(0, ngene*neff), sd = rep(1, ngene*neff), fdr_q=c(0, 0, 0, 1, 1), primerid = pid, method='test', walltime=1, coretime=4)

truthrf = TruthRanefScalar(fixef, sd = rep(1, ngene*neff), primerid = pid)

test_that('Can create FitteRanef', {
    expect_is(frf, 'FittedRanefScalar')
})

test_that('Can create TruthRanef', {
    expect_is(truthrf, 'TruthRanefScalar')
})

numeric.equal = function(x, y) all(abs(x - y) < 1e-5)

test_that('Check power', {
    power = check_fdr_power(frf, truthrf)  
    expect_true(numeric.equal(power$obs_power,.5))
    expect_true(numeric.equal(power$obs_fdr, 1/3))
})

test_that('Check bias', {
    eff = check_bias(frf, truthrf)
    expect_true(numeric.equal(eff$err_fixef, 0))
    expect_true(numeric.equal(eff$err_sd, 0))
    truthrf$sd[] = 0
    eff2 = check_bias(frf, truthrf)
    expect_true(numeric.equal(eff2$err_sd, frf$sd))
})
