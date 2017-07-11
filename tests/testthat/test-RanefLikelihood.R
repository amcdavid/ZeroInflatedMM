context('RanefLikelihood')
test_that('Can pull out random effects and model matrices', {
    mf <- data.frame(g=clusters, X)
    parsed <- fixed_design_and_re(mf, ~treat + (1|g))
    expect_equal(parsed$block, clusters)
    mmtreat <- model.matrix(~treat, mf)
    expect_equivalent(parsed$design, mmtreat)
    mm1 <- model.matrix(~1, mf)
    expect_equal(parsed$r_design, mm1)
})

test_that('Handle no intercept model', {
    mf <- data.frame(g=clusters, X)
    parsed <- fixed_design_and_re(mf, ~0+(1|g))
    expect_equivalent(parsed$design, model.matrix(~0, mf))
})

test_that('Handle scrambled covariates', {
    mf <- data.frame(g=clusters, X)
    mf <- mf[sample(nrow(mf)),]
    parsed <- fixed_design_and_re(mf, ~ (1|g))
    expect_equivalent(parsed$block, mf$g)
})
