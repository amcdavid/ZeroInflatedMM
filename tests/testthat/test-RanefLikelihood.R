context('RanefLikelihood')
test_that('Can pull out random effects and model matrices', {
    mf <- data.frame(g=clusters, X)
    parsed <- fixed_design_and_re(mf, ~treat + 1|g)
    expect_equal(parsed$block, cluster)
    expect_equal(parsed$design, model.matrix(~treat, mf))
    expect_equal(parsed$r_design, model.matrix(~1, mf))
})

test_that('Handle no intercept model', {
    mf <- data.frame(g=clusters, X)
    parsed <- fixed_design_and_re(mf, ~0+1|g)
    expect_equal(parsed$design, model.matrix(~0, mf))
})

test_that('Handle scrambled covariates', {
    mf <- data.frame(g=clusters, X)
    mf <- mf[sample(nrow(mf)),]
    parsed <- fixed_design_and_re(mf, ~1|g)
    expect_equal(parsed$block, mf$g)
})
