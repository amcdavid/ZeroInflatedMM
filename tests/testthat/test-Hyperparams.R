context('Hyperparameters')

test_that('Sigma Hyperparameter generator', {
    ss = GenerateSigma(constant_sigma, constant_sigma_control(Sigma0=1))
    expect_is(ss, 'function')
    expect_is(ss, 'argumentative')
    expect_equal(ss(), list(Sigma = 1, args = list()))
})
