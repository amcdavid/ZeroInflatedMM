context('Simulation')
library(Matrix)


sng = simulateOne(X, beta, r_dm, nr)

test_that('Simulate from one gene, scalar RF', {
    expect_equal(nrow(sng$Y), nrow(X))
    expect_equal(nrow(sng$b), length(r_dm$Z_list))
})

test_that('Simulate from one gene, vector RF', {
    r_dm = RandomDesignMatrix(X, f = clusters)  
    nr = RanefLikelihood(normal_ranef, Sigma = diag(nrow = 2))
    sng = simulateOne(X, beta, r_dm, nr, obs_likelihood = ObsLikelihood(nb_likelihood))
    expect_equal(length(sng$Y), nrow(X))
    expect_equal(nrow(sng$b), length(r_dm$Z_list))
    expect_false(any(is.na(sng$Y)))
})



test_that('Simulate from many genes, scalar RF', {
    sm = simulateMany(500,  X,
                      GenerateBeta(constant_beta, constant_beta_control(beta0=c(1, 1))),
                      r_dm,
                      GenerateSigma(constant_sigma, constant_sigma_control(Sigma0=3)),
                      normal_ranef,
                      phi_constant(1), 
                      ObsLikelihood(normal_likelihood))
    expect_equal(dim(sm$Y), c(nrow(X), 500))
    expect_true(all(sapply(sm$sigma_list, '==', 3)))
    expect_false(any(is.na(sm$Y)))
})