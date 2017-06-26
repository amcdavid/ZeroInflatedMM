context('Likelihoods')

test_likelihood = function(eta, sigma, args){
    return(args)   
}

test_likelihood_control = function(A, B=10){
    llist(A, B)
}

test_likelihood = FunPair(test_likelihood, test_likelihood_control)

test_that('Obs likelihood class behaves', {
    ol = ObsLikelihood(test_likelihood, list(A = 1))
    expect_is(ol, 'function')
    expect_is(ol, 'argumentative')
    expect_equal(ol$args, llist(A = 1, B = 10))
    expect_equal(ol(1, 1), llist(A = 1, B = 10))
})

ol = ObsLikelihood(normal_likelihood)
test_that("Can sample normal", {
  oo = ol(eta = rep(100, 1000), sigma = 1:10)
  tt = t.test(oo, mu = 100)
  expect_lt(abs(tt$statistic), 3)
})

eta_seq = seq(0, 10, length.out = 1000)
test_that('Can sample zero-inflated normal', {
    ol_zin = ObsLikelihood(zin_likelihood, zin_control(b = 0))
    oo = ol_zin(eta = rep(100, 1000), sigma = 1)
    tt = t.test(oo, mu = 100)
    expect_lt(abs(tt$statistic), 3)
    
    ol_zin = ObsLikelihood(zin_likelihood, zin_control(b = .1, pi = .8))
    oo = ol_zin(eta = rep(100, 1000), sigma = 1)
    tt = t.test(oo, mu = 100)
    expect_lt(abs(tt$statistic), 3)
    
    oo = ol_zin(eta = eta_seq, sigma = 1)
    ll = lm(oo ~ eta_seq)
    ll0 = lm(oo ~ offset(eta_seq))
    expect_lt(anova(ll, ll0)$F[2], 7) # pf(7, 1, 999) > .991
})

test_that('Can sample NB', {
    ol_nb = ObsLikelihood(nb_likelihood)  
    oo = ol_nb(eta_seq, sigma = .5)
    ll = lm(oo ~ eta_seq)
    ll0 = lm(oo ~ offset(eta_seq))
    expect_lt(anova(ll, ll0)$F[2], 7) # pf(7, 1, 999) > .991
})