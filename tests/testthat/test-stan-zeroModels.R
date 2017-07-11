Y <- matrix(1:8, nrow=8, ncol=7)
Y <- t(Y)
Y[!upper.tri(Y)] <- 0
obj <- MAST::FromMatrix(Y, cData = data.frame(X = 1, g=gl(4, 2)))
model <- ~1|g

test_that('Can fit stylized model', {
    test <- stan_frs(obj, model, modelargs = list(debug=0, marginal=1, ranefs=1), stanargs = list(iter=500, chains=1))
})
