Y <- matrix(1:8, nrow=8, ncol=7)
Y <- t(Y)
Y[!upper.tri(Y)] <- 0
obj <- MAST::FromMatrix(Y, cData = data.frame(X = 1, g=gl(4, 2)))
model <- ~1|g

test_that('Can compile model', {
    rstan::stanc(stan_modelfile)  
})

test_that('Can fit stylized model', {
    test <- fit_FxCHM(obj, model, model_control = fxc_model_control(debug = 0, marginal = 1, ranefs = 0, center = FALSE), stan_control = fxc_stan_control(iter=50, chains=1, control=list(max_treedepth = 12)), contrasts = '(Intercept)', returnfit = TRUE)
})

