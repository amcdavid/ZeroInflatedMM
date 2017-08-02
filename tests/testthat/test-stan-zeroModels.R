Y <- matrix(1:8, nrow=8, ncol=7)
Y <- t(Y)
Y[!upper.tri(Y)] <- 0
cData = data.frame(X = 1, g=gl(4, 2))
obj <- MAST::FromMatrix(Y, cData = cData)
model <- ~1|g

test_that('Can compile model', {
    rstan::stanc(stan_modelfile)  
})

test_that('Can fit stylized model', {
    test <- fit_FxCHM(obj, model, model_control = fxc_model_control(debug = 0, marginal = 1, ranefs = 0, center = FALSE), stan_control = fxc_stan_control(iter=50, chains=1), contrasts = '(Intercept)', returnfit = TRUE)
})


test_coefs = function(stan, expected, par, idx){
    mean_se = summary(fit(stan), par = par)$summary[idx,]
    Z = abs(mean_se[,'mean'] - expected)/mean_se[,'sd']
    expect_true(all(Z<2))
}

test_that('Converge with large data', {
    nr = 30
    Y = matrix(0, nrow = 3, ncol = 8)
    Y[1,] = 1
    Y[2,5:8] = 2
    Y[3,7:8] = 4
    cd_rep = do.call(rbind, rep(list(cData), nr))
    y_rep = do.call(cbind, rep(list(Y), nr))
    y_rep[y_rep!=0] = y_rep[y_rep!=0] + rnorm(sum(y_rep!=0))*.1
    y_rep[1,1] = 0 #to keep us off the parameter boundary
    obj = MAST::FromMatrix(y_rep, cData=cd_rep)
    obj = obj[,order(colData(obj)$g)]
    test_conditional <- fit_FxCHM(obj, model, model_control = fxc_model_control(debug = 0, marginal = 0, ranefs = 0, center = FALSE), stan_control = fxc_stan_control(iter=500, cores = 4), contrasts = '(Intercept)', returnfit = TRUE, method = 'hmc')
    test_coefs(test_conditional, c(1, 2, 4), 'beta_Tf_G', c(2, 4, 6))
    test_marginal <- fit_FxCHM(obj, model, 
                               model_control = fxc_model_control(debug = 0, marginal = 1, ranefs = 0, center = FALSE), stan_control = fxc_stan_control(iter=500, cores = 4), contrasts = '(Intercept)', returnfit = TRUE, method = 'hmc')
    test_coefs(test_marginal, c(1, 1, 1), 'beta_Tf_G', c(2, 4, 6))
})

