Y <- matrix(1:8, nrow=8, ncol=7)
Y <- t(Y)
Y[!upper.tri(Y)] <- 0
obj <- MAST::FromMatrix(Y, cData = data.frame(g=gl(4, 2)))
model <- ~1|g
test <- stan_frs(obj, model)
