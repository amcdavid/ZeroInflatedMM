N = 6
Nc = 5
Ntot = N*Nc
X = cbind(matrix(1, nrow = Ntot), treat=rep(c(0, 1), each = Nc, times = N/2))
beta = c(1,1)
clusters = gl(N, Nc)
# Using indicator matrices
#Ji <- t(as(clusters, Class = "sparseMatrix"))
#Zi <- t(KhatriRao(t(Ji), t(Xi)))

Ji  = split(seq_along(clusters), clusters)
Zlist = lapply(Ji, . %>% matrix(data = 1, ncol = 1, nrow = length(.)))
r_dm = RandomDesignMatrix(Z_list = Zlist, cluster_idx = Ji)

nr = RanefLikelihood(normal_ranef, Sigma = 1, normal_ranef_control())
