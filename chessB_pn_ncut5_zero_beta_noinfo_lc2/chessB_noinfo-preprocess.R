source("R/chessB-preprocess.R")

mvar = readr::read_csv("input/mvar.csv", col_names=FALSE) %>% as.matrix()
I = mvar[1,1]; N = mvar[1,2]; C = mvar[1,3]; G = mvar[1,4];

## lambda
l_c = 2
l_m = ((G - 1) / (quantile(time, 0.75)  * (G - 1:G + 0.5)))
a_lambda = matrix(l_m/l_c, I, G, byrow=T)
b_lambda = matrix(1 / l_c,I,G)
## a_lambda = matrix(0.0001, I, G, byrow=T)
## b_lambda = matrix(0.0001, I, G)
jump_lambda = matrix(0.5,I,G)

mu_beta = matrix(0.0,I,2)
sigma_beta = matrix(sqrt(10000.0),I,2)
jump_beta = matrix(0.25,I,2)

mu_theta = matrix(0.0,N,2)
sigma_theta = matrix(sqrt(10000.0),N,2)
jump_theta = matrix(0.2,N,2)

a_sigma = 0.0001
b_sigma = 0.0001

mu_gamma = matrix(0.0,1,2)
sigma_gamma = matrix(sqrt(2.0),1,2)
jump_gamma = matrix(0.01,1,2)

mu_z = matrix(0.0,N,2)
sigma_z = matrix(sqrt(1.0),N,2)
jump_z = matrix(0.5,N,2)

mu_w = matrix(0.0,I,2)
sigma_w = matrix(sqrt(1.0),I,2)
jump_w = matrix(0.15,I,2)

source("R/write-prior.R")
