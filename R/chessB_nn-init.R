mvar = as.matrix( readr::read_csv("input/mvar.csv", col_names=F) )
I = mvar[1]
N = mvar[2]
G = mvar[4]

set.seed(as.numeric(Sys.time()))

init_lambda = matrix(1, nrow = 2*I, ncol = G);
init_beta = 0*matrix(rnorm(2*I), ncol = 2);
init_theta = 0*matrix(rnorm(2*N), ncol = 2);

init_gamma = 1*matrix(c(-1, -1), ncol = 1);

init_w = 0*matrix(rnorm(2*2*I), ncol = 2);
init_z = 0*matrix(rnorm(2*2*N), ncol = 2);

source("R/write-init.R")
