mvar = as.matrix( readr::read_csv("input/mvar.csv", col_names=F) )
I = mvar[1]
N = mvar[2]
G = mvar[4]

set.seed(as.numeric(Sys.time()))

init_lambda = matrix(1, nrow = 2*I, ncol = G);
init_beta = 0*matrix(rnorm(2*I), ncol = 2);
init_theta = 0*matrix(rnorm(2*N), ncol = 2);

init_gamma = 1*matrix(c(-1, 1), ncol = 1);

init_w = 1*matrix(rnorm(2*2*I), ncol = 2);
init_z = 1*matrix(rnorm(2*2*N), ncol = 2);

readr::write_csv(as.data.frame(init_lambda), "input/init_lambda.csv", col_names = FALSE)
readr::write_csv(as.data.frame(init_beta), "input/init_beta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(init_theta), "input/init_theta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(init_gamma), "input/init_gamma.csv", col_names = FALSE)
readr::write_csv(as.data.frame(init_w), "input/init_w.csv", col_names = FALSE)
readr::write_csv(as.data.frame(init_z), "input/init_z.csv", col_names = FALSE)
