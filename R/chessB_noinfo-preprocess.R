source("R/chessB-preprocess.R")

mvar = readr::read_csv("input/mvar.csv", col_names=FALSE) %>% as.matrix()
I = mvar[1,1]; N = mvar[1,2]; C = mvar[1,3]; G = mvar[1,4];

## lambda
a_lambda = matrix(0.01,I,G)
b_lambda = matrix(0.01,I,G)
jump_lambda = matrix(0.5,I,G)

mu_beta = matrix(0.0,I,2)
sigma_beta = matrix(sqrt(100.0),I,2)
jump_beta = matrix(0.25,I,2)

mu_theta = matrix(0.0,N,2)
sigma_theta = matrix(sqrt(100.0),N,2)
jump_theta = matrix(1.0,N,2)

a_sigma = 0.01
b_sigma = 0.01

mu_gamma = matrix(0.0,1,2)
sigma_gamma = matrix(sqrt(100.0),1,2)
jump_gamma = matrix(0.01,1,2)

mu_z = matrix(0.0,N,2)
sigma_z = matrix(sqrt(100.0),N,2)
jump_z = matrix(0.5,N,2)

mu_w = matrix(0.0,I,2)
sigma_w = matrix(sqrt(100.0),I,2)
jump_w = matrix(0.25,I,2)

readr::write_csv(as.data.frame(rbind(a_lambda,b_lambda,jump_lambda)),"input/pj_lambda.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_beta,sigma_beta,jump_beta)),"input/pj_beta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_theta,sigma_theta,jump_theta)),"input/pj_theta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(a_sigma,b_sigma)),"input/pj_sigma.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_gamma,sigma_gamma,jump_gamma)),"input/pj_gamma.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_z,sigma_z,jump_z)),"input/pj_z.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_w,sigma_w,jump_w)),"input/pj_w.csv", col_names = FALSE)
