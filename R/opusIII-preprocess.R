num_person = 504 ## max 504: do not change
num_item =35 ## max 35

setwd("~/Dropbox/research/lsjm-art/lsjm-code")

library(dplyr)

source("R/art-functions.R")
source("R/opusIII-preprocess.R")

opusIII = readr::read_delim("data/opusIII-matrices-data.dat"," ")
df = as_tibble(opusIII[,-1])[,1:4] ## drop row names
name_item = unique(df$item)
pick_item = df$item %in% name_item[1:num_item]

df = df[pick_item, ]

pick_person = df$person %in% unique(df$person)[1:num_person]
df = df[pick_person,]

di = df[,-4]
dt = df[,-3]

nitem = length(unique(df$item))
nperson = length(unique(df$person))

time = pull(dt,RT)
ncut = 5
## interval <- seq(from=0, to = max(time)+1,length.out = 8)
pseq =  seq(from=0, to = 1, length.out = ncut + 1)
sj = quantile(time, probs = pseq) %>% round()
sj[1] = 0; sj[length(sj)] = sj[length(sj)] + 1

library(survival)
status = rep(1, nrow(dt))

tdf = data.frame(item = dt$item, person = dt$person, time = dt$RT, response = di$resp, status = 1)
tmp <- survSplit(formula = Surv(time, status) ~ ., data = tdf, cut = sj, episode ="seg_g") %>%
  mutate(seg = factor(tstart),
         seg_g = seg_g - 2,
         len = time - tstart
         ) %>% filter(status == 1) %>%
  as_tibble

item = pull(tmp, item)
person = pull(tmp, person)
seg_g = pull(tmp,seg_g)
H = pull(tmp,len)

mi = reshape2::dcast(tmp %>% select(item, person, response), item ~ person)[,-1]
mt = reshape2::dcast(tmp %>% select(item, person, time), item ~ person)[,-1]
mNA = mi; mNA[!is.na(mNA)] = 1;  mNA[is.na(mNA)] = 0
mseg = reshape2::dcast(tmp %>% select(item, person, seg_g), item ~ person)[,-1]
mh = reshape2::dcast(tmp %>% select(item, person, len), item ~ person)[,-1]
mlen = sj[2:(ncut+1)] - sj[1:(ncut)]

mi[is.na(mi)] = -99
mt[is.na(mt)] = -99
mseg[is.na(mseg)] = -99
mh[is.na(mh)] = -99

## data and fixed parameters
I = nrow(mt)
N = ncol(mt)
C = 2
G = ncut #
readr::write_csv(data.frame(I=I, N=N, C=C, G=G), "mvar.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mlen),"mlen.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mseg),"mseg.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mh),"mh.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mt),"mt.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mi),"mi.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mNA),"mNA.csv", col_names = FALSE)

mvar = readr::read_csv("mvar.csv", col_names=FALSE) %>% as.matrix()
I = mvar[1,1]; N = mvar[1,2]; C = mvar[1,3]; G = mvar[1,4];

## lambda
a_lambda = matrix(0.001,I,G)
b_lambda = matrix(0.001,I,G)
jump_lambda = matrix(1.0,I,G)

mu_beta = matrix(0.0,I,2)
sigma_beta = matrix(1.0,I,2)
jump_beta = matrix(0.25,I,2)

mu_theta = matrix(0.0,N,2)
sigma_theta = matrix(1.0,N,2)
jump_theta = matrix(1.0,N,2)

a_sigma = 1.0
b_sigma = 1.0

mu_gamma = matrix(0.0,1,2)
sigma_gamma = matrix(1.0,1,2)
jump_gamma = matrix(1.0,1,2)

mu_z = matrix(0.0,N,2)
sigma_z = matrix(1.0,N,2)
jump_z = matrix(1.0,N,2)

mu_w = matrix(0.0,I,2)
sigma_w = matrix(1.0,I,2)
jump_w = matrix(0.1,I,2)

readr::write_csv(as.data.frame(rbind(a_lambda,b_lambda,jump_lambda)),"pj_lambda.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_beta,sigma_beta,jump_beta)),"pj_beta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_theta,sigma_theta,jump_theta)),"pj_theta.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(a_sigma,b_sigma)),"pj_sigma.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_gamma,sigma_gamma,jump_gamma)),"pj_gamma.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_z,sigma_z,jump_z)),"pj_z.csv", col_names = FALSE)
readr::write_csv(as.data.frame(rbind(mu_w,sigma_w,jump_w)),"pj_w.csv", col_names = FALSE)
