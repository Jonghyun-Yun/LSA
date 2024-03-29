num_person = 504 ## max 504: do not change
num_item =35 ## max 35

setwd("~/Dropbox/research/lsjm-art/lsjm-code")

library(dplyr)
library(magrittr)

source("R/art-functions.R")

opusIII = readr::read_delim("data/opusIII-matrices-data.dat"," ")
df = as_tibble(opusIII[,-1])[,1:4] ## drop row names
name_item = unique(df$item)
name_person = unique(df$person)
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

source("R/write-data.R")

mtab_sj = t( apply(mseg, 1, function(x) tab_sj(x,G)) )

tmp_0 = mseg; tmp_0[mi==1] = -99;
tmp_1 = mseg; tmp_1[mi==0] = -99;

mIY = rbind( t( apply(tmp_0, 1, function(x) tab_IY(x,G)) ), t( apply(tmp_1, 1, function(x) tab_IY(x,G)) ))

readr::write_csv(as.data.frame(mtab_sj),"input/mtab_sj.csv", col_names = FALSE)
readr::write_csv(as.data.frame(mIY),"input/mIY.csv", col_names = FALSE)

mvar = readr::read_csv("input/mvar.csv", col_names=FALSE) %>% as.matrix()
I = mvar[1,1]; N = mvar[1,2]; C = mvar[1,3]; G = mvar[1,4];

## lambda
a_lambda = matrix(0.1,I,G)
b_lambda = matrix(0.1,I,G)
jump_lambda = matrix(0.5,I,G)

mu_beta = matrix(0.0,I,2)
sigma_beta = matrix(sqrt(1.0),I,2)
jump_beta = matrix(0.25,I,2)

mu_theta = matrix(0.0,N,2)
sigma_theta = matrix(sqrt(1.0),N,2)
jump_theta = matrix(1.0,N,2)

a_sigma = 1.0
b_sigma = 1.0

mu_gamma = matrix(0.0,1,2)
sigma_gamma = matrix(sqrt(1.0),1,2)
jump_gamma = matrix(1.0,1,2)

mu_z = matrix(0.0,N,2)
sigma_z = matrix(sqrt(1.0),N,2)
jump_z = matrix(1.0,N,2)

mu_w = matrix(0.0,I,2)
sigma_w = matrix(sqrt(1.0),I,2)
jump_w = matrix(0.5,I,2)

source("R/write-prior.R")
