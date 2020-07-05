setwd("~/Dropbox/research/lsjm-art/lsjm-code")

source("Rfunction.R")

library(tidyverse)
load("data/pisa2015/US_PISA2015.rdata")
df = as_tibble(US_2015)
info = readr::read_csv("data/pisa2015/ItemTimeInfo.csv")
polytomous = c("DS519Q01C","DS498Q04C","DS465Q01C","CS635Q01S", "CS635Q04S","DS635Q05C","DS605Q04C","DS607Q03C","CS634Q02S", "CS645Q01S","DS657Q04C","DS629Q01C","CS637Q02S")
pdx = which(colnames(df) %in% polytomous)
pdx = c(pdx,pdx + 184)
df[,-pdx]

di = df[,1:187]
dt = df[,c(1:3,188:371)]
its = pullit(info,"S01")
dt01 = dt %>% select(CNTSCHID,CNTSTUID,ST004D01T,any_of(its[,2])) %>% na.omit
di01 = di %>% select(CNTSCHID,CNTSTUID,ST004D01T,any_of(its[,1])) %>% filter(CNTSTUID %in% dt01$CNTSTUID)

di01[di01 == 2] = 1
di01[is.na(di01)] = 999

di01 = di01[pick_person,pick_item]
dt01 = dt01[pick_person,pick_item]

tab_schid = tabulate_id(di01$CNTSCHID)
tab_stuid = tabulate_id(di01$CNTSTUID)
tab_item = tabulate_id(colnames(di01)[-(1:3)])

di01$schid = to_numID(dt01$CNTSCHID, tab_schid)
di01$stuid = to_numID(dt01$CNTSTUID, tab_stuid)
dt01$schid = to_numID(dt01$CNTSCHID, tab_schid)
dt01$stuid = to_numID(dt01$CNTSTUID, tab_stuid)

##di01_long <- melt(di01, id.vars=c("CNTSCHID","CNTSTUID","ST004D01T"))
##dt01_long <- melt(dt01, id.vars=c("CNTSCHID","CNTSTUID","ST004D01T"))

di01_long = di01 %>% select(- CNTSCHID, - CNTSTUID)
dt01_long = dt01 %>% select(- CNTSCHID, - CNTSTUID)

di01_long <- reshape2::melt(di01_long, id.vars=c("schid","stuid","ST004D01T"))
dt01_long <- reshape2::melt(dt01_long, id.vars=c("schid","stuid","ST004D01T"))

identical(di01_long[,1],dt01_long[,1])
identical(di01_long[,2],dt01_long[,2])
identical(di01_long[,3],dt01_long[,3])

dit01 = cbind(di01_long, dt01_long[,5])
colnames(dit01)[4:6] = c("item","res","time")

dit01$item = to_numID(dit01$item, tab_item)

sdi01 = di01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID) %>% as.matrix()
sdt01 = dt01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID) %>% as.matrix()

time = pull(dit01,time)
## ncut = 5

## interval <- seq(from=0, to = max(time)+1,length.out = 8)
pseq =  seq(from=0, to = 1, length.out = ncut + 1)
sj = quantile(time, probs = pseq) %>% round()
sj[1] = 0; sj[length(sj)] = sj[length(sj)] + 1
N = ncol(sdi01)
msj = array(0,dim=c(N,ncut+1,2))
for (i in 1:N) {
  msj[i,1,1]  = msj[i,1,1] = 0
  msj[i,2:(ncut+1),2] = quantile(sdt01[sdi01[,i]==1,i], probs = pseq[-1]) %>% round()
  msj[i,2:(ncut+1),1] = quantile(sdt01[sdi01[,i]==0,i], probs = pseq[-1]) %>% round()
}

library(survival)
tt01 = dt01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID)
ti01 = di01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID)
nitem = ncol(tt01)

tf01 = data.frame(time = c(as.matrix(tt01)), status = 1)
tmp <- survSplit(formula = Surv(time, status) ~ ., data = tf01, cut = sj, episode ="seg_g") %>%
    mutate(seg = factor(tstart),
           seg_g = seg_g - 2,
           len = time - tstart
           ) %>% filter(status == 1) %>%
    as_tibble

mseg = matrix(pull(tmp,seg_g),ncol=nitem) %>% t()
mh = matrix(pull(tmp,len),ncol=nitem) %>% t()
mlen = sj[2:(ncut+1)] - msj[1:(ncut)]
mt = dt01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID) %>% as.matrix %>% t()
mi = di01 %>% select(-ST004D01T,-schid,-stuid) %>% select(- CNTSCHID, - CNTSTUID) %>% as.matrix %>% t()

## data and fixed parameters
I = nrow(mt)
N = ncol(mt)
C = 2
G = ncut # 
write_csv(data.frame(I=I, N=N, C=C, G=G), "mvar.csv", col_names = FALSE)
write_csv(as.data.frame(mlen),"mlen.csv", col_names = FALSE)
write_csv(as.data.frame(mseg),"mseg.csv", col_names = FALSE)
write_csv(as.data.frame(mh),"mh.csv", col_names = FALSE)
write_csv(as.data.frame(mt),"mt.csv", col_names = FALSE)
write_csv(as.data.frame(mi),"mi.csv", col_names = FALSE)

mvar = readr::read_csv("mvar.csv", col_names=FALSE) %>% as.matrix()
I = mvar[1,1]; N = mvar[1,2]; C = mvar[1,3]; G = mvar[1,4];

## lambda
a_lambda = matrix(0.001,I,G)
b_lambda = matrix(0.001,I,G)
jump_lambda = matrix(1.0,I,G)

mu_beta = matrix(0.0,I,2)
sigma_beta = matrix(1.0,I,2)
jump_beta = matrix(0.5,I,2)

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
jump_w = matrix(0.5,I,2)

write_csv(as.data.frame(rbind(a_lambda,b_lambda,jump_lambda)),"pj_lambda.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(mu_beta,sigma_beta,jump_beta)),"pj_beta.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(mu_theta,sigma_theta,jump_theta)),"pj_theta.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(a_sigma,b_sigma)),"pj_sigma.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(mu_gamma,sigma_gamma,jump_gamma)),"pj_gamma.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(mu_z,sigma_z,jump_z)),"pj_z.csv", col_names = FALSE)
write_csv(as.data.frame(rbind(mu_w,sigma_w,jump_w)),"pj_w.csv", col_names = FALSE)
