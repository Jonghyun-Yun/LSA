## out_dir <- "chessB_pn/"
## out_dir <- "chessB_np/"
out_dir <- "chessB_no_latent/"

## out_dir <- "chessB_pn_ncut2_zero_beta/"
out_dir <- "chessB_no_latent_ncut2_zero_beta/"
## out_dir <- "chessB_np_ncut2_zero_beta/"
## out_dir <- "chessB_no_latent_ncut2_zero_beta/"
out_dir <- "chessB_swdz_pp_ncut5_zero_beta/"
## out_dir <- "chessB_double_pp_ncut2_zero_beta/"

out_dir <- "chessB_double_nn_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_double_no_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_double_pp_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_nn_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_np_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_pn_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_pp_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_swdz_nn_ncut5_zero_beta_noinfo_lc2/"
out_dir <- "chessB_swdz_pp_ncut5_zero_beta_noinfo_lc2/"

num_chain <- 2
double_w <- 0
double_z <- 1
HAS_REF <- 0

## library(art)
library(coda)
library(dplyr)
library(stringr)
library(magrittr)
library(bayesplot)
library(foreach)
library(doParallel)
## library(timereg)
## registerDoParallel(cores = detectCores() - 1)
stopImplicitCluster()
registerDoParallel(2)

## setwd("/Users/yunj/Dropbox/research/lsjm-art/lsjm-code")

source("R/art-functions.R")
source("R/load-outputs.R")

## see get_hazitem() for "out"
source("R/chessB_noinfo-preprocess.R")
mH <- reshape2::melt(t(mh))
chain = 1
ii = 1
out <- gethaz_item(mylist[[chain]], cname, ii, N, theta = NULL, double_w, double_z)

seg_i = unlist(mseg[ii,])
H_i = mh[ii,]

resp_i = mi[ii,]
time_i = mt[ii,]

num_iter = nrow(mylist[[chain]])
  loglike = numeric(num_iter)

get_like <- function(loglike, out, resp_i, time_i, sj, H_i, seg_i) {

  rr = out$rr
  log_rr = out$log_rr
  lambda = out$lambda
  
  slam = hlam = list()
  hlam[[cc]] = t( t(lambda[[cc]]) * mlen)
  slam[[cc]] = rbind(0, apply(hlam[[cc]],  1,  cumsum))

  for (nn in 1:nim_iter) {
    loglike[nn] = loglike[nn] - sum((slam[[cc]][(seg_i+1), nn] + H_i * out$lambda[[cc]][nn, (seg_i+1)]) * rr[[cc]][nn, ])
    
  }
  clam[[cc]] = 
}

set.seed(1)
num_iter <- nrow(mylist[[1]])
sim_data <- list()
for (item in 1:I) {
  start_time <- proc.time()
  stk_tt <- foreach(chain = 1:num_chain, .combine = "rbind") %dopar% {
    out <- gethaz_item(mylist[[chain]], cname, item, N, theta = NULL, double_w, double_z)
    time <- foreach(nn = 1:num_iter, .combine = "rbind") %do% {
      gen_surv_time(out, sj, H, N, nn)
    }
    pp <- gen_surv_pp(out, time, sj) %>% array(dim = dim(time))
    cbind(time, pp)
  }
  sim_data[[item]] <- stk_tt
  cat("\nelapsed time to simulate item", item, "\n")
  print(proc.time() - start_time)
}

save(sim_data, file = paste0(out_dir, "surv_sim.RData"))
