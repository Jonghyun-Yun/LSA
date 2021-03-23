## out_dir <- "chessB_pn/"
## out_dir <- "chessB_np/"
out_dir <- "chessB_no_latent/"

## out_dir <- "chessB_pn_ncut2_zero_beta/"
out_dir <- "chessB_no_latent_ncut2_zero_beta/"
## out_dir <- "chessB_np_ncut2_zero_beta/"
## out_dir <- "chessB_no_latent_ncut2_zero_beta/"
out_dir <- "chessB_swdz_pp_ncut2_zero_beta/"
out_dir <- "chessB_double_pp_ncut2_zero_beta/"

num_chain <- 2
double_w <- 1
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
