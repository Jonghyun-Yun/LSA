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
out_dir <- "chessB_double_pp_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_nn_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_pn_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_pp_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_swdz_nn_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_swdz_pp_ncut5_zero_beta_noinfo_lc2/"
## out_dir <- "chessB_no_latent_ncut5_zero_beta_noinfo_lc2/"
out_dir = "chessB_np_ncut5_zero_beta_noinfo_lc2/"
out_dir="chessB_np_ncut5_zero_beta_noinfo_lc2/"

num_chain <- 2
HAS_REF <- 0

source("R/Renviron.R")
source("R/load-outputs.R")

library(art)
ddpp = get_loglike(lambda, theta, z, w, gamma, param)

## see get_hazitem() for "out"
chain = 1
sam = mylist[[chain]][seq(1, nrow(mylist[[chain]]), 50), ]
num_iter = nrow(sam)

mloglike = foreach(ii = 1:I, .combine="cbind") %dopar% {
out <- gethaz_item(sam, cname, ii, N, theta = NULL, double_w, double_z)

seg_i = unlist(mseg[ii,])
H_i = mh[ii,]

resp_i = mi[ii,]
time_i = mt[ii,]

get_ll(out, resp_i, time_i, sj, H_i, seg_i, loglike = NULL)
}

loglike = rowSums(mloglike)

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
