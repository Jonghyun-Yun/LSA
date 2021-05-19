out_dir = "chessB_pn_ncut5_zero_beta_noinfo_lc2/"
num_chain = 1; HAS_REF = 0
source("R/Renviron.R")
source("R/load-outputs.R")
source("R/wrap_param.R")

sim_data <- list()
len = diff(sj)
num_chain = 1
mll  <- mauc <- matrix(0, num_iter, I)

for(ii in 1:I){
start_time <- proc.time()
  stk_tt <- foreach(chain = 1:num_chain, .combine = "rbind") %do% {
    simtp <- art::rcpp_gen_surv(lambda, theta, z, w, gamma, param, ii)
    cbind(simtp$time, simtp$pp)
    }
  sim_data[[ii]] <- stk_tt
  pp = simtp$pp
  for (nn in 1:num_iter) {
    mll[nn, ii] = MLmetrics::LogLoss(pp[nn, ], mi[ii, ])
    mauc[nn, ii] = MLmetrics::AUC(pp[nn, ], mi[ii, ])
    }
  cat("\nelapsed time to simulate item", ii, "\n")
  print(proc.time() - start_time)
}
