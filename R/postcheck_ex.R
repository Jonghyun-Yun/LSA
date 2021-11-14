out_dir = "chessB_pn_ncut5_zero_beta_noinfo_lc2/"
out_dir = "duolingo_pn_ncut5_zero_beta_noinfo_lc2/"
num_chain = 1; HAS_REF = 0

## simulating samples
source("R/Renviron.R")
source("R/load-outputs.R")
source("R/wrap_param.R")

sim_data <- list()
len = diff(sj)
num_chain = 1
num_iter = min(unlist(lapply(mydf,  nrow)))
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

## visualizing performance metrics
dll = t(mll)
dauc = t(mauc)
colnames(mauc) = colnames(mll) = row.names(mi) = 1:I
## row.names(dauc) = row.names(dll) = colnames(mi) = 1:num_iter
dauc = reshape2::melt(mauc)
dll <- reshape2::melt(mll)
## mi_long <- reshape2::melt(t(mi))
dd = plyr::join(dll, dauc, by = c("Var1","Var2"))
names(dd) = c("iter","item","LogLoss","AUC")

ll_boxp <- ggplot(dd, aes(x=item,y=LogLoss,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 1)

auc_boxp <- ggplot(dd, aes(x=item,y=AUC,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 1)

pdf(paste0(out_dir,"figure/sim_cmetrics.pdf"))
print(ll_boxp)
print(auc_boxp)
dev.off()
