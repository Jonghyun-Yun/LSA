## out_dir = "chessB_pn_ncut5_zero_beta_noinfo_lc2/"
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

mytheme = theme(
  legend.position = "none",
  text = element_text(size = 35),
  ## plot.title =  element_text(size = 20, hjust = 0.5),
  ## axis.title = element_text(size = 15),
  axis.text = element_text(size = 25)
  ## panel.grid.minor = element_blank(),
  ## panel.background = element_blank()
  ##panel.border = element_blank()
  ##panel.grid.major = element_blank()
)

ll_boxp <- ggplot(dd, aes(x=item,y=LogLoss,fill=factor(item))) +
  geom_boxplot() +
  ## ylim(0, 1) +
  mytheme +
  scale_x_continuous(breaks = c(1,max(dd$item))) + scale_y_continuous(limits=c(0,1), breaks = c(0,0.5,1)) +
  labs(x = "Item", y = "Log-loss")

auc_boxp <- ggplot(dd, aes(x=item,y=AUC,fill=factor(item))) +
  geom_boxplot() +
  ## ylim(0, 1) +
  mytheme +
  scale_x_continuous(breaks = c(1,max(dd$item))) +  scale_y_continuous(limits=c(0,1), breaks = c(0,0.5,1)) +
  labs(x = "Item", y = "AUC")

pdf(paste0(out_dir,"figure/sim_cmetrics.pdf"))
print(ll_boxp)
print(auc_boxp)
dev.off()
