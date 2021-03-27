source("R/Renviron.R")

mt = readr::read_csv(file="chessB_pn/input/mt.csv", col_names=F) %>% as.matrix()
mi = readr::read_csv(file="chessB_pn/input/mi.csv", col_names=F) %>% as.matrix()

source("R/load-outputs.R")

## save(mylist, file = paste0(out_dir,"mylist"))

source("R/wrap_param.R")

set.seed(1)
num_iter <- nrow(mylist[[1]])
sim_data <- list()
mll  <- mauc <- matrix(0, num_iter, I)

for (ii in 1:I) {
  start_time <- proc.time()
  ## stk_tt <- foreach(chain = 1:1, .combine = "rbind") %do% {
    ## out <- gethaz_item(mylist[[chain]], cname, ii, N, double_w, double_z)
    ## time <- matrix(mt[ii, ], nrow = num_iter, ncol = N, byrow=T)
    ## pp <- gen_surv_pp(out, time, sj) %>% array(dim = dim(time))
    pp <- art::rcpp_gen_surv_pp(lambda, theta, z, w, gamma, param, ii)
    ## cbind(time, pp)

  for (nn in 1:num_iter) {
    mll[nn, ii] = MLmetrics::LogLoss(pp[nn, ], mi[ii, ])
    mauc[nn, ii] = MLmetrics::AUC(pp[nn, ], mi[ii, ])
    }

  ## }
  ## sim_data[[item]] <- stk_tt
  cat("\nelapsed time to simulate item", ii, "\n")
  print(proc.time() - start_time)
}

mlike = art::get_loglike(lambda, theta, z, w, gamma, param)
mmet = list(mll=mll,mauc=mauc,mlike=mlike)
## save(mmet, file = paste0(out_dir,"mll_auc.RData"))

dll = t(mll)
dauc = t(mauc)
colnames(mauc) = colnames(mll) = row.names(mi) = 1:I
## row.names(dauc) = row.names(dll) = colnames(mi) = 1:num_iter
dauc = reshape2::melt(mauc)
dll <- reshape2::melt(mll)
## mi_long <- reshape2::melt(t(mi))
dd = plyr::join(dll, dauc, by = c("Var1","Var2"))
names(dd) = c("iter","item","LogLoss","AUC")

save(dd, file = paste0(out_dir,"dll_auc.RData"))

ll_boxp <- ggplot(dd, aes(x=item,y=LogLoss,fill=factor(item))) + 
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 1)

auc_boxp <- ggplot(dd, aes(x=item,y=AUC,fill=factor(item))) + 
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 1)

llike_boxp <- ggplot(data.frame(LogLike = mlike), aes(y=LogLike)) + geom_boxplot() + theme(legend.position = "none") + ggtitle(concat_summary(mlike,0)) + theme(plot.title = element_text(size=10))

pdf(paste0(out_dir,"figure/cmetrics.pdf"))
print(ll_boxp)
print(auc_boxp)
print(llike_boxp)
dev.off()
