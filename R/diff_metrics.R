out_dir="chessB_double_no_ncut5_zero_beta_noinfo_lc2/"
## out_dir="chessB_nn_ncut2_zero_beta_noinfo_lc2/"

num_chain <- 2
double_z <- 0
double_w <- 0
HAS_REF <- 0

source("R/test_simcode.R")
ddd = dd

out_dir="chessB_double_no_ncut5_zero_beta_noinfo_lc2/"
## out_dir="chessB_double_pp_ncut5_zero_beta_noinfo_lc2/"
out_dir="chessB_double_nn_ncut5_zero_beta_noinfo_lc2/"

num_chain <- 2
double_z <- 1
double_w <- 1
HAS_REF <- 0

source("R/test_simcode.R")
ddd$dLogLoss = ddd$LogLoss - dd$LogLoss
ddd$dAUC = ddd$AUC - dd$AUC

## logloss (smaller, better)
ll_boxp <- ggplot(ddd, aes(x=item,y=dLogLoss,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(-1, 1)

## AUC (larger, better)
auc_boxp <- ggplot(ddd, aes(x=item,y=dAUC,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(-1, 1)

pdf(paste0(out_dir,"dmetrics.pdf"))
print(ll_boxp)
print(auc_boxp)
dev.off()

cat("Plot(s) are saved in",out_dir,"\n")
cat("=================================\n")
