out_dir="chessB_no_latent_ncut5_zero_beta_noinfo_lc2/"
## out_dir="chessB_nn_ncut2_zero_beta_noinfo_lc2/"

num_chain <- 1
HAS_REF <- 0
source("R/test_simcode.R")

## fname = paste0(out_dir,"dll_auc.RData")
## if (file.exists(fname)) load(file = fname) else source("R/test_simcode.R")

ddd = dd # baseline to ddd

lout = NULL
lout=c(lout,"chessB_double_pp_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_double_nn_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_double_np_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_double_pn_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_np_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_pn_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_pp_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_nn_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_swdz_nn_ncut5_zero_beta_noinfo_lc2/")
lout=c(lout,"chessB_swdz_pp_ncut5_zero_beta_noinfo_lc2/")

## fname = paste0(out_dir,"dll_auc.RData")
## if (file.exists(fname)) load(file = fname) else source("R/test_simcode.R")

num_chain <- 1
HAS_REF <- 0

for (out_dir in lout) {

source("R/test_simcode.R")
ddd$dLogLoss = ddd$LogLoss - dd$LogLoss
ddd$dAUC = ddd$AUC - dd$AUC

## logloss (smaller, better)
ll_boxp <- ggplot(ddd, aes(x=item,y=dLogLoss,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(-0.25, 0.25)

## AUC (larger, better)
auc_boxp <- ggplot(ddd, aes(x=item,y=dAUC,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ylim(-0.25, 0.25)

pdf(paste0(out_dir,"dmetrics.pdf"))
print(ll_boxp)
print(auc_boxp)
dev.off()

cat("Plot(s) are saved in",out_dir,"\n")
cat("=================================\n")

sink("loglikelihood.txt", append=TRUE)
cat(out_dir,"\n")
summary(mlike)
cat("=============================\n")

sink()

sink("DIC.txt", append=TRUE)
cat("out_dir, DIC1, DIC2\n")
cat(paste0(out_dir,","),paste0(round(mm$DIC1,2),","),round(mm$DIC2,2),"\n")

sink()
}
