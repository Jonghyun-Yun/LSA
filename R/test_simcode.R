## out_dir <- "chessB_pn/"
## out_dir <- "chessB_no_latent/"

mt = readr::read_csv(file="chessB_pn/input/mt.csv", col_names=F) %>% as.matrix()
mi = readr::read_csv(file="chessB_pn/input/mi.csv", col_names=F) %>% as.matrix()

num_chain <- 3
double_z <- 0
double_w <- 0
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
registerDoParallel(2)
## stopImplicitCluster()

## setwd("/Users/yunj/Dropbox/research/lsjm-art/lsjm-code")

source("R/art-functions.R")
source("R/load-outputs.R")

save(mylist, file = paste0(out_dir,"mylist"))

set.seed(1)
num_iter <- nrow(mylist[[1]])
sim_data <- list()
mll  <- 99 * mi

for (item in 1:I) {
  start_time <- proc.time()
  stk_tt <- foreach(chain = 1:1, .combine = "rbind") %do% {
    out <- gethaz_item(mylist[[chain]], cname, item, N)
    time <- matrix(mt[item, ], nrow = num_iter, ncol = N, byrow=T)
    pp <- gen_surv_pp(out, time, sj) %>% array(dim = dim(time))
    ## cbind(time, pp)

  for (k in 1:N) {
    mll[item, k] = MLmetrics::LogLoss(pp[, k], mi[item, k])
    }

  }
  sim_data[[item]] <- stk_tt
  cat("\nelapsed time to simulate item", item, "\n")
  print(proc.time() - start_time)
}

dll = t(mll)
colnames(dll) = 1:I
row.names(dll) = 1:N
d_long <- reshape2::melt(dll)
names(d_long) = c("person","item","LogLoss")
ll_boxp <- ggplot(d_long, aes(x=item,y=LogLoss,fill=factor(item)))+
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 10.5)

pdf(paste0(out_dir,"LogLoss.pdf"))
print(ll_boxp)
dev.off()
