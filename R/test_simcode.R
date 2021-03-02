out_dir <- "chessB_pn/"
## out_dir <- "chessB_no_latent/"
## out_dir <- "chessB_np/"

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

## save(mylist, file = paste0(out_dir,"mylist"))

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

save(mll,file = paste0(out_dir,"mll.RData"))

dll = t(mll)
colnames(dll) = row.names(mi) = 1:I
row.names(dll) = colnames(mi) = 1:N
d_long <- reshape2::melt(dll)
mi_long <- reshape2::melt(t(mi))
d_long = plyr::join(d_long, mi_long, by = c("Var1","Var2"))
names(d_long) = c("person","item","LogLoss","res")
d_long$res = factor(d_long$res, labels=c("incorrect","correct"))
ll_boxp <- ggplot(d_long, aes(x=item,y=LogLoss,fill=factor(item))) +  facet_wrap(~res) +
  geom_boxplot() + theme(legend.position = "none") + ylim(0, 10.5)

pdf(paste0(out_dir,"LogLoss_by.pdf"))
print(ll_boxp)
dev.off()

out_dir <- "chessB_pn/"
## out_dir <- "chessB_no_latent/"
## out_dir <- "chessB_np/"

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
stopImplicitCluster()
registerDoParallel(2)

## setwd("/Users/yunj/Dropbox/research/lsjm-art/lsjm-code")

source("R/art-functions.R")
## source("R/load-outputs.R")

## save(mylist, file = paste0(out_dir,"mylist"))

out_dir = "chessB_np/"
## out_dir = "chessB_pn/"
## out_dir = "chessB_no_latent/"


## log loss using simulated time
simll  <- 99 * mi

load(paste0(out_dir, "surv_sim.RData"))


for (item in 1:I) {
  pp = sim_data[[i]][, (N+1): (2*N)]
  for (k in 1:N) {
    simll[item, k] = MLmetrics::LogLoss(pp[, k], mi[item, k])
    }
  }

save(simll,file = paste0(out_dir,"simll.RData"))
