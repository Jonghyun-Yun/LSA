setwd("~/research/lsjm-art/lsjm-code")

library(dplyr)
library(magrittr)

source("R/art-functions.R")

library(foreign)
chess <- read.spss("data/opennkweb.sav", to.data.frame = TRUE)
dataset.labels <- as.data.frame(attr(chess, "variable.labels"))
PPNR = chess[,1]

dt = cbind(PPNR, chess[,grepl("BR[1-9]+", names(chess))]) %>% na.omit
di = cbind(PPNR, chess[,grepl("B[1-9]+", names(chess))]) %>% filter(PPNR %in% dt$PPNR)

delo = cbind(PPNR, chess['ELO']) %>% filter(PPNR %in% dt$PPNR)

nitem = 40
nperson = nrow(di)

di_long <- reshape2::melt(di, id.vars=c("PPNR"))
dt_long <- reshape2::melt(dt, id.vars=c("PPNR"))

identical(di_long[,1],dt_long[,1])

dit = cbind(di_long, dt_long[,3])
colnames(dit) = c("person","item","resp","RT")

tab_item = data.frame(chr = colnames(di)[-1], num = 1:nitem)

dit$item = to_numID(dit$item, tab_item)

time = pull(dit, RT)
item = pull(dit, item)
resp = pull(dit, resp)

ncut = 5
## interval <- seq(from=0, to = max(time)+1,length.out = 8)
pseq =  seq(from=0, to = 1, length.out = ncut + 1)
sj = quantile(time, probs = pseq) %>% round()
sj[1] = 0; sj[length(sj)] = sj[length(sj)] + 1

library(survival)
status = rep(1, nrow(dt))

tdf = data.frame(item = dit$item, person = dit$person, time = dit$RT, response = dit$resp, status = 1)
tmp <- survSplit(formula = Surv(time, status) ~ ., data = tdf, cut = sj, episode ="seg_g") %>%
  mutate(seg = factor(tstart),
         seg_g = seg_g - 2,
         len = time - tstart
         ) %>% filter(status == 1) %>%
  as_tibble
