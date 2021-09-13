out_dir <- "duolingo_pn_ncut5_zero_beta_noinfo_lc2/"
num_chain <- 2
HAS_REF <- 0

source("R/Renviron.R")
source("R/load-outputs.R")

z0 <- matched$z0
w0 <- matched$w0

readr::write_csv(as.data.frame(w0), paste0(out_dir, "w0.csv"))
readr::write_csv(as.data.frame(z0), paste0(out_dir, "z0.csv"))
