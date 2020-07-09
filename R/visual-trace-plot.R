pdf("figure/lambda_mcmc_interval_plot.pdf")
mcmc_intervals(
  mylist,
  regex_pars = "^lambda\\.0\\.1\\.",
  transformations = "log"
)
##mcmc_areas(
##  lambda0.sam,
##  prob = 0.8, # 80% intervals
##  prob_outer = 0.99, # 99%
##  point_est = "mean"
##)

mcmc_intervals(
  mylist,
  regex_pars = "^lambda\\.1\\.1\\.",
  transformations = "log"
)
dev.off(which = dev.cur())

pdf("figure/z_pairs_plot.pdf")
mcmc_pairs(mylist,
           regex_pars = "^z.[0-1]\\.1\\.",
           off_diag_args = list(size = 0.75))
dev.off(which = dev.cur())
