pdf("figure/lambda_mcmc_interval_plot.pdf")
p0 = mcmc_intervals(
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

p1 = mcmc_intervals(
  mylist,
  regex_pars = "^lambda\\.1\\.1\\.",
  transformations = "log"
)
print(p0)
print(p1)
dev.off(which = dev.cur())

pdf("figure/z_pairs_plot.pdf")
p = mcmc_pairs(mylist,
               regex_pars = "^z.[0-1]\\.1\\.",
               off_diag_args = list(size = 0.75))
print(p)
dev.off(which = dev.cur())
