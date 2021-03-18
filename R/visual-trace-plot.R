pdf("figure/lambda_mcmc_interval_plot.pdf")
p0 <- mcmc_intervals(
  mylist,
  regex_pars = "^lambda\\.0\\.1\\.",
  transformations = "log"
)
## mcmc_areas(
##  lambda0.sam,
##  prob = 0.8, # 80% intervals
##  prob_outer = 0.99, # 99%
##  point_est = "mean"
## )

p1 <- mcmc_intervals(
  mylist,
  regex_pars = "^lambda\\.1\\.1\\.",
  transformations = "log"
)
print(p0)
print(p1)
dev.off(which = dev.cur())

pdf("figure/beta_parcoord.pdf")
p <- bayesplot::mcmc_parcoord(mylist,
  regex_pars = "^beta\\.[0-9]\\."
)
print(p)
dev.off(which = dev.cur())

pdf("figure/theta_parcoord.pdf")
p <- bayesplot::mcmc_parcoord(mylist,
  regex_pars = "^theta\\.[0-9]\\."
)
print(p)
dev.off(which = dev.cur())

pdf("figure/lambda_parcoord.pdf")
p <- bayesplot::mcmc_parcoord(mylist,
  regex_pars = "^lambda\\.[0-1]\\.1\\.",
  transformations = "log"
)
print(p)
dev.off(which = dev.cur())

pdf("figure/w_parcoord_plot.pdf")
p <- mcmc_parcoord(mylist,
  regex_pars = "^w\\.[0-1]\\.[1-5]\\."
)
print(p)
dev.off(which = dev.cur())

## mcmc_intervals(mylist, pars=c("lambda.0.1.1")
pdf("figure/dist_mcmc_trace_plot.pdf")
p0 <- mcmc_trace(mylist,
  regex_pars = "^dist_z.[0-1]\\.[0-2]_w",
  ## transformations = "log",
  facet_args = list(nrow = 2, labeller = label_parsed)
)
print(p <- p + facet_text(size = 15))
p <- mcmc_trace(mylist,
  regex_pars = "^dist_z.[0-1]\\.[3-6]_w",
  ## transformations = "log",
  facet_args = list(nrow = 2, labeller = label_parsed)
)
print(p <- p + facet_text(size = 15))
p <- mcmc_trace(mylist,
  regex_pars = "^dist_z.[0-1]\\.[7-9]_w",
  ## transformations = "log",
  facet_args = list(nrow = 2, labeller = label_parsed)
)
print(p <- p + facet_text(size = 15))
dev.off(which = dev.cur())

colnames(df)[grepl("^w\\.1\\.", colnames(df))]


pdf("figure/lambda_mcmc_trace_plot.pdf")
color_scheme_set("mix-blue-pink")
for (item in 1:I) {
  p <- mcmc_trace(mylist,
    regex_pars = paste0("^lambda\\.0\\.", item, "\\."),
    transformations = "log",
    facet_args = list(nrow = 2, labeller = label_parsed)
  )
  print(p <- p + facet_text(size = 15) + ggtitle(paste0("lambda.1 trace for item ", item, ".")))

  p <- mcmc_trace(mylist,
    regex_pars = paste0("^lambda\\.1\\.", item, "\\."),
    transformations = "log",
    facet_args = list(nrow = 2, labeller = label_parsed)
  )
  print(p <- p + facet_text(size = 15) + ggtitle(paste0("lambda.1 trace for item ", item, ".")))
}
dev.off(which = dev.cur())

pdf("figure/beta_mcmc_trace_plot.pdf")
for (i in 1:I) {
  color_scheme_set("mix-blue-pink")
  p <- mcmc_trace(mylist,
    regex_pars = paste0("^beta\\.", i, "\\."),
    facet_args = list(nrow = 2, labeller = label_parsed)
  )
  p + facet_text(size = 15)
  print(p <- p + facet_text(size = 15))
}
dev.off(which = dev.cur())

pdf("figure/theta_mcmc_trace_plot.pdf")
color_scheme_set("mix-blue-pink")
for (k in 1:50) {
  p <- mcmc_trace(mylist,
    regex_pars = paste0("^theta\\.", k, "\\."),
    facet_args = list(nrow = 2, labeller = label_parsed)
  )
  print(p <- p + facet_text(size = 15))
}
dev.off(which = dev.cur())

pdf("figure/gamma_mcmc_trace_plot.pdf")
color_scheme_set("mix-blue-pink")
p <- mcmc_trace(mylist,
  regex_pars = "^gamma\\.",
  facet_args = list(nrow = 2, labeller = label_parsed)
)
print(p <- p + facet_text(size = 15))
dev.off(which = dev.cur())

pdf("figure/lp_sigma_mcmc_trace_plot.pdf")
color_scheme_set("mix-blue-pink")
p <- mcmc_trace(mylist,
  pars = c("sigma", "lp_"),
  facet_args = list(nrow = 2, labeller = label_parsed)
)
print(p <- p + facet_text(size = 15))
dev.off(which = dev.cur())
