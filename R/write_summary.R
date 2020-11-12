ss = summary(mylist)
mm = ss$statistics[,"Mean"]
rr = c(grep("^beta",  names(mm)), grep("^theta",  names(mm)))
dout = data.frame(vname = names(mm[rr]), mean = mm[rr])
readr::write_csv(dout, paste0(out_dir, "beta_theta_mean.csv"))
