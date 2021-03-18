mvar <- readr::read_csv(paste0(out_dir, "mvar.csv"), col_names = F) %>% as.matrix()
I <- mvar[1]
N <- mvar[2]
G <- mvar[4]

sj <- readr::read_csv(paste0(out_dir, "mlen.csv"), col_names = F) %>% as.matrix()
sj <- c(0, cumsum(sj))
H <- sj[2:(G + 1)] - sj[1:G]

cnames <- c(".chain", ".iteration")
for (c in 0:1) {
  for (i in 1:I) {
    for (g in 1:G) {
      cnames <- c(cnames, paste0("lambda.", c, ".", i, ".", g))
    }
  }
}

for (k in 1:N) {
  for (c in 0:1) {
    cnames <- c(cnames, paste0("theta.", k, ".", c))
  }
}
for (i in 1:I) {
  for (c in 0:1) {
    cnames <- c(cnames, paste0("beta.", i, ".", c))
  }
}

for (c in 0:1) {
  for (k in 1:N) {
    for (d in 1:2) {
      cnames <- c(cnames, paste0("z.", c, ".", k, ".", d))
    }
  }
}
for (c in 0:1) {
  for (i in 1:I) {
    for (d in 1:2) {
      cnames <- c(cnames, paste0("w.", c, ".", i, ".", d))
    }
  }
}

for (c in 0:1) {
  cnames <- c(cnames, paste0("gamma.", c))
}

cnames <- c(cnames, "sigma", "lp_")

## mythin = 10
## mystart = 5001
## myend = 25000

no_z1 <- !grepl("^z\\.1\\.", cnames)
no_w1 <- !grepl("^w\\.1\\.", cnames)

dlist <- list()

for (cid in 1:num_chain) {
  ## pisa-KR KR-sci data should skip 1000
  dlist[[cid]] <- readr::read_csv(paste0(out_dir, "sample_chain", cid, ".csv"), col_names = F, skip = 0) %>% as.data.frame()
  ## dlist[[cid]] = readr::read_csv(paste0(out_dir,"sample_chain",cid,".csv"), col_names=F) %>% as.data.frame()
  colnames(dlist[[cid]]) <- cnames
  if (!double_z && !double_w) {
    dlist[[cid]] <- dlist[[cid]][, no_z1 & no_w1] ## remove duplicates when single_z and single_w
  } else if (!double_w) {
    dlist[[cid]] <- dlist[[cid]][, no_w1] ## remove duplicates when single_w
  } else if (!double_z) {
    dlist[[cid]] <- dlist[[cid]][, no_z1] ## remove duplicates when single_z
  }
}

## mylist[[cid]] = mcmc(df, start = mystart, end = myend, thin = mythin)

if (!HAS_REF) {
  Xstar <- find_xstar_inlist(dlist)
  readr::write_csv(as.data.frame(Xstar$z.0), paste0(out_dir, "z0star.csv"), col_names = FALSE)
  readr::write_csv(as.data.frame(Xstar$w.0), paste0(out_dir, "w0star.csv"), col_names = FALSE)
}

if (HAS_REF) {
  Xstar <- list()
  Xstar$z.0 <- readr::read_csv(paste0(ref_dir, "z0star.csv"), col_names = FALSE) %>% as.matrix()
  Xstar$w.0 <- readr::read_csv(paste0(ref_dir, "w0star.csv"), col_names = FALSE) %>% as.matrix()
}

matched <- my_procrustes(Xstar, dlist, is_list = T)
## , my_translation = TRUE, my_scale = FALSE)
mydf <- matched$dlist
mdf <- bind_rows(matched$dlist, .id = "column_label")

mylist <- mcmc.list()
item <- 1
cname <- names(mydf[[1]])
mylist <- mcmc.list()
for (cid in 1:num_chain) {
  for (c in 0:1) {
    for (k in 1:N) {
      z <- mydf[[cid]][, stringr::str_which(cname, paste0("z\\.", c * double_z, "\\.", k, "\\.[1-2]"))]
      w <- mydf[[cid]][, stringr::str_which(cname, paste0("w\\.", c * double_w, "\\.", item, "\\."))]
      mydf[[cid]][[paste0("dist_z.", c, ".", k, "_", "w.", c, ".", item)]] <- sqrt(rowSums((z - w)^2))
    }
  }
  mylist[[cid]] <- mcmc(mydf[[cid]])
}
