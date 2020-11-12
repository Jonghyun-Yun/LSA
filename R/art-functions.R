pullit = function(info,cl) {
  it = info %>% filter(Cluster_A == cl)# %>% dplyr::select(Item,Time)
  item = pull(it,Item)
  time = pull(it,Time)
  return(cbind(item,time))
}

tabulate_id = function(chrid) {
  ## reference table of charactor and numeric id
  chr = sort(unique(chrid))
  out = data.frame(chr = chr, num = 1:length(chr))
  return(out)
}
to_numID = function(x, tab) {
  sapply(x, function(x) tab$num[which(tab$chr == x)])
}

to_chrID = function(x, tab) {
  sapply(x, function(x) tab$chr[which(tab$num == x)])
}

tab_sj = function(seg_g, G) {
  res = NULL
  for (m in 0:(G-1)) {
    res = c(res, sum(seg_g >= m))
  }
  return(res)
}

tab_IY = function(seg_g, G) {
  res = NULL
  for (m in 0:(G-1)) {
    res = c(res, sum(seg_g == m))
  }
  return(res)
}

tab_sj = function(seg_g, G) {
  res = NULL
  for (m in 0:(G-1)) {
    res = c(res, sum(seg_g >= m))
  }
  return(res)
}

tab_IY = function(seg_g, G) {
  res = NULL
  for (m in 0:(G-1)) {
    res = c(res, sum(seg_g == m))
  }
  return(res)
}

tabulate_id = function(chrid) {
  ## reference table of charactor and numeric id
  chr = sort(unique(chrid))
  out = data.frame(chr = chr, num = 1:length(chr))
  return(out)
}
to_numID = function(x, tab) {
  sapply(x, function(x) tab$num[which(tab$chr == x)])
}

to_chrID = function(x, tab) {
  sapply(x, function(x) tab$chr[which(tab$num == x)])
}

my_procrustes = function(Xstar, dlist, is_list = FALSE, translation = TRUE, scale = FALSE, reflect = TRUE) {
  posm = 0
  if (is_list == TRUE) {
    num_chain = length(dlist)
  } else { num_chain = 1 }
  for (i in 1:num_chain) {
    if (is_list == TRUE) { df = dlist[[i]]
    } else { df = dlist }

    num_samples = nrow(df)

    z0dx = grepl("^z\\.0\\.", colnames(df))
    z1dx = grepl("^z\\.1\\.", colnames(df))
    w0dx = grepl("^w\\.0\\.", colnames(df))
    w1dx = grepl("^w\\.1\\.", colnames(df))

    no_z1 = sum(z1dx) == 0
    no_w1 = sum(w1dx) == 0

    num_w = sum(w0dx) / 2;
    num_z = sum(z0dx) / 2;
    w0 = aperm( array(unlist( t(df[,w0dx])), dim = c(2, num_w, num_samples)), c(2,1,3))
    z0 = aperm( array(unlist( t(df[,z0dx])), dim = c(2, num_z, num_samples)), c(2,1,3))

    w0star = Xstar$w.0

    if (no_w1) {
      w1 = NULL
    } else {
      w1star = MCMCpack::procrustes(Xstar$w.1, w0star)$X.new
      w1 = aperm( array(unlist( t(df[,w1dx])), dim = c(2, num_w, num_samples)), c(2,1,3))
    }
    if (no_z1) {
      z1 = NULL
    } else {
      z1 = aperm( array(unlist( t(df[,z1dx])), dim = c(2, num_z, num_samples)), c(2,1,3))
    }

    adx = z0dx | z1dx | w0dx | w1dx

    mm = foreach (k = 1:num_samples) %dopar% {
      pout = MCMCpack::procrustes(w0[,,k], w0star)
      w0[,,k] = pout$X.new
      z0[,,k] = z0[,,k] %*% pout$R

      if (!no_w1) {
        pout = MCMCpack::procrustes(w1[,,k], w1star)
        w1[,,k] = pout$X.new
      }
      if (!no_z1) {
        z1[,,k] = z1[,,k] %*% pout$R
      }
      rbind(z0[,,k], z1[,,k], w0[,,k], w1[,,k])
    }
    tmm = lapply(mm,t)
    df[,adx] = t( matrix(unlist(tmm), nrow = sum(adx)) )

    posm = posm + Reduce("+",mm) / num_samples
    if (is_list == TRUE) { dlist[[i]] = df
    } else { dlist = df }
  }

  posm = posm / num_chain
  z0 = posm[1:num_z,]
  if (no_z1 && no_w1) {
    w0 = posm[num_z + (1:num_w),]
  } else if (no_z1 && !no_w1) {
    w0 = posm[num_z + (1:num_w),]
    w1 = posm[num_z + num_w + (1:num_w),]
  }

  if (!no_z1 && no_w1) {
    z1 = posm[num_z + (1:num_z),]
    w0 = posm[2*num_z + (1:num_w),]
  } else if (!no_z1 && !no_w1) {
    z1 = posm[num_z + (1:num_z),]
    w0 = posm[2*num_z + (1:num_w),]
    w1 = posm[2*num_z + num_w + (1:num_w),]
  }
  return(list(dlist = dlist, z0=z0, z1=z1, w0=w0, w1=w1))
}
