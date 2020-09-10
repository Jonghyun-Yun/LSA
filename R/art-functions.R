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

getparam = function(posm, sj, i, k) {
  cname = names(posm)
  z = posm[str_which(cname, paste0("z\\.[0-1]\\.",k,"\\.[1-2]"))] %>% matrix(nrow = 2, ncol = 2) %>% t()
  w = posm[str_which(cname, paste0("w\\.[0-1]\\.",i,"\\.[1-2]"))] %>% matrix(nrow = 2, ncol = 2) %>% t()
  gamma = posm[str_which(cname, paste0("gamma"))]
  beta = posm[str_which(cname, paste0("beta\\.",i,"\\."))]
  theta = posm[str_which(cname, paste0("theta\\.",k,"\\."))]
  lambda = posm[str_which(cname, paste0("lambda\\.[0-1]\\.",i,"\\."))] %>% matrix(ncol = 2) %>% t()
  H = sj[2:(G+1)] - sj[1:G]
  res = list(lambda=lambda,beta=beta,theta=theta,gamma=gamma,z=z,w=w,sj=sj,H=H)
  if (any(unlist(lapply(res, is.na)))) stop("Index out of range")
  ## should if be G+1? or G?
  ## if (ncol(lambda) != (G+1)) stop("ncol(lambda) != G+1")
  if (ncol(lambda) != (G)) stop("ncol(lambda) != G")
  else return(res)
}

fun_hazard_surv = function(t,i,k,posm,cname,sj) {
  z = posm[str_which(cname, paste0("z\\.[0-1]\\.",k,"\\.[1-2]"))] %>% matrix(ncol = 2)
  w = rep(posm[str_which(cname, paste0("w\\.",i,"\\."))], 2)  %>% matrix(ncol = 2)
  gamma = posm[str_which(cname, paste0("gamma"))] %>% matrix(ncol = 2)
  beta = posm[str_which(cname, paste0("beta\\.",i,"\\."))] %>% matrix(ncol = 2)
  theta = posm[str_which(cname, paste0("theta\\.",k,"\\."))] %>% matrix(ncol = 2)
  lambda = posm[str_which(cname, paste0("lambda\\.[0-1]\\.",i,"\\."))] %>% matrix(ncol = 2)

  G = length(lambda[,1])
  H = sj[2:(G+1)] - sj[1:G]

  seg = 0
  for (g in 1:G) {
    seg = seg + 1 * (t > sj[g])
  }
  out = lambda[seg,2] * exp(beta[,2] + theta[,2] - gamma[,2] * sqrt(sum((z[,2]-w[,2])^2)))
  if (seg == 1) {
    for (c in 1:2) {
      out = out * exp(
                    - ((t - sj[seg]) *lambda[seg,c]) * exp(beta[,c] + theta[,c] - gamma[,c] * sqrt( sum((z[,c]-w[,c])^2))))
    }
  } else {
    for (c in 1:2) {
      out = out * exp(
                    - ((t - sj[seg]) *lambda[seg,c] + sum(H[1:(seg-1)] * lambda[1:(seg-1),c])) * exp(beta[,c] + theta[,c] - gamma[,c] * sqrt(sum((z[,c]-w[,c])^2))))
    }
  }
  names(out) = NULL
  return(out)
}

fun_hazard_ick = function(t,i,c,k,posm,cname,sj) {
  z = posm[str_which(cname, paste0("z\\.",c,"\\.",k,"\\.[1-2]"))]
  w = posm[str_which(cname, paste0("w\\.",i))]
  gamma = posm[str_which(cname, paste0("gamma\\.",c))]
  beta = posm[str_which(cname, paste0("beta\\.",i,"\\.",c))]
  theta = posm[str_which(cname, paste0("theta\\.",k,"\\.",c))]
  lambda = posm[str_which(cname, paste0("lambda\\.",c,"\\.",i,"\\."))]

  G = length(lambda)
  seg = 0 * t
  for (g in 1:G)
    seg = seg + 1 * (t > sj[g])
  hazard = lambda[seg] * exp(beta + theta - gamma * sum((z-w)^2))
  names(hazard) = NULL
  return(hazard)
}

fun_accuracy_ick = function(t,i,k,posm,cname,sj) {
  fun_hazard_ick(t,i,1,k,posm,cname,sj) / (fun_hazard_ick(t,i,1,k,posm,cname,sj) + fun_hazard_ick(t,i,0,k,posm,cname,sj))
}

library(ggplot2)
library(ggrepel)

lsjmplot <- function( z, w, myname = NULL, xlim=NA, ylim=NA, lab = "Coordinate") {

  ## extract objects

  x = rbind(z,w)
  idx = rep("w", nrow(x))
  idx[1:nrow(z)] = "z"
  position <- as.data.frame(x)
  ndim <- dim(x)[2]

  colnames(position) <- paste("position",1:ndim,sep="")

  padding = 1.05
  if (any(is.na(xlim))) {
    x1 <- -max(abs(position[,1]))*padding
    x2 <- max(abs(position[,1]))*padding
  } else {
    x1 <- xlim[1]
    x2 <- xlim[2]
  }
  if (any(is.na(ylim))) {
    y1 <- -max(abs(position[,2]))*padding
    y2 <- max(abs(position[,2]))*padding
  } else {
    y1 <- ylim[1]
    y2 <- ylim[2]
  }

  mytheme = theme(axis.line = element_line(colour = "black"),
                  ##panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  ##panel.border = element_blank(),
                  panel.background = element_blank()
                  )

  ## plot
  pp = ggplot(position,aes(x=position1,y=position2,colour=idx)) +
    theme(text=element_text(size=20)) +
    ## geom_point()+
    xlim(x1,x2) + ylim(y1,y2) +
    xlab(paste(lab," 1",sep="")) + ylab(paste(lab," 2",sep="")) +
    ##xlab("Position 1") + ylab("Position 2") +
    geom_hline(yintercept = 0, color = "gray70", linetype=2) +
    geom_vline(xintercept = 0, color = "gray70", linetype=2)
  ##  pp = pp + geom_text_repel(label=rownames(position), segment.color = "grey50", size=6)
  if (!is.null(myname)) {
    pp = pp + geom_text(label=myname,
                        ## segment.color = "grey50",
                        check_overlap = FALSE, show.legend=FALSE,size = 2)
  } else pp = pp + geom_point()
  pp + mytheme
}

find_xstar = function(df) {
  num_samples = nrow(df)

  z0dx = grepl("^z\\.0\\.", colnames(df))
  z1dx = grepl("^z\\.1\\.", colnames(df))
  w0dx = grepl("^w\\.0\\.", colnames(df))
  w1dx = grepl("^w\\.1\\.", colnames(df))
  adx = z0dx | z1dx | w0dx | w1dx

  mlp_ = max(df$lp_)
  star = min(which.max(df$lp_))
  lpos = df[,adx]
  Xstar = list()
  Xstar$z.0 = matrix(unlist(df[star,z0dx]), byrow = T, ncol = 2)
  if (sum(z1dx) == 0) {
    Xstar$z.1 = NULL
  } else {
    Xstar$z.1 = matrix(unlist(df[star,z1dx]), byrow = T, ncol = 2)
  }
  if (sum(w1dx) == 0) {
    Xstar$w.1 = NULL
  } else {
    Xstar$w.1 = matrix(unlist(df[star,w1dx]), byrow = T, ncol = 2)
  }
  Xstar$w.0 = matrix(unlist(df[star,w0dx]), byrow = T, ncol = 2)
  return(list(lp_ = mlp_, Xstar=Xstar))}

find_xstar_inlist = function(mydf) {
  num_chain = length(mydf)
  mlp = -Inf
  for (i in 1:num_chain) {
    slist = find_xstar(mydf[[i]])
    if (slist$lp_ > mlp) Xstar = slist$Xstar
  }
  return(Xstar)}

do_procrustes = function(Xstar, mydf, is_list = FALSE, translation = TRUE, scale = FALSE, reflect = TRUE) {
  posm = 0
  if (is_list == TRUE) {
    num_chain = length(mydf)
  } else { num_chain = 1 }
  for (i in 1:num_chain) {
    if (is_list == TRUE) { df = mydf[[i]]
    } else { df = mydf }

    num_samples = nrow(df)

    z0dx = grepl("^z\\.0\\.", colnames(df))
    z1dx = grepl("^z\\.1\\.", colnames(df))
    wdx = grepl("^w", colnames(df))
    adx = z0dx | z1dx | wdx
    N = sum(z0dx) / 2
    nall = sum(adx)/2

    mlp_ = max(df$lp_)
    star = min(which.max(df$lp_))
    lpos = df[,adx]

    ## mm = list()
    ## for (k in 1:num_samples) {
    ##   X = matrix(unlist(lpos[k,]), nrow = 2) %>% t()
    ##   ## mm[[k]] = MCMCpack::procrustes(X, Xstar, translation, dilation)$X.new #MCMCpack
    ##   mm[[k]] = vegan::procrustes(X, Xstar, scale = scale)$Yrot #vegan
    ##   df[k,adx] = mm[[k]] %>% t() %>% c()
    ## }

    mm = foreach (k = 1:num_samples) %dopar% {
      ## X = matrix(unlist(lpos[k,]), nrow = 2) %>% t()
      ## mm[[k]] = MCMCpack::procrustes(X, Xstar, translation, dilation = scale)$X.new #MCMCpack
      ## vegan::procrustes(Xstar, t( matrix(unlist(lpos[k,]), nrow = 2) ), scale = scale)$Yrot #vegan
      shapes::procOPA(Xstar, t( matrix(unlist(lpos[k,]), nrow = 2)) , scale = scale, reflect = reflect)$Bhat #shapes
    }
    tmm = lapply(mm,t)
    df[,adx] = t( matrix(unlist(tmm), nrow = sum(adx)) )


    posm = posm + Reduce("+",mm) / num_samples
    if (is_list == TRUE) { mydf[[i]] = df
    } else { mydf = df }
  }

  posm = posm / num_chain
  z0= posm[1:N,]
  if (sum(z1dx) == 0) {
    w = posm[-(1:N),]
    z1 = NULL
  } else {
    z1 = posm[(N + 1):(2*N),]
    w = posm[-(1:(2*N)),]
  }
  return(list(mydf = mydf, z0=z0, z1=z1, w=w))
}
