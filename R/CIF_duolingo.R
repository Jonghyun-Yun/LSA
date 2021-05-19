########### CIF for Duolingo data ############
out_dir <- "duolingo_pn_ncut5_zero_beta_noinfo_lc2/"
num_chain <- 2
HAS_REF <- 0
source("R/Renviron.R")
source("R/load-outputs.R")

tmp <- foreach(v = 1:num_chain, .combine = "rbind") %dopar% apply(mydf[[v]], 2, mean)
if (num_chain > 1) {
  tmp <- tmp[, 1:(which(colnames(mydf[[1]]) == "lp_"))]
  posm <- apply(tmp, 2, mean)
} else {
  posm <- tmp[1:(which(colnames(mydf[[1]]) == "lp_"))]
}

cname <- names(posm)
param <- getparam(posm, sj, i, k)

myN <- max(100, N) # 151
# sj <- c(0, 5, 8, 12, 19, 201)
maxt <- sj[G + 1] + 10
num_seg <- 100
time <- seq(0, maxt, (maxt) / num_seg)

# red item
item_red_ind <- c(1, 5, 7, 8, 10, 12, 13, 14, 15, 17)
person_red_sample <- c(1, 2, 7, 19, 34, 52, 70, 72, 85, 111) # row index

# blue item
item_blue_ind <- c(2, 3, 4, 6, 9, 11, 16, 18)
person_blue_sample <- c(10, 11, 12, 17, 24, 29, 47, 101, 141, 144) # row index

# visualizing latent position of selected subjects
z0 <- matched$z0
z1 <- matched$z1
w0 <- matched$w0
w1 <- matched$w1
xmin <- min(z0[, 1], z1[, 1], w0[, 1], w1[, 1])
ymin <- min(z0[, 2], z1[, 2], w0[, 2], w1[, 2])
xmax <- max(z0[, 1], z1[, 1], w0[, 1], w1[, 1])
ymax <- max(z0[, 2], z1[, 2], w0[, 2], w1[, 2])

z0 = z0[c(person_red_sample,person_blue_sample), ]

myname <- c(c(item_red_ind,item_blue_ind), paste0("I.", 1:I))
lsjmplot(z0, w0, xlim = c(xmin, xmax), ylim = c(ymin, ymax), myname)


library(ggplot2)
library(dplyr)

# temp
for (item in 1:2) {
  CIF_T <- foreach(k = c(person_red_sample, person_blue_sample), .combine =
                                                                   "rbind") %do% {
    c(0, cumcicurve(getparam(posm, sj, item, k), 1, 0, maxt, num_seg))
  }
  CIF_F <- foreach(k = c(person_red_sample, person_blue_sample), .combine =
                                                                   "rbind") %do% {
    c(0, cumcicurve(getparam(posm, sj, item, k), 0, 0, maxt, num_seg))
  }
  CIF_T <- data.frame(t(CIF_T), time)
  CIF_F <- data.frame(t(CIF_F), time)
  plot_T <- reshape::melt(CIF_T, id.vars = "time")
  plot_T$group <- c(rep('red',1010), rep('blue',1010))
  colnames(plot_T) <- c("time", "person", "CIF", "group")
  p <- ggplot(data = plot_T, aes(x = time, y = CIF, group = person)) +
    geom_line(aes(color = factor(group)), show.legend = FALSE) +
    scale_color_manual("plot_T",
                       breaks = c("red",'blue'),
                       values = c("red" = 'red', "blue" ='blue')) +
    ylim(0, 1) + ggtitle(paste0("Cumulative incidence function for TRUE response", item, ")."))
  print(p)
  
  plot_F <- reshape::melt(CIF_F, id.vars = c("time"))
  plot_F$group <- c(rep('red',1010),rep('blue',1010))
  colnames(plot_F) <- c("time", "person", "CIF", "group")
  p <- ggplot(data = plot_F, aes(x = time, y = CIF, group = person)) +
    geom_line(aes(color = factor(group)), show.legend = FALSE) +
    scale_color_manual("plot_T",
                       breaks = c("red", 'blue'),
                       values = c("red" = 'red', "blue" = "blue")) +
    ylim(0, 1) + ggtitle(paste0("Cumulative incidence function for FALSE response", item, ")."))
  print(p)
}


# CIF_total
pdf("cif_item_total_sub.pdf")
for (item in 1:18) {
  CIF_T <- foreach(k = c(person_red_sample, person_blue_sample), .combine = "rbind") %do% {
    c(0, cumcicurve(getparam(posm, sj, item, k), 1, 0, maxt, num_seg))
  }
  CIF_F <- foreach(k = c(person_red_sample, person_blue_sample), .combine = "rbind") %do% {
    c(0, cumcicurve(getparam(posm, sj, item, k), 0, 0, maxt, num_seg))
  }
  CIF_T <- data.frame(t(CIF_T), time)
  CIF_F <- data.frame(t(CIF_F), time)
  plot_T <- reshape::melt(CIF_T, id.vars = "time")
  plot_T$group <- c(rep('red',1010), rep('blue',1010))
  colnames(plot_T) <- c("time", "person", "CIF", "group")
  p <- ggplot(data = plot_T, aes(x = time, y = CIF, group = person)) +
    geom_line(aes(color = factor(group)), show.legend = FALSE) +
    scale_color_manual("plot_T",
                       breaks = c("red", 'blue'),
                       values = c("red" = 'red', "blue" ='blue')) +
    ylim(0, 1) + ggtitle(paste0("Cumulative incidence function for TRUE response (Item", item, ")."))
  print(p)
  
  plot_F <- reshape::melt(CIF_F, id.vars = c("time"))
  plot_F$group <- c(rep('red',1010), rep('blue',1010))
  colnames(plot_F) <- c("time", "person", "CIF", "group")
  p <- ggplot(data = plot_F, aes(x = time, y = CIF, group = person)) +
    geom_line(aes(color = factor(group)), show.legend = FALSE) +
    scale_color_manual("plot_T",
                       breaks = c("red", 'blue'),
                       values = c("red" = 'red', "blue" = "blue")) +
    ylim(0, 1) + ggtitle(paste0("Cumulative incidence function for FALSE response (Item", item, ")."))
  print(p)
}
dev.off(which = dev.cur())




