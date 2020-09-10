## choose a few people adn items from latent space

which_z = function(w, z) {
  which.min(colSums((t(z) - w)^2))
}

## myN = NULL
## for (i in 1:length(myI)) {
## myN = c(myN, which_z(w0[myI[i],],z0))
## }
## myN = c(myN, which.min(rowSums(z0^2)), which.max(rowSums(z0^2)))

maxt = sj[G] + 500
num_seg = 100
time = seq(0, maxt, (maxt) / num_seg)

pdf(paste0("figure/CIF_cl1.pdf"))

for (item in myI) {

  CIF_T = foreach (k = myN, .combine='rbind') %dopar% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 1, 0, maxt, num_seg))
  }

  CIF_F = foreach (k = myN, .combine='rbind') %dopar% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 0, 0, maxt, num_seg))
  }

  CIF_T = data.frame(t(CIF_T), c(time))
  colnames(CIF_T) = c(myN,"time")

  plot_T <- reshape2::melt(CIF_T, id.vars="time")
  colnames(plot_T) = c("time", "person", "CIF")

  p = ggplot(data=plot_T, aes(x=time, y=CIF, group=person)) +
    geom_line(aes(color=factor(person)), show.legend=TRUE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ##scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for TRUE response (item ",item,")."))
  print(p)

  CIF_F = data.frame(t(CIF_F), c(time))
  colnames(CIF_F) = c(myN,"time")

  plot_F <- reshape2::melt(CIF_F, id.vars=c("time"))
  colnames(plot_F) = c("time", "person", "CIF")

  p = ggplot(data=plot_F, aes(x=time, y=CIF, group=person)) +
    geom_line(aes(color=factor(person)), show.legend=TRUE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ## scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for FALSE response (item ",item,")."))
  print (p)

}
dev.off(which = dev.cur())

system(paste0("rsync -v figure/CIF_cl1.pdf ", out_dir, "figure/"))
system(paste0("rsync -v figure/latent_space_cl1.pdf ", out_dir, "figure/"))
