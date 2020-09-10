## choose a few people adn items from latent space
z0 = matched$z0
w0 = matched$w0
## myI = sort(c(35,7,6,24,28,29,8,27,13,15)) ## opus
## myI = sort(c(32,24,30,8,2,15,6,12,5,20,7,26,21,31)) ## verbal
## myI = sort(c(17,14,3,2,15,9,10,13,11)) ## pisa
## myI = sort(c(8,9,5,4,13,2,3,15)) ## pisa KR sci
myI = sort(c(8,12,14,17,1,6,13,18)) ## pisa KR

which_z = function(w, z) {
  which.min(colSums((t(z) - w)^2))
}

myN = NULL
for (i in 1:length(myI)) {
  myN = c(myN, which_z(w0[myI[i],],z0))
}
myN = c(myN, which.min(rowSums(z0^2)), which.max(rowSums(z0^2)))

maxt = sj[G+1] + 10
num_seg = 100
time = seq(0, maxt, (maxt) / num_seg)

pdf(paste0("figure/CIF_person.pdf"))

for (k in myN) {

  CIF_T = foreach (item=myI, .combine='rbind') %do% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 1, 0, maxt, num_seg))
  }

  CIF_F = foreach (item=myI, .combine='rbind') %dopar% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 0, 0, maxt, num_seg))
  }

  CIF_T = data.frame(t(CIF_T), c(time))
  colnames(CIF_T) = c(paste0("I",myI),"time")

  plot_T <- reshape2::melt(CIF_T, id.vars="time")
  colnames(plot_T) = c("time", "item", "CIF")

  p = ggplot(data=plot_T, aes(x=time, y=CIF, group=item)) +
    geom_line(aes(color=factor(item)), show.legend=TRUE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ##scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for TRUE response (person ",k,")."))
  print(p)

  CIF_F = data.frame(t(CIF_F), c(time))
  colnames(CIF_F) = c(paste0("I",myI),"time")

  plot_F <- reshape2::melt(CIF_F, id.vars=c("time"))
  colnames(plot_F) = c("time", "item", "CIF")

  p = ggplot(data=plot_F, aes(x=time, y=CIF, group=item)) +
    geom_line(aes(color=factor(item)), show.legend=TRUE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ## scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for FALSE response (person ",k,")."))
  print (p)

}
dev.off(which = dev.cur())

pdf("figure/latent_space_person.pdf")

z1 = matched$z1
w1 = matched$w1
z0 = matched$z0[myN,]
w0 = matched$w0[myI,]

xmin = min(z0[,1],z1[,1],w0[,1],w1[,1])
ymin = min(z0[,2],z1[,2],w0[,2],w1[,2])
xmax = max(z0[,1],z1[,1],w0[,1],w1[,1])
ymax = max(z0[,2],z1[,2],w0[,2],w1[,2])

myname = c(myN,paste0("I.",myI))
print(lsjmplot(z0,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (double_z && !double_w) {
  print(lsjmplot(z1,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (double_z && double_w) {
  print(lsjmplot(z1,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!double_z && double_w) {
  print(lsjmplot(z0,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}
dev.off(which = dev.cur())

system(paste0("rsync -v figure/CIF_person.pdf ", out_dir, "figure/"))
system(paste0("rsync -v figure/latent_space_person.pdf ", out_dir, "figure/"))
