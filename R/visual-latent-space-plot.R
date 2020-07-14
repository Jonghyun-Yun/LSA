z0.star = Xstar$z.0
z1.star = Xstar$z.1
w0.star = Xstar$w.0
w1.star = Xstar$w.1

xmin = min(z0.star[,1],z1.star[,1],w0.star[,1],w1.star[,1])
ymin = min(z0.star[,2],z1.star[,2],w0.star[,2],w1.star[,2])
xmax = max(z0.star[,1],z1.star[,1],w0.star[,1],w1.star[,1])
ymax = max(z0.star[,2],z1.star[,2],w0.star[,2],w1.star[,2])

myname = c(1:N,paste0("I.",1:I))
pdf("figure/star_latent.pdf")
print(lsjmplot(z0.star,w0.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (!single_z && single_w) {
  print(lsjmplot(z1.star,w0.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!single_z && !single_w) {
  print(lsjmplot(z1.star,w1.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}

dev.off(which = dev.cur())

pdf("figure/latent_position_plot_vegan.pdf")

z0 = matched$z0
z1 = matched$z1
w0 = matched$w0
w1 = matched$w1
xmin = min(z0[,1],z1[,1],w0[,1],w1[,1])
ymin = min(z0[,2],z1[,2],w0[,2],w1[,2])
xmax = max(z0[,1],z1[,1],w0[,1],w1[,1])
ymax = max(z0[,2],z1[,2],w0[,2],w1[,2])

myname = c(1:N,paste0("I.",1:I))
print(lsjmplot(z0,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (!single_z && single_w) {
  print(lsjmplot(z1,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!single_z && !single_w) {
  print(lsjmplot(z1,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}
dev.off(which = dev.cur())
