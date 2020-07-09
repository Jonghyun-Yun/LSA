z0.star = Xstar[1:N,]
 z1.star = Xstar[(N+1):(2*N),]
 w.star = Xstar[(2*N + 1):nrow(Xstar),]
 xmin = min(z0.star[,1],z1.star[,1],w.star[,1])
 ymin = min(z0.star[,2],z1.star[,2],w.star[,2])
 xmax = max(z0.star[,1],z1.star[,1],w.star[,1])
 ymax = max(z0.star[,2],z1.star[,2],w.star[,2])

myname = c(1:N,paste0("I.",1:I))
pdf("figure/star_latent.pdf")
 lsjmplot(z0.star,w.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
 lsjmplot(z1.star,w.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
dev.off(which = dev.cur())

z0.star = Xstar[1:N,]
 w.star = Xstar[(N + 1):nrow(Xstar),]
 xmin = min(z0.star[,1],w.star[,1])
 ymin = min(z0.star[,2],w.star[,2])
 xmax = max(z0.star[,1],w.star[,1])
 ymax = max(z0.star[,2],w.star[,2])

myname = c(1:N,paste0("I.",1:I))
pdf("figure/star_latent.pdf")
 lsjmplot(z0.star,w.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
dev.off(which = dev.cur())

pdf("figure/latent_position_plot_vegan.pdf")

z0 = matched$z0
z1 = matched$z1
w = matched$w %>% as.matrix()
xmin = min(z0[,1],z1[,1],w[,1])
ymin = min(z0[,2],z1[,2],w[,2])
xmax = max(z0[,1],z1[,1],w[,1])
ymax = max(z0[,2],z1[,2],w[,2])

myname = c(1:N,paste0("I.",1:I))
lsjmplot(z0,w,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
lsjmplot(z1,w,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
dev.off(which = dev.cur())

pdf("figure/latent_position_plot_vegan.pdf")

z0 = matched$z0
w = matched$w %>% as.matrix()
xmin = min(z0[,1],w[,1])
ymin = min(z0[,2],w[,2])
xmax = max(z0[,1],w[,1])
ymax = max(z0[,2],w[,2])

myname = c(1:N,paste0("I.",1:I))
lsjmplot(z0,w,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname)
dev.off(which = dev.cur())
