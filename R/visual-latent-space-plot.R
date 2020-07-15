mcdf = mcmc(mdf[,-1])
tn = nrow(mcdf)
hpd = apply(HPDinterval(mcdf, prob = 1/tn), 1, median)

z0dx = grepl("^z\\.0\\.", names(hpd))
z1dx = grepl("^z\\.1\\.", names(hpd))
w0dx = grepl("^w\\.0\\.", names(hpd))
w1dx = grepl("^w\\.1\\.", names(hpd))

num_w = sum(w0dx) / 2;
num_z = sum(z0dx) / 2;

no_z1 = sum(z1dx) == 0
no_w1 = sum(w1dx) == 0
w0 = matrix(hpd[w0dx], byrow = T, nrow = num_w, ncol = 2)
if (no_w1) {
  w1 = NULL
} else w1 = matrix(hpd[w1dx], byrow = T, nrow = num_w, ncol = 2)
if (no_z1) {
  z1 = NULL
} else z1 = matrix(hpd[z1dx], byrow = T, nrow = num_z, ncol = 2)

z0 = matrix(hpd[z0dx], byrow = T, nrow = num_z, ncol = 2)
xmin = min(z0[,1],z1[,1],w0[,1],w1[,1])
ymin = min(z0[,2],z1[,2],w0[,2],w1[,2])
xmax = max(z0[,1],z1[,1],w0[,1],w1[,1])
ymax = max(z0[,2],z1[,2],w0[,2],w1[,2])
myname = c(1:N,paste0("I.",1:I))

pdf("figure/latent_space_hpd.pdf")
print(lsjmplot(z0,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (double_z && !double_w) {
  print(lsjmplot(z1,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (double_z && double_w) {
  print(lsjmplot(z1,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!double_z && double_w) {
  print(lsjmplot(z0,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}
dev.off(which = dev.cur())

med = apply(mdf[,-1] , 2, median)

z0dx = grepl("^z\\.0\\.", names(med))
z1dx = grepl("^z\\.1\\.", names(med))
w0dx = grepl("^w\\.0\\.", names(med))
w1dx = grepl("^w\\.1\\.", names(med))

num_w = sum(w0dx) / 2;
num_z = sum(z0dx) / 2;

no_z1 = sum(z1dx) == 0
no_w1 = sum(w1dx) == 0
w0 = matrix(med[w0dx], byrow = T, nrow = num_w, ncol = 2)
if (no_w1) {
  w1 = NULL
} else w1 = matrix(med[w1dx], byrow = T, nrow = num_w, ncol = 2)
if (no_z1) {
  z1 = NULL
} else z1 = matrix(med[z1dx], byrow = T, nrow = num_z, ncol = 2)

z0 = matrix(med[z0dx], byrow = T, nrow = num_z, ncol = 2)
xmin = min(z0[,1],z1[,1],w0[,1],w1[,1])
ymin = min(z0[,2],z1[,2],w0[,2],w1[,2])
xmax = max(z0[,1],z1[,1],w0[,1],w1[,1])
ymax = max(z0[,2],z1[,2],w0[,2],w1[,2])
myname = c(1:N,paste0("I.",1:I))

pdf("figure/latent_space_med.pdf")
print(lsjmplot(z0,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (double_z && !double_w) {
  print(lsjmplot(z1,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (double_z && double_w) {
  print(lsjmplot(z1,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!double_z && double_w) {
  print(lsjmplot(z0,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}
dev.off(which = dev.cur())

z0.star = Xstar$z.0
z1.star = Xstar$z.1
w0.star = Xstar$w.0
w1.star = Xstar$w.1

xmin = min(z0.star[,1],z1.star[,1],w0.star[,1],w1.star[,1])
ymin = min(z0.star[,2],z1.star[,2],w0.star[,2],w1.star[,2])
xmax = max(z0.star[,1],z1.star[,1],w0.star[,1],w1.star[,1])
ymax = max(z0.star[,2],z1.star[,2],w0.star[,2],w1.star[,2])

myname = c(1:N,paste0("I.",1:I))
pdf("figure/latent_space_mode.pdf")
print(lsjmplot(z0.star,w0.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
if (double_z && !double_w) {
  print(lsjmplot(z1.star,w0.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (double_z && double_w) {
  print(lsjmplot(z1.star,w1.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
} else if (!double_z && double_w) {
  print(lsjmplot(z0.star,w1.star,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
}

dev.off(which = dev.cur())

pdf("figure/latent_space_mean.pdf")

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
 if (double_z && !double_w) {
   print(lsjmplot(z1,w0,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
 } else if (double_z && double_w) {
   print(lsjmplot(z1,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
 } else if (!double_z && double_w) {
   print(lsjmplot(z0,w1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),myname))
 }
dev.off(which = dev.cur())
