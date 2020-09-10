pdf("figure/latent_space_cl1.pdf")

## myI = sort(c(35,7,6,24,28,29,8,27,13,15)) ## opus
## myI = sort(c(32,24,30,8,2,15,6,12,5,20,7,26,21,31)) ## verbal
## myI = sort(c(17,14,3,2,15,9,10,13,11)) ## pisa
## myI = sort(c(8,9,5,4,13,2,3,15)) ## pisa KR sci

## myN = c(23,37,20,27,33,36) ## marketing
## myI = 1:I # marketing

## myN = c(62,589,14,219,524,509,161,435,252,595) ## pisa KR
## myI = 1:I # pisa KR

## myN = c(518, 584, 370, 418, 79, 276, 572, 358, 73, 274) ## pisa
## myI = 1:I # pisa

myN = c(270,367,712,17,615,265,653,484,155,290) ## pisa KR sci
myI = 1:I # pisa KR sci

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
