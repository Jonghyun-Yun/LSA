mvar = read_csv("mvar.csv", col_names=F) %>% as.matrix()
I = mvar[1]
N = mvar[2]
G = mvar[4]

sj = read_csv("mlen.csv",col_names=F) %>% as.matrix()
sj = c(0, cumsum(sj))
H = sj[2:(G+1)] - sj[1:G]

cnames = c(".chain", ".iteration")
for (c in 0:1) {
  for (i in 1:I) {
    for (g in 1:G) {
      cnames = c(cnames, paste0("lambda.",c,".",i,".",g))
    }}}

for (i in 1:I) {
  for (c in 0:1) {
    cnames = c(cnames, paste0("beta.",i,".",c))
  }}

for (k in 1:N) {
  for (c in 0:1) {
    cnames = c(cnames, paste0("theta.",k,".",c))
  }}

for(c in 0:1)
  for (k in 1:N) {
    for (d in 1:2) {
      cnames = c(cnames, paste0("z.",c,".",k,".",d))
    }}
for (i in 1:I) {
  for (d in 1:2) {
    cnames = c(cnames, paste0("w.",i,".",d))
  }}

for (c in 0:1) {
  cnames = c(cnames, paste0("gamma.",c))
}

cnames = c(cnames, "sigma", "lp_")

## mythin = 10
## mystart = 5001
## myend = 25000
mydf = list()
mylist = mcmc.list()
for (cid in 1:5) {
mydf[[cid]] = read_csv(paste0("output/sample_chain",cid,".csv"), col_names=F) %>% as.data.frame()
colnames(mydf[[cid]]) = cnames
}

## mylist[[cid]] = mcmc(df, start = mystart, end = myend, thin = mythin)

Xstar = find_xstar_inlist(mydf)
temp = do_procrustes(Xstar, mydf, is_list = TRUE)

mydf = temp$mydf
z0 = temp$z0
z1 = temp$z1
w = temp$w
xmin = min(z0[,1],z1[,1],w[,1])
ymin = min(z0[,2],z1[,2],w[,2])
xmax = max(z0[,1],z1[,1],w[,1])
ymax = max(z0[,2],z1[,2],w[,2])

mylist = mcmc.list()
for (cid in 1:5) {
mylist[[cid]] = mcmc(mydf[[cid]])
}
