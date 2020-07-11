mvar = readr::read_csv("input/mvar.csv", col_names=F) %>% as.matrix()
I = mvar[1]
N = mvar[2]
G = mvar[4]

sj = readr::read_csv("input/mlen.csv",col_names=F) %>% as.matrix()
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

for(c in 0:1) {
  for (k in 1:N) {
    for (d in 1:2) {
      cnames = c(cnames, paste0("z.",c,".",k,".",d))
    }}
  }
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
for (cid in 1:num_chain) {
mydf[[cid]] = readr::read_csv(paste0("output/sample_chain",cid,".csv"), col_names=F) %>% as.data.frame()
colnames(mydf[[cid]]) = cnames
if (single_z) {
  mydf[[cid]] = mydf[[cid]][,!grepl("^z\\.1\\.", cnames)] ## remove duplicates when single_z
}
}

## mylist[[cid]] = mcmc(df, start = mystart, end = myend, thin = mythin)

Xstar = find_xstar_inlist(mydf)
matched = do_procrustes(Xstar, mydf, is_list = TRUE)
mydf = matched$mydf

item = 1
 c = 0
 cname = names(mydf[[1]])
 mylist = mcmc.list()
 for (cid in 1:num_chain) {
    for (k in 1:N) {
       z = mydf[[cid]][,str_which(cname, paste0("z\\.",c,"\\.",k,"\\.[1-2]"))]
       w = mydf[[cid]][,str_which(cname, paste0("w\\.",item,"\\."))]
       mydf[[cid]][[paste0("dist_z.",c,".",k,"_","w.",item)]] = sqrt(rowSums((z-w)^2))
}
   mylist[[cid]] = mcmc(mydf[[cid]])
 }
