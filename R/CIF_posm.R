tmp = foreach(v=1:num_chain, .combine='rbind') %dopar% apply(mydf[[v]], 2, mean)
if (num_chain > 1) {
  tmp = tmp[,1:( which( colnames(mydf[[1]]) == "lp_") )]
  posm = apply(tmp, 2, mean)
} else {
  posm = tmp[1:( which( colnames(mydf[[1]]) == "lp_") )]
}

cname = names(posm)
param = getparam(posm,sj,i,k)
