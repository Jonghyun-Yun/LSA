## accuracy = foreach(k=1:N, .combine='rbind') %dopar% fun_accuracy_ick(t,i,k,posm,cname,sj)
pdf(paste0("figure/tradeoff.pdf"))
for (item in 1:I) {
  time = 1:(sj[G] + 10)
  accuracy = foreach(k=1:N, .combine='rbind') %dopar%
    {
      param = getparam(posm,sj,item,k)
      eval_accuracy(param, time)
    }

  plotdf = reshape::melt(accuracy, id.vars=c("time"))
  colnames(plotdf) = c("person", "time", "accuracy")

  p = ggplot(data=plotdf, aes(x=time, y=accuracy, group=person)) +
    geom_line(aes(color=factor(person)), show.legend=FALSE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ##scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Plot of speed-accuracy tradeoff for item ",item,"."))
  print(p)
}
dev.off(which = dev.cur())

## system(paste0("open figure/tradeoff_", item, ".pdf")

myN = min(100, N)
maxt = sj[G+1] + 10
num_seg = 100
time = seq(0, maxt, (maxt) / num_seg)

pdf(paste0("figure/CIF.pdf"))

for (item in 1:I) {

  CIF_T = foreach (k=1:myN, .combine='rbind') %do% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 1, 0, maxt, num_seg))
  }

  CIF_F = foreach (k=1:myN, .combine='rbind') %dopar% {
    c(0,cumcicurve(getparam(posm, sj, item, k), 0, 0, maxt, num_seg))
  }

  CIF_T = data.frame(t(CIF_T), time)
  CIF_F = data.frame(t(CIF_F), time)
  plot_T <- reshape::melt(CIF_T, id.vars="time")
  colnames(plot_T) = c("time", "person", "CIF")
  p = ggplot(data=plot_T, aes(x=time, y=CIF, group=person)) +
    geom_line(aes(color=factor(person)), show.legend=FALSE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ##scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for TRUE response (item ",item,")."))
  print(p)

  plot_F <- reshape::melt(CIF_F, id.vars=c("time"))
  colnames(plot_F) = c("time", "person", "CIF")
  p = ggplot(data=plot_F, aes(x=time, y=CIF, group=person)) +
    geom_line(aes(color=factor(person)), show.legend=FALSE) +
    ylim(0,1) +
    ## geom_line(group="3", col="red") +
    ## scale_colour_grey() +
    theme_bw() +
    ggtitle(paste0("Cumulative incidence function for FALSE response (item ",item,")."))
  print (p)
}
dev.off(which = dev.cur())
