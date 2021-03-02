out_dir = "chessB_np/"
load(paste0(out_dir,fname))
dll = t(mll)
colnames(dll) = row.names(mi) = 1:I
row.names(dll) = colnames(mi) = 1:N
d_long <- reshape2::melt(dll)
mi_long <- reshape2::melt(t(mi))
d_long = plyr::join(d_long, mi_long, by = c("Var1","Var2"))

out_dir = "chessB_no_latent/"
load(paste0(out_dir,fname))
dll = t(mll)
colnames(dll) = 1:I
row.names(dll) = 1:N
no_long <- reshape2::melt(dll)
d_long = plyr::join(d_long, no_long, by = c("Var1","Var2"))

out_dir = "chessB_pn/"
load(paste0(out_dir,fname))
dll = t(mll)
colnames(dll) = 1:I
row.names(dll) = 1:N
np_long <- reshape2::melt(dll)
d_long = plyr::join(d_long, np_long, by = c("Var1","Var2"))

names(d_long) = c("person","item","LogLoss_np","res", "LogLoss_no", "LogLoss_pn")
d_long$LogLoss_diff_np = d_long$LogLoss_no - d_long$LogLoss_np
d_long$LogLoss_diff_pn = d_long$LogLoss_no - d_long$LogLoss_pn

d_long$res = factor(d_long$res, labels=c("incorrect","correct"))
ll_boxp_np_by <- ggplot(d_long, aes(x=item,y=LogLoss_diff_np,fill=factor(item))) + facet_wrap(~res) +
  geom_boxplot() + theme(legend.position = "none") + ggtitle("baseline - (-,+) log loss")+
  ylab("Difference in LogLoss")

ll_boxp_pn_by <- ggplot(d_long, aes(x=item,y=LogLoss_diff_pn,fill=factor(item))) + facet_wrap(~res) +
  geom_boxplot() + theme(legend.position = "none") + ggtitle("baseline - (+,-) log loss") +
    ylab("Difference in LogLoss")

ll_boxp_np <- ggplot(d_long, aes(x=item,y=LogLoss_diff_np,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ggtitle("baseline - (-,+) log loss")+
  ylab("Difference in LogLoss")

ll_boxp_pn <- ggplot(d_long, aes(x=item,y=LogLoss_diff_pn,fill=factor(item))) +
  geom_boxplot() + theme(legend.position = "none") + ggtitle("baseline - (+,-) log loss") +
    ylab("Difference in LogLoss")

pdf(pname)
print(ll_boxp_np_by)
print(ll_boxp_pn_by)
print(ll_boxp_np)
print(ll_boxp_pn)
dev.off()
