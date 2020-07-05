library(ggplot2)
library(ggrepel)

lsjmplot <- function( z, w, myname = NULL, xlim=NA, ylim=NA, lab = "Coordinate") {

  ## extract objects

  x = rbind(z,w)
  idx = rep("w", nrow(x))
  idx[1:nrow(z)] = "z"
  position <- as.data.frame(x)
  ndim <- dim(x)[2]

  colnames(position) <- paste("position",1:ndim,sep="")

  padding = 1.05
  if (any(is.na(xlim))) {
    x1 <- -max(abs(position[,1]))*padding
    x2 <- max(abs(position[,1]))*padding
  } else {
    x1 <- xlim[1]
    x2 <- xlim[2]
  }
  if (any(is.na(ylim))) {
    y1 <- -max(abs(position[,2]))*padding
    y2 <- max(abs(position[,2]))*padding
  } else {
    y1 <- ylim[1]
    y2 <- ylim[2]
  }

  mytheme = theme(axis.line = element_line(colour = "black"),
                  ##panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  ##panel.border = element_blank(),
                  panel.background = element_blank()
                  )

  ## plot
  pp = ggplot(position,aes(x=position1,y=position2,colour=idx)) +
    theme(text=element_text(size=20)) +
    ## geom_point()+
    xlim(x1,x2) + ylim(y1,y2) +
    xlab(paste(lab," 1",sep="")) + ylab(paste(lab," 2",sep="")) +
    ##xlab("Position 1") + ylab("Position 2") +
    geom_hline(yintercept = 0, color = "gray70", linetype=2) +
    geom_vline(xintercept = 0, color = "gray70", linetype=2)
  ##  pp = pp + geom_text_repel(label=rownames(position), segment.color = "grey50", size=6)
  if (!is.null(myname)) {
    pp = pp + geom_text(label=myname,
                        ## segment.color = "grey50",
                        check_overlap = FALSE, show.legend=FALSE,size = 2)
  } else pp = pp + geom_point()
  pp + mytheme
}
