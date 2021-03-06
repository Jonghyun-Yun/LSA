#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 2) {
  # default output file
  num_chain <- args[2]
  HAS_REF <- 0
  ref_dir <- ""
} else if (length(args) == 1) {
  # default output file
  num_chain <- 3
  HAS_REF <- 0
  ref_dir <- ""
}
out_dir <- args[1]
double_z <- 0
double_w <- 0
## out_dir <- "chessB-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("rsync -rv figure/*.pdf ", out_dir, "figure/"))
