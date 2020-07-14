num_chain = 3
single_z = 0
single_w = 0

library(art)
library(coda)
library(dplyr)
library(stringr)
library(bayesplot)
library(foreach)
library(doParallel)
registerDoParallel(cores = detectCores() - 1)

setwd("/Users/yunj/Dropbox/research/lsjm-art/lsjm-code")

source("R/art-functions.R")

source('R/load-outputs.R')
source('R/visual-latent-space-plot.R')
source('R/visual-trace-plot.R')
source('R/CIF.R')
