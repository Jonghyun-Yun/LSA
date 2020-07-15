num_chain = 3
double_z = 0
double_w = 0

out_dir = "pisa-singleW-singleZ/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p figure/",out_dir))
system(paste0("mv figure/*.pdf figure/",out_dir))

out_dir = "opusIII-singleW-singleZ/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p figure/",out_dir))
system(paste0("mv figure/*.pdf figure/",out_dir))

out_dir = "verbal-singleW-singleZ/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p figure/",out_dir))
system(paste0("mv figure/*.pdf figure/",out_dir))
