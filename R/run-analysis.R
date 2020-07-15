num_chain = 3
double_z = 0
double_w = 0

out_dir = "pisa-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mv figure/*.pdf ", out_dir,"figure/"))

out_dir = "opusIII-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("mv figure/*.pdf figure/",out_dir))

out_dir = "verbal-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("mv figure/*.pdf figure/",out_dir))
