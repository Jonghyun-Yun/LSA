num_chain = 3; double_z = 0; double_w = 0; HAS_REF = 0;
out_dir = "pisa-KR-sci2018-singleZ-singleW/"

out_dir = "pisa-KR-sci2018-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("rsync -rv figure/*.pdf ", out_dir,"figure/"))

ref_dir = "pisa-singleZ-singleW/"
HAS_REF = 1

out_dir = "pisa-KR-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("rsync -rv figure/*.pdf ", out_dir,"figure/"))

HAS_REF = 0

out_dir = "marketing-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("mv figure/*.pdf ", out_dir,"figure/"))
