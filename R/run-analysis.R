num_chain = 3; double_z = 0; double_w = 0; HAS_REF = 0;
out_dir = "chessB-singleZ-singleW/"
system(paste0("rm figure/*.pdf"))
source("R/art-analysis.R")
system(paste0("mkdir -p ", out_dir, "figure/"))
system(paste0("rsync -rv figure/*.pdf ", out_dir,"figure/"))
