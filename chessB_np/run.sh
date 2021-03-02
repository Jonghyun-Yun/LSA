#!/usr/bin/env bash
export STAN_NUM_THREADS=3
mkdir -p output
rm output/*
Rscript "R/chessB-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for v in {1..3}
do
    Rscript "R/chessB_np-init.R"
    ./main initialize parallel single_w single_z full latent no_gamma true $v 10000 10000 10
done
mv output chessB_np
Rscript R/run-analysis.R chessB_np/

# rm output/*
# Rscript "R/chessB-preprocess.R"
# cp input/{mvar,mlen}.csv output/
# cp run.sh output/

# for v in {1..3}
# do
#     Rscript "R/chessB_no_latent-init.R"
#     ./main initialize parallel single_w single_z full latent no_gamma true $v 10000 10000 10
# done
# mv output chessB_no_latent
