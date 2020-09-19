#!/usr/bin/env bash
export STAN_NUM_THREADS=6
mkdir -p output
# rm output/*
# Rscript "R/chessB-preprocess.R"
# cp input/{mvar,mlen}.csv output/
# cp run.sh output/

# for v in {1..3}
# do
#     Rscript "R/chessB_pn-init.R"
#     ./main initialize parallel single_w single_z full latent no_gamma true $v 10000 10000 10
# done
# mv output chessB_pn
# Rscript R/run-analysis.R chessB_pn/

rm output/*
Rscript "R/chessB-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for v in {1..3}
do
    Rscript "R/chessB_nn-init.R"
    ./main initialize parallel single_w single_z full latent no_gamma true $v 10000 10000 10
done
mv output chessB_nn
Rscript R/run-analysis.R chessB_nn/
