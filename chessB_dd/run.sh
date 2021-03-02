#!/usr/bin/env bash
export STAN_NUM_THREADS=2
mkdir -p output
rm output/*
Rscript "R/chessB-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for v in {1..3}
do
    Rscript "R/chessB-init.R"
    ./main initialize parallel double_w double_z full latent no_gamma true $v 10000 10000 10
done
mv output chessB_dd
Rscript R/run-analysis.R chessB_dd/
