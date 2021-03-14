#!/usr/bin/env bash
export STAN_NUM_THREADS=2
mkdir -p output
rm output/*
Rscript "R/chessB-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for v in {1..3}
do
    Rscript "R/chessB_pf-init.R"
    ./main initialize parallel single_w single_z full latent gamma true correct $v 10000 10000 10
done

cp R/chessB_pf-init.R chessB_pf/
mv output chessB_pf
Rscript R/run-analysis.R chessB_pf/
