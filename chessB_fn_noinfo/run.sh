#!/usr/bin/env bash
out_dir="chessB_fn_noinfo/"

export STAN_NUM_THREADS=2
mkdir -p output
rm output/*
Rscript "R/chessB_fn-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for v in {1..3}
do
    Rscript "R/chessB_fn-init.R"
    ./main initialize parallel single_w single_z full latent gamma true incorrect $v 10000 10000 10
done

mkdir -p $out_dir
cp R/chessB_nf-init.R $out_dir
mv output/* $out_dir
Rscript R/run-analysis.R $out_dir
