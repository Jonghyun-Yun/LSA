#!/usr/bin/env bash
out_dir="chessB_no_latent/"
n_chain=1

export STAN_NUM_THREADS=2
mkdir -p output
rm output/*
rm input/*

Rscript "R/chessB-preprocess.R"
cp input/{mvar,mlen}.csv output/
cp run.sh output/

for ((v = 1; v <= $n_chain; v++))
do
    Rscript "R/chessB-init.R"
    ./main initialize parallel single_w single_z full no_latent no_gamma true false no_ars $v 10000 10000 10
done

 # continue initialize
 # parallel serial
 # single_w double_w
 # single_z double_z
 # full sparse
 # latent no_latent
 # gamma no_gamma
 # true false
 # correct incorrect false
 # do_ars no_ars

mkdir -p $out_dir
cp R/chessB-init.R $out_dir
mv output/* $out_dir
Rscript R/run-analysis.R $out_dir $n_chain
