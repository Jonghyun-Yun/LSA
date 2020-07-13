#!/usr/bin/env bash
export STAN_NUM_THREADS=11
mkdir -p output
rm output/*
Rscript "R/pisa-preprocess.R"
for v in {1..3}
do
Rscript "R/pisa-init.R"
./main parallel double_z full latent no_gamma $v 10000 10000 10
done
mv output pisa-new

rm output/*
mkdir -p output
Rscript "R/verbal-preprocess.R"
for v in {1..3}
do
Rscript "R/verbal-init.R"
./main parallel double_z sparse latent no_gamma $v 10000 10000 10
done
mv output verbal-new

rm output/*
mkdir -p output
Rscript "R/opusIII-preprocess.R"
for v in {1..3}
do
Rscript "R/opusIII-init.R"
./main parallel double_z sparse latent no_gamma $v 10000 10000 10
done
mv output opusIII-new
