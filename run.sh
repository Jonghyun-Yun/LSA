#!/usr/bin/env bash
export STAN_NUM_THREADS=11
mkdir -p output
rm output/*
Rscript "R/verbal-preprocess.R"
for v in {1..5}
do
Rscript "R/verbal-init.R"
./main parallel single_z sparse latent no_gamma $v 10000 10000 10
done
Rscript "R/art-analysis.R"
mv output verbal-nothetabeta
mkdir -p figure/verbal-nothetabeta
mv figure/*.pdf figure/verbal-nothetabeta/

rm output/*
Rscript "R/opusIII-preprocess.R"
for v in {1..5}
do
Rscript "R/opusIII-init.R"
./main parallel single_z sparse latent no_gamma $v 10000 10000 10
done
Rscript "R/art-analysis.R"
