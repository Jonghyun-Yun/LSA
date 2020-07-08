#!/usr/bin/env bash
export STAN_NUM_THREADS=8
mkdir -p output
rm output/*
for v in {1..5}
do
./main parallel full latent no_gamma $v 1000 1000 10
done
