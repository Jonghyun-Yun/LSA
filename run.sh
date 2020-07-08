#!/usr/bin/env bash
export STAN_NUM_THREADS=8
mkdir -p output
rm output/*
for v in {1..5}
do
./main parallel sparse no_latent no_gamma $v 5000 5000 10
done
