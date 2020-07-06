#!/usr/bin/env bash
export STAN_NUM_THREADS=10
mkdir -p output
rm output/*
for v in {1..5}
do
./main parallel latent no_gamma $v 2500 2500 10
done
