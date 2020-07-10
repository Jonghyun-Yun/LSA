#!/usr/bin/env bash
export STAN_NUM_THREADS=11
mkdir -p output
rm output/*
for v in {1..5}
do
./main parallel single_z sparse latent gamma $v 10000 5000 5
done
