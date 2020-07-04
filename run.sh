#!/usr/bin/env bash
export STAN_NUM_THREADS=11
mkdir -p output
rm output/*
for v in {1..5}
do
./main parallel $v 5000 5000 10
done
