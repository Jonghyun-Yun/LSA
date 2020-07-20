#!/usr/bin/env bash
export STAN_NUM_THREADS=10

mkdir -p output
rm -r output/*

rsync -rv pisa-KR-sci2018-singleZ-singleW/ output

Rscript "R/pisa-KR-sci2018-preprocess.R"
cp input/{mvar,mlen}.csv output/

touch output/run_hist.sh
cat run.sh >> output/run_hist.sh

for v in {1..3}
do
    ./main continue parallel single_w single_z full latent no_gamma $v 10000 10 10
done
rsync -rv output/ pisa-KR-sci2018-singleZ-singleW

mkdir -p output
rm -r output/*

rsync -rv pisa-KR-singleZ-singleW/ output

Rscript "R/pisa-KR-preprocess.R"
cp input/{mvar,mlen}.csv output/

touch output/run_hist.sh
cat run.sh >> output/run_hist.sh

for v in {1..3}
do
    ./main continue parallel single_w single_z full latent no_gamma $v 10000 10 10
done
rsync -rv output/ pisa-KR-singleZ-singleW
