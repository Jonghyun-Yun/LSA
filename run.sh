#!/usr/bin/env bash
out_dir="chessB_pn_ncut2_zero_beta_noinfo/"
prename="R/chessB_noinfo-preprocess.R"
initname="R/chessB_pn-init.R"
n_chain=2

echo "================================"
echo "Output dir:" $out_dir
echo "preprocessing:" $prename
echo "initializing:" $initname
echo "n_chain:" $n_chain
echo "================================"

export STAN_NUM_THREADS=2
mkdir -p output
rm output/*
rm input/*

Rscript $prename
cp input/{mvar,mlen}.csv output/
cp run.sh output/
cp $prename output/
cp $initname output/

for ((v = 1; v <= $n_chain; v++))
do
    Rscript $initname
    ./main initialize parallel single_w single_z full latent no_gamma true false no_ars nonzero_theta zero_beta lambda_free $v 1000 1000 10
done
## explain comandline arguments:
 # continue initialize -> start new chains?
 # parallel serial -> parallel computation?
 # single_w double_w
 # single_z double_z
 # full sparse -> missing or not
 # latent no_latent -> update latent space?
 # gamma no_gamma -> update gamma?
 # true false -> can I play with the gamma sign?
 # correct incorrect false -> gamma for what process?
 # do_ars no_ars -> ARS for gamma
 # zero_theta nonzero_theta -> theta.k.0 are zero
 # double_beta single_beta

mkdir -p $out_dir
mv output/* $out_dir
Rscript R/run-analysis.R $out_dir $n_chain
echo "Outputs are moved to" $out_dir"."
echo "================================"

# sh run_run.sh
# sh run_run_run.sh
# sh run_run_run_run.sh
# sh run_run_run_run_run.sh
# sh run_run_run_run_run_run.sh
