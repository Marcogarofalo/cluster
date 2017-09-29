#!/bin/sh`ZZ
########infile
kappas=" 0.130600   0.130900 0.131000 0.131100 0.131180  0.131200 0.131220  0.131240 0.131260 0.131280  0.131300 0.131320 0.131400  0.131500 0.131600 0.131700  0.131800 0.132000  0.132400"
kappas="   0.129280  "
kappas=" 0.131850 "
#kappas="0.132283"
#renormalization scan
#kappas=" 0.131400  0.131500 0.131600 0.131700  0.131800 0.131900 0.132000  0.132400"
#kappas="0.131680    0.131700 0.131720 "
#kappas="0.130000"
#kappas=" 0.131690 0.131695  0.131700"
#kappas="0.131695"
#bastian
#lambda=0.145235

#simulation
#b=0.123
#lambda=0.145175
#b=0.093
#lambda=0.148388
#0.01031
#O(N) b0.093
#lambda=0.1485525

#beta=5.95
lambda=0.150458

L=24
T=48  #$(( $L*2 ))

replicas=" "
for rep in `seq 1 100`
do
for kappa in $kappas
do
infile="infile_L"$L"T"$T"_k"$kappa"_lambda"$lambda"_rep"$rep".in"
output="L"$L"T"$T"_k"$kappa"_lambda"$lambda".out"

echo L = $T $L $L $L > $infile
#NOcontinumm
echo formulation = continuum >> $infile
echo kappa = $kappa >> $infile
echo lambda = $lambda >>$infile

echo metropolis_local_hits = 10 >> $infile
echo metropolis_global_hits = 1 >> $infile
echo metropolis_delta = 4.7  >> $infile

echo cluster_hits = 32 >>$infile
echo cluster_min_size = 0.1  >>$infile

echo seed = $RANDOM >> $infile
echo restart = 0 >> $infile
echo replica = $rep >>$infile
echo start_measure = 100 >> $infile
echo total_measure = 20101 >> $infile
echo measure_every_X_updates = 10 >> $infile
echo save_config = no >> $infile
echo save_config_rotated = yes >> $infile
echo save_config_every_X_updates = 10 >> $infile

echo outpath ="./beta5.95" >> $infile

sed -i "8s/.*/#$ -o "$output" /" qsub_cluster.sh
# Run the program


#./run_cluster -i $infile  &
qsub qsub_cluster.sh   $infile $output
done
done
