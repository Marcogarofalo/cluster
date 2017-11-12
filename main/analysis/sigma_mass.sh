#!/bin/bash

bash compile.sh

kappas=" 0.131870 "
kappas="   0.132000  "

#bastian
#clambda=0.145235

#BSMM
clambda=0.147925
clambda=0.150458
clambda=0.150450
#clambda=0.150375
#clambda=0.145175
#b0.093
#clambda=0.148388
#clambda=0.1485525
#beta=5.95
clambda=0.15045
clambda=0.1504375
clambda=0.150425
clambda=0.15070
clambda=0.147925

L=16
T=40  #$(( $L*2 ))
replicas=150


for kappa in $kappas
do

lambda=`./lambda  $kappa $clambda `
conf=0
set=`seq 1 $replicas`
#set=`echo $set | sed "s/148//"`
#set=`echo $set | sed "s/149//"`
#set=`echo $set | sed "s/150//"`
for rep in $set
do

infile="../beta5.95/G2t.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
infile1="../beta5.95/mag.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
outfile="sigma_mass.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
confs=` wc -l $infile | awk '{print $1;}'  `
#confs=$((  confs / T  ))
echo $((confs/$T))   rep $rep
if [ $confs -gt 0 ]
then
tail -n $(( $confs - (10*$T) ))  $infile >> tmp
tail -n $(( $confs/$T - 10 ))  $infile1 >> tmp1
conf=$(( conf+ confs/$T - 10  ))
fi
done

#infile="../cluster/main/data/T"$T".X"$L".Y"$L".Z"$L".kap"$kappa".lam"$lambda".rep_1.rot_no.conf"
#infile="../cluster/main/data_b093/L"$L"T"$T"_k"$kappa"_lambda"$clambda".out"
#infile="../cluster/main/data_b093/mag.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_1.dat"

#cat $infile >> ../cluster/main/data_b093/mag.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_2.dat
#echo $infile

blocking=1000
blocks=`cat tmp | wc -l`
blocks=$(( blocks/($T*$blocking)  ))
echo blocking $blocks
for i in `seq 1 $blocks `
do
echo block number $i
head -n  $(( $i*$T*$blocking )) tmp  | tail -n $(($T*$blocking)) > tmp2
head -n $((i*$blocking)) tmp1  | tail -n $blocking > tmp3
./analysis -i tmp2 -o tmp4 -L $L -T $T -N $blocking
done



res="$kappa"
#res1=`./magnetization -L $L -T $T -N $conf -e 1 -i tmp -o bo`
#res1=`./vev -L $L -T $T -N $conf -e 1 -i tmp -o bo`
echo $conf  $confs   $lambda
echo $kappa  
./effective_mass_gauge_log -i tmp4 -o $outfile -L $L -T $T -N $((conf/$blocking)) -i tmp5
./jacknife -i tmp4 -o jacknife_$outfile -L $L -T $T -N $((conf/$blocking)) -i tmp5

./Z_phi -i tmp5 -o bo -L $L -T $T -N $((conf/blocking)) -e 1 > wolff_$outfile


rm tmp2 
rm tmp3
rm tmp4
rm tmp5
rm tmp1
rm tmp
#echo "$res    $res1" >> susceptivity.dat

done
