#!/bin/bash
L=16
T=40

array=(047:160  048:260   049:360  050:110 051:210 052:310 053:410 054:510 055:610 056:710 057:810 058:910 059:120 060:220 061:320 062:420 063:520 064:620 065:720 066:820 067:920 068:150 069:250 070:350 071:450 072:550 073:650 074:750 075:850 076:950  )
#array=( 055:610 )

kappa=.0132
eta="-0.9"
M02="0.0000 "
mu03="0.0224"
#0.0224 0.0316 0.0387
rho="2.015326"

ScalarsXoneGauge=8
nruns=4
npergauge=$(( $ScalarsXoneGauge/$nruns   ))

correlator="P1P1TRIVIAL"



etam=`echo ${eta/\-/m}`
M02_decimal=`echo ${M02:2}`
mu03_decimal=`echo ${mu03:2}`

folder="eta_"$etam"_M02_"$M02_decimal"_mu03_"$mu03_decimal"_rho"$rho"_phi_double_smeared"
imax=${#array[@]} 
imax=$((imax-1))

g++ analysis.cpp -o analysis
mkdir DATA_tmp

for  f in $folder
do
	confs=` ls  ../$folder | grep conf.0`

	for conf in $confs
	do
		runs=` ls  ../"$folder"/"$conf"/ | grep run `
		for run in $runs
		do
		sXgs=` ls  ../"$folder"/"$conf"/"$run" | grep bsmcontractions`
               
			for sXg in $sXgs
			do
			infile="../$f/$conf/$run/$sXg"
			grep "$correlator"  $infile >> DATA_tmp/"$correlator"_"$conf"
                	full=`grep "$correlator"  $infile | wc -l  `
			done	
		done
	
	done
done

file="DATA_tmp/"$correlator"_eta"$etam"_M02_"$M02_decimal"_mu03_"$mu03_decimal"_rho"$rho"_phi_double_smeared1.tmp"
out="DATA/"$correlator"_eta"$etam"_M02_"$M02_decimal"_mu03_"$mu03_decimal"_rho"$rho"_phi_double_smeared1.dat"
if [ -f $file ]
then
rm $file
fi

confs=`ls DATA_tmp/ | grep $correlator `
gauges=0
for conf in $confs
do
N=`cat DATA_tmp/$conf | wc -l`
N=$(( N/$T ))
echo $conf
./analysis -i DATA_tmp/$conf -o $file -L $L -T $T -N $N

gauges=$((gauges+1))
done

./effective_mass_gauge_log  -i $file -o $out -L $L -T $T -N $gauges 

rm -r DATA_tmp 


