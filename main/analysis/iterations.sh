#!/bin/bash
L=16
T=40

array=(047:160  048:260   049:360  050:110 051:210 052:310 053:410 054:510 055:610 056:710 057:810 058:910 059:120 060:220 061:320 062:420 063:520 064:620 065:720 066:820 067:920 068:150 069:250 070:350 071:450 072:550 073:650 074:750 075:850 076:950  )

kappa=.0132
eta="-0.60"
M02="0.0000 "
mu03="0.0316"
#0.0224 0.0316 0.0387
rho="1.00766"

ScalarsXoneGauge=8
nruns=4
npergauge=$(( $ScalarsXoneGauge/$nruns   ))

correlator="DDTAU0TAU0 "



etam=`echo ${eta/\-/m}`
M02_decimal=`echo ${M02:2}`
mu03_decimal=`echo ${mu03:2}`

folder="eta_"$etam"_M02_"$M02_decimal"_mu03_"$mu03_decimal"_phi_double_smeared"


iter=" "
rm tmp
touch tmp

for i in `seq 0 29`
do

IFS=':' read -a split <<< ${array[$i]}
echo ${split[0]}  ${split[1]}



for ii in `seq 1 $nruns`
do

for is in `seq 0 $(( $npergauge-1)) `
do

#grep "$correlator" ../$folder/conf.0${split[0]}/run.N${split[1]}_$ii/bsmcontractions.0${split[0]}.0.0000000$is | head  -n $T  >> tmp
line=` grep cg_her_bi  ../$folder/conf.0${split[0]}/run.N${split[1]}_$ii/eta_"$etam"_M02_0000_mu03_"$mu03_decimal"_L16_rho_"$rho".log | tail -n 1 `

IFS=' ' read -a split_line <<< $line

iter=`echo $iter  ${split_line[6]}`

done
done

done
#echo $iter

aaa=`./mean $iter`
echo $eta $mu03 $aaa
