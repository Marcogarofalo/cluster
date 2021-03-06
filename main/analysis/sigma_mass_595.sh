kappas="0.130000 0.130400   0.130800 0.131000 0.131200 0.131400 0.131600 0.131800  0.132000  0.132200 0.132400 0.132600 0.133000 "

#kappas="0.130000 0.130400 0.130800 0.131000 0.131200 0.131400 0.131600 0.131800 0.132000 0.132200 0.132400 0.132600  0.133000 0.133500 0.134000"
#kappas="0.130000 0.130600  0.130800 0.130900 0.131000 0.131100  0.131200 0.131300 0.131400  0.131500 0.131600 0.131700  0.131800 0.132000 0.132200 0.132400"
kappas=" 0.130600   0.130900 0.131000 0.131100 0.131180  0.131200 0.131220  0.131240 0.131260 0.131280  0.131300 0.131320 0.131400  0.131500 0.131600 0.131700  0.131800 0.132000  0.132400"


#renormalization scan
kappas=" 0.131400  0.131500 0.131600 0.131700  0.131800 0.131900 0.132000  0.132400"
kappas=" 0.131720"
kappas=" 0.131900 "

#kappas="0.130000 0.132000"
#kappas=" 0.131690 0.131695  0.131700"
#bastian
#clambda=0.145235

#BSMM
clambda=0.147925
clambda=0.150458
#clambda=0.145175
#b0.093
#clambda=0.148388
#clambda=0.1485525
L=24
T=48  #$(( $L*2 ))
replicas=1


for kappa in $kappas
do

lambda=`./lambda  $kappa $clambda `
conf=0
for rep in `seq 1 $replicas `
do
infile="../beta5.95/G2t.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
infile1="../beta5.95/mag.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
outfile="sigma_mass.T"$T"X"$L"Y"$L"Z"$L"kap"$kappa"lam"$lambda".rep_"$rep".dat"
confs=` wc -l $infile | awk '{print $1;}'  `
#confs=$((  confs / T  ))
echo $((confs/$T))
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

blocking=1
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

./Z_phi -i tmp5 -o bo -L $L -T $T -N $((conf/blocking)) -e 1 > wolff_$outfile
./jacknife -i tmp4 -o bo -L $L -T $T -N $((conf/$blocking)) -i tmp5


rm tmp2 
rm tmp3
rm tmp4
rm tmp5
rm tmp1
rm tmp
#echo "$res    $res1" >> susceptivity.dat

done
