
fs=`ls scalars`
cat > eo2lex.in << EOL
L = 40 20 20 20
formulation = no
kappa = 0.132
lambda = 0.010310
metropolis_local_hits = 10
metropolis_global_hits = 1
metropolis_delta = 4.7
cluster_hits = 32
cluster_min_size = 0.1
seed = 15302
restart = 0
replica = 1
start_measure = 10000
total_measure = 10001
measure_every_X_updates = 10
save_config = yes
save_config_rotated = yes
save_config_every_X_updates = 10
outpath =./
EOL

#for f in $fs
#do

#echo ""
#echo convertin to lex scalar $f 

#cp scalars/$f in
./reorder_scalars/main/reorder_scalar  -i eo2lex.in
#mv out scalars_lex/$f

#rm in


#done


