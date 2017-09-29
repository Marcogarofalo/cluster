#!/bin/sh
# Grid Engine options

#$ -N cluster

#$ -cwd

#$ -o L24T48_k0.131850_lambda0.150458.out 

#$ -e cluster_error

#$  -l h_rt=48:00:00

### #$ -l h_vmem=0G

#### #$ -pe sharedmem 1


./run_cluster -i $1  
