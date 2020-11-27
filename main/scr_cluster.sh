#!/bin/sh

#SBATCH --partition=skl_usr_prod
#SBATCH --time=20:00:00
#SBATCH --account=INF20_lqcd123


/marconi/home/userexternal/mgarofal/cluster_eddie/main/run_cluster -i infile.in  
