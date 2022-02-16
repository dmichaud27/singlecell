#!/bin/bash

#SBATCH -p general
#SBATCH -t 20:00:00
#SBATCH --mem=80g
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -o R.slurm.out.%J

module add r/4.1.0
Rscript Cluster_ID2.R