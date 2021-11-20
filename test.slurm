#!/bin/bash

# testing some scripts

#SBATCH --job-name=tcline_test
#SBATCH -n 4
#SBATCH -N 2
#SBATCH -p normal
#SBATCH --account=norock
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcline@usgs.gov
#SBATCH -o test.out

module load gnu8/8.3.0 R/4.0.1

prun Rscript script.R
