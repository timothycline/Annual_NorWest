#!/bin/bash

# testing some scripts

#SBATCH --job-name=SpokootSSN_Refit
#SBATCH -n 7
#SBATCH -N 7
#SBATCH -p normal
#SBATCH --account=norock
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcline@usgs.gov
#SBATCH -o DayMetRefitParallel.out

module load gnu8/8.3.0 R/1.1.1 gdal/3.1.0 proj/7.0.1

prun Rscript ProcessUnits_FitModels_DayMetRefit_Parallel.R
