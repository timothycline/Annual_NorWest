#!/bin/bash

# testing some scripts

#SBATCH --job-name=UpMissMariasSSN_Refit
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --account=norock
#SBATCH --time=2-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcline@usgs.gov
#SBATCH -o UpMissMariasSSN_Refit.out

module load gnu8/8.3.0 R/4.1.1 gdal/3.1.0 proj/7.0.1

prun Rscript ProcessUnits_FitModels_DayMetRefit_UpMissMarias.R
