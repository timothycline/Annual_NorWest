#!/bin/bash

# testing some scripts

#SBATCH --job-name=SpokootSSN_Refit
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --account=norock
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcline@usgs.gov
#SBATCH -o SpokootDayMet.out

module load gnu8/8.3.0 R/1.1.1 gdal/3.1.0 proj/7.0.1

prun Rscript ProcessUnits_FitModels_DayMetRefit_R.R
