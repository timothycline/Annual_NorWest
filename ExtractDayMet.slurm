#!/bin/bash

#SBATCH --job-name=ExtractDayMet
#SBATCH --array=1-8
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothy.cline2@montana.edu
#SBATCH -o ExtractDayMet_%a.out

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

srun Rscript ProcessUnits_Extract_DayMet.R ${SLURM_ARRAY_TASK_ID}
