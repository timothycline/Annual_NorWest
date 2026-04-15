#!/bin/bash
#SBATCH --job-name=ExtractDayMet
#SBATCH --partition=priority
#SBATCH --array=1-8
#SBATCH --cpus-per-task=12
#SBATCH --account=priority-timothycline
#SBATCH --mem=32G
#SBATCH --time=0-12:0:00
#SBATCH --output=ExtractDayMet-%a-%j.out
#SBATCH --error=ExtractDayMet-%a-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothy.cline2@montana.edu

date
echo "partition = ${SLURM_JOB_PARTITION}"
echo "Hello from $(hostname -s)"

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a
module load UDUNITS/2.2.28-GCCcore-13.3.0

Rscript ProcessUnits_Extract_DayMet.R ${SLURM_ARRAY_TASK_ID}

date
