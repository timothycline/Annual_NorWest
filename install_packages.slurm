#!/bin/bash
#SBATCH --job-name=InstallPackages
#SBATCH --partition=priority
#SBATCH --cpus-per-task=1
#SBATCH --account=priority-timothycline
#SBATCH --mem=8G
#SBATCH --time=0-1:0:00
#SBATCH --output=InstallPackages-%j.out
#SBATCH --error=InstallPackages-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothy.cline2@montana.edu

date
echo "partition = ${SLURM_JOB_PARTITION}"
echo "Hello from $(hostname -s)"

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

Rscript install_packages.R

date
