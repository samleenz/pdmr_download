#! /bin/bash
# slurm_dps.sh
# sam lee
# 2022-04-11
# -------EDIT HERE------
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lee.sa@wehi.edu.au 
#SBATCH -J pdmr_dps
#SBATCH -t 8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
# ----------------------

# do some things
# module load R
activate R41

Rscript R/download_paired_samples.R
