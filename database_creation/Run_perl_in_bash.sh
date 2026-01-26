#!/bin/bash
#SBATCH --job-name=perl
#SBATCH --time=48:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=9000MB          # Memory in MB
#SBATCH --partition=long
#SBATCH --cpus-per-task=5
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/perl03553_output%A.out
#SBATCH --error=slurm_output/perl03553_errror%A.err

perl download_PF03600.pl
