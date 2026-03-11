#!/bin/bash
#SBATCH --job-name=Run_pythonred
#SBATCH --time=336:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=5000MB          # Memory in MB
#SBATCH --partition=long
#SBATCH --cpus-per-task=5
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/Run_pythonred_output%A.out
#SBATCH --error=slurm_output/Run_pythonred_error%A.err

module purge && module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#module load Python/3.11.6-foss-2023a
#conda create -n redvals_env python=3.12 biopython pandas tqdm 

conda activate redvals_env

python run_red.py
