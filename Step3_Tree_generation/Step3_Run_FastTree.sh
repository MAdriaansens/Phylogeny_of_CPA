#!/bin/bash -e
#SBATCH --job-name=FT  # job name (shows up in the queue)
#SBATCH --time=168:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=4GB          # Memory in MB
#SBATCH --cpus-per-task=3
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/FT_output%A.out
#SBATCH --error=slurm_output/FT_error%A.err


module load FastTree/2.1.11-GCC-12.3.0
 export OMP_NUM_THREADS=3
FastTree -lg -cat 20 -gamma /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/CPA_April30_merged_alignedPF0999_67392.faa  > /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/CPA_28April_tree_seeds_lg_cat_gamma20.treefile
