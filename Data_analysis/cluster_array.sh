#!/bin/bash
#SBATCH --job-name=Bac_cluster_Homology_MMseq_part1
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=6GB
#SBATCH --array=0-7
#SBATCH --account=uc04105 
#SBATCH --cpus-per-task=6
#SBATCH --error=slurm_output/Bac_cluster_seq_%A-%a.err
#SBATCH --output=slurm_output/Bac_cluster_seq_%A-%a.out
module load MMseqs2/15-6f452-gompi-2023a
declare -a array=(0.6 0.65 0.7 0.75 0.8 0.85 0.9 1)

MMseq=/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences
cd /nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences
mmseqs easy-cluster --threads 3 -c 0.0 --min-seq-id ${array[$SLURM_ARRAY_TASK_ID]} BAC_hmmscanned_hmmaligned.fasta ${MMseq}/MMseq/Bacteria_allhmmscanned_final_clustered_at${array[$SLURM_ARRAY_TASK_ID]}.fasta tmp_${array[$SLURM_ARRAY_TASK_ID]}
