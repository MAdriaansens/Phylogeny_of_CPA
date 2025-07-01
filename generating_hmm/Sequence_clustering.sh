#!/bin/bash
#SBATCH --job-name=Bac_cluster_Homology_MMseq_part1
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=55GB
#SBATCH --array=0-3
#SBATCH --account=uc04105 
#SBATCH --cpus-per-task=32
#SBATCH --error=slurm_output/Bac_cluster_seq_%A-%a.err
#SBATCH --output=slurm_output/Bac_cluster_seq_%A-%a.out
module load MMseqs2/15-6f452-gompi-2023a
declare -a array=("0.6" "0.7" "0.8" "0.9")

MMseq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/PF00999


mmseqs easy-cluster --threads 30 -c 0.0 --min-seq-id ${array[$SLURM_ARRAY_TASK_ID]} ${MMseq}/Archaea_Manual_seq_e05_cov30.fasta ${MMseq}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa tmp_${array[$SLURM_ARRAY_TASK_ID]}
