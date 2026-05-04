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

mmseqs easy-cluster --threads 12 -c 0.0 --min-seq-id 0.7 /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Archaea_April_hmmaligned_e03_mmseq.fasta.fa /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Archaea_April_hmmaligned_e03_mmseq_treeinput_clusterd_at_${array[$SLURM_ARRAY_TASK_ID]}.faa  tmp_${array[$SLURM_ARRAY_TASK_ID]}

mmseqs easy-cluster --threads 12 -c 0.0 --min-seq-id 0.7 /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Eukarya_April_hmmaligned_e03_mmseq.fasta.fa /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Eukarya_April_hmmaligned_e03_mmseq_treeinput_clusterd_at_${array[$SLURM_ARRAY_TASK_ID]}.faa  tmp_${array[$SLURM_ARRAY_TASK_ID]}

mmseqs easy-cluster --threads 12 -c 0.0 --min-seq-id 0.7 /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Bacteria_April_hmmaligned_e03_mmseq.fasta.fa /home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/final_april28/tree_input/Bacteria_April_hmmaligned_e03_mmseq_treeinput_clusterd_at_${array[$SLURM_ARRAY_TASK_ID]}.faa  tmp_${array[$SLURM_ARRAY_TASK_ID]}
