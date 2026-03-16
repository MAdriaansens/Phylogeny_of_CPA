#!/bin/bash -e
#SBATCH --account=uc04105
#SBATCH --job-name=parse_filter_cat
#SBATCH --time=72:00:00
#SBATCH --mem=12GB
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --error=slurm_output/parse_seq_%A-%a.err
#SBATCH --output=slurm_output/parse_seq_%A-%a.out
#SBATCH --partition=large

HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign
merged=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/merged_HMMsearch_MMseq
HMMdir=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch
MMseqdir=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq
declare -a array=("03600" "06450" "03553")

module load Python/3.11.6-foss-2023a

cat HMMsearch/Archaea_PF${array[$SLURM_ARRAY_TASK_ID]}HMM_e03.faa MMseq/Archaea_Pfam${array[$SLURM_ARRAY_TASK_ID]}_seq_e03.faa > ${merged}/Archaea_merged_PF${array[$SLURM_ARRAY_TASK_ID]}.fasta

python get_uniq.py ${merged}/Archaea_merged_PF${array[$SLURM_ARRAY_TASK_ID]}.fasta 

module load HMMER/3.3.2-GCC-12.3.0

hmmalign -o ${HMMalign}/Archaea_merged_PF${array[$SLURM_ARRAY_TASK_ID]}_unique_aligned_PF${array[$SLURM_ARRAY_TASK_ID]}.sthk /nesi/nobackup/uc04105/results/HMM/PF${array[$SLURM_ARRAY_TASK_ID]}.hmm ${merged}/Archaea_merged_PF${array[$SLURM_ARRAY_TASK_ID]}_unique.fasta
