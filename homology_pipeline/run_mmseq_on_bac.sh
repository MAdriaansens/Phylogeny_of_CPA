#!/bin/bash -e

#SBATCH --account       uc04105
#SBATCH --job-name      00999prokka
#SBATCH --time          72:00:00
#SBATCH --mem           50GB
#SBATCH --array         0-61
#SBATCH --cpus-per-task 14
#SBATCH --error         slurm_output_00999/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_00999/slurm_prokka_%A-%a.out

PFAM_fasta=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria
DB=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/DB/fasta
TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/DB/tsv
MMSEQS=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea


module load MMseqs2/15-6f452-gompi-2023a

declare -a array=($(seq 1 61))


#mmseqs easy-search -e 1.00E-03 --threads 12 /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF06450/PF06450_merged_Bac_fl.fasta ${DB}/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF06450/PF06450Bac_own_merged_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}.tsv  tmp_${array[$SLURM_ARRAY_TASK_ID]}

module load Python/3.11.6-foss-2023a
#echo 's_May/GTDB_226/results/MMseq/Bacteria/PF06450/PF06450Bac_own_merged_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}.tsv'
python getting_fasta_from_hit_extra.py /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF06450/PF06450Bac_own_merged_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}.tsv MMSEQ ${TSV}/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF06450/PF06450Bac_own_merged_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
