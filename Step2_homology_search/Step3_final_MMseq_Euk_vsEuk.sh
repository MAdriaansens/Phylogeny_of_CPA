#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      EukA
#SBATCH --time          24:00:00
#SBATCH --mem           15GB
#SBATCH --cpus-per-task 10
#SBATCH --error         slurm_outputE/slurm_EukA_%A-%a.err
#SBATCH --output        slurm_outputE/slurm_EukA_%A-%a.out
#SBATCH --array         0-15

declare -a array=($(seq 0 15))

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM/PF00999.hmm
Euk_TSV=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/Euk_DB/tsv/Euk_db_May_protein_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Euk_db=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/Euk_DB/fasta/Euk_db_May_protein_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
HMMsearch=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/HMMsearch/PF00999
MMseqs=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/MMseq/PF00999/part3
Seq=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/MMseq/PF00999/PF00999_cross_vsEukarya_dupes_removed_sequences.fasta

#module load MMseqs2/15-6f452-gompi-2023a

#mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Euk_db} ${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp


if [ -f "${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta" ]; then
    echo "${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta"
else
    echo "${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta doesnt exist"
    module load MMseqs2/15-6f452-gompi-2023a
    mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 10 ${Seq}  ${Euk_db} ${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp
    module load Python/3.11.6-foss-2023a
    python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/Eukfoundseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta
    rm -r ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp
fi
#


