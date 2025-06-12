#!/bin/bash -e

#SBATCH --account       uc04105
#SBATCH --job-name      03600prokka
#SBATCH --time          72:00:00
#SBATCH --mem           50GB
#SBATCH --array         0-2
#SBATCH --cpus-per-task 14
#SBATCH --error         slurm_output_03600/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_03600/slurm_prokka_%A-%a.out
PFAM_fasta=/nesi/nobackup/uc04105/database/Pfam_download
DB=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta
TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
MMSEQS=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/MMseq
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign
module load MMseqs2/15-6f452-gompi-2023a

declare -a array=("PF06450" "PF03553" "PF03600")

module load HMMER/3.3.2-GCC-12.3.0

mmseqs easy-search -e 1.00E-03 --threads 12 ${PFAM_fasta}/${array[$SLURM_ARRAY_TASK_ID]}_sequences_size_excluded_cd09.fasta ${DB} ${MMSEQS}/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}MMSEQ_vsEukarya.tsv ${MMSEQS}/tmp_${array[$SLURM_ARRAY_TASK_ID]}

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${MMSEQS}/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}MMSEQ_vsEukarya.tsv MMSEQ ${TSV} ${MMSEQS}/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}MMSEQ_vsEukarya.fa

hmmalign --amino --trim -o ${HMMalign}/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}MMseq_vsEukarya_aligned_${array[$SLURM_ARRAY_TASK_ID]}.tsv /nesi/nobackup/uc04105/results/HMM/${array[$SLURM_ARRAY_TASK_ID]}.hmm ${MMSEQS}/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}MMSEQ_vsEukarya.fa
~            
