#!/bin/bash
#SBATCH --account       uc04105
#SBATCH --job-name      HMMprokka
#SBATCH --time          72:00:00
#SBATCH --mem           3GB
#SBATCH --cpus-per-task 8
#SBATCH --array=0-2
#SBATCH --error         pipeline_hmm/slurm_prokka_%A-%a.err
#SBATCH --output        pipeline_hmm/slurm_prokka_%A-%a.out

PDB=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta
TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMsearch_all_vsPF00999 
HMMdir=/nesi/nobackup/uc04105/results/HMM/PF00999.hmm
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain

declare -a array=('Arc_7juli_merged_MMSEQ_vsEukarya_PF00999_aligned_fulllength' 'Bac_7juli_merged_MMSEQ_vsEukarya_PF00999_aligned_fulllength' 'Self_6juli_merged_MMSEQ_vsEukarya_PF00999hmmaligned_fulllength')

#module load Python/3.11.6-foss-2023a

module load HMMER/3.3.2-GCC-12.3.0
hmmsearch --noali --cpu 8 -E 0.001 --domtblout ${HMMsearch}/${array[$SLURM_ARRAY_TASK_ID]}_fl_domtblout_PF009999.tsv ${HMMdir} ${HMMalign}/${array[$SLURM_ARRAY_TASK_ID]}.fasta
