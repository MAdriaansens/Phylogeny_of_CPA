#!/bin/bash
#SBATCH --job-name=00999_Homology_Euk
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=18500MB          # Memory in MB
#SBATCH --partition=milan
#SBATCH --cpus-per-task=11
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/00999_Homology_Euk_output%A.out
#SBATCH --error=slurm_output/00999_Homology_Euk_error%A.err

PWD=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMsearch
HMMdir=/nesi/nobackup/uc04105/results/HMM
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign

Euk_TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
Euk_DB=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta


module load Python/3.11.6-foss-2023a
module load HMMER/3.3.2-GCC-12.3.0

#07
hmmsearch --noali --cpu 10 -E 0.0000001 --tblout ${HMMsearch}/PF00999_07_vsEuk.tsv ${HMMdir}/PF00999.hmmm ${Euk_DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/PF00999_07_vsEuk.tsv HMM ${Euk_TSV} ${HMMsearch}/PF000999_07_vsEuk.fasta EukM6

#hmmalign --amino --trim -o /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF00999_07_vsEuk_AlignedPF00999.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm /nesi/nobackup/uc04105//new_databases_May/Euk_database_May/results/HMMsearch/PF00999_07_vsEuk.fasta

python parse_stockholm_filter.py ${HMMalign}/PF00999_07_vsEuk_AlignedPF00999.sthk ${HMMalign}/PF00999_07__vsEuk_AlignedPF00999 257
