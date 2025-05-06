#!/bin/bash
#SBATCH --job-name=00999_Homology_MMseq_part1
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=40900MB          # Memory in MB
#SBATCH --partition=milan
#SBATCH --cpus-per-task=25
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/00999_Homology_MMseq_Euk_part1_output%A.out
#SBATCH --error=slurm_output/00999_Homology_MMseq_part1_Euk_error%A.err

module load MMseqs2/12-113e3-gimkl-2020a
module load Python/3.11.6-foss-2023a

PWD=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMsearch
HMMdir=/nesi/nobackup/uc04105/results/HMM
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign

Euk_TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
Euk_DB=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May.fasta
MMseqs=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/MMseq
PFAM_seq=/nesi/nobackup/uc04105/database/Pfam_download/Pfam00999_seq/Pfam00999_sequences_size_excluded_cdd
09.faa

#run HMMsearch direct
sbatch Run_HMMsearch_direct.sh

#run MMseq and retrieve sequences

#Pfam sequences

mmseqs easy-search -e 1.00E-03 -c 0.5 --threads 22 ${PFAM_seq} ${Euk_DB} ${MMseqs}/Pfam00999_seq_cov50_e00
3.tsv tmp

python getting_fasta_from_hit_extra.py ${MMseqs}/Pfam00999_seq_cov50_e03.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PP
fam00999_seq_cov50_e03.faa  EukM6

#directly run HMMalign on the retrieved MMseq sequences
sbatch hmmalign_direct.sh


#continue with creating an inhouse HMM
mmseqs easy-cluster --threads 22 -c 0.0 --min-aln-len 0.8 ${MMseqs}/Pfam00999_seq_cov50_e03.faa ${MMseqs}}
/Pfam00999_seq_cov50_e03_seqid08cluster.faa  tmp

mafft --auto --maxiterate 1000 --anysymbol ${MMseqs}/Pfam00999_seq_cov50_e05_seqid09cluster.faa > ${MMseqq
s}/Pfam00999_seq_cov50_e05_seqid09cluster_aligned.faa
