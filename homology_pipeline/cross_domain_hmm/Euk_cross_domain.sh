#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      EUK00999INH
#SBATCH --time          72:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 16
#SBATCH --error         slurm_output/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output/slurm_prokka_%A-%a.out

DB=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta
TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF00999
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMsearch/PF00999
HMMdir=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMM/PF00999
ARCHMM=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Archaea/Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned.hmm
BACHMM=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Bacteria/Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned.hmm


module load HMMER/3.3.2-GCC-12.3.0
module load Python/3.11.6-foss-2023a


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/BAC_HMM_e03vsEukarya.tsv ${BACHMM} ${DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/BAC_HMM_e03vsEukarya.tsv HMM ${TSV} ${HMMsearch}/BAC_HMM_e03vsEukarya.fa

#mkdir -p ${HMMalign}/HMMsearch

hmmalign --amino --trim -o ${HMMalign}/BAC_HMM_e03vsEukarya.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/BAC_HMM_e03vsEukarya.fa

python parse_stockholm_filter.py ${HMMalign}/BAC_HMM_e03vsEukarya.sthk ${HMMalign}/BAC_HMM_e03vsEukarya.fa 258


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/ARC_HMM_e03vsEukarya.tsv ${ARCHMM} ${DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/ARC_HMM_e03vsEukarya.tsv HMM ${TSV} ${HMMsearch}/ARC_HMM_e03vsEukarya.fa

#mkdir -p ${HMMalign}/HMMsearch

hmmalign --amino --trim -o ${HMMalign}/ARC_HMM_e03vsEukarya.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/ARC_HMM_e03vsEukarya.fa

python parse_stockholm_filter.py ${HMMalign}/ARC_HMM_e03vsEukarya.sthk ${HMMalign}/ARC_HMM_e03vsEukarya.fa 258
