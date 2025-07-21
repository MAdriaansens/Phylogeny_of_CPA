#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      ARC00999
#SBATCH --time          72:00:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         slurm_output_cross/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_cross/slurm_prokka_%A-%a.out

DB=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.faa
TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv
HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/PF00999
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea
HMMdir=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Archaea
HMMscan=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea
BACHMM=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Bacteria/Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned.hmm
EUKHMM=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMM/PF00999/Manual_seq_cov30_e05_seqid0.7_genafpair_aligned.hmm
module load HMMER/3.3.2-GCC-12.3.0
module load Python/3.11.6-foss-2023a


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/BAC_HMM_e03vsArchaea.tsv ${BACHMM} ${DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/BAC_HMM_e03vsArchaea.tsv HMM ${TSV} ${HMMsearch}/BAC_HMM_e03vsArchaea.fa

hmmalign --amino --trim -o ${HMMalign}/BAC_HMM_e03vsArchaea.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/BAC_HMM_e03vsArchaea.fa

python parse_stockholm_filter.py ${HMMalign}/BAC_HMM_e03vsArchaea.sthk ${HMMalign}/BAC_HMM_e03vsArchaea.fa 258

python retrieve_fl.py BAC_HMM_e03vsArchaea.fa.fasta BAC_HMM_e03vsArchaea_fulllength.fa.fasta ${HMMalign} ${TSV}


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/EUK_HMM_e03vsArchaea.tsv ${EUKHMM} ${DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/EUK_HMM_e03vsArchaea.tsv HMM ${TSV} ${HMMsearch}/EUK_HMM_e03vsArchaea.fa

hmmalign --amino --trim -o ${HMMalign}/EUK_HMM_e03vsArchaea.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/EUK_HMM_e03vsArchaea.fa

python parse_stockholm_filter.py ${HMMalign}/EUK_HMM_e03vsArchaea.sthk ${HMMalign}/EUK_HMM_e03vsArchaea.fa 258

python retrieve_fl.py EUK_HMM_e03vsArchaea.fa.fasta EUK_HMM_e03vsArchaea_fulllength.fa.fasta ${HMMalign} ${TSV}
