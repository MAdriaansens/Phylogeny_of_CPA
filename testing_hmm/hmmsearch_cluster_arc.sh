#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      00999INH
#SBATCH --time          72:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 16
#SBATCH --error         slurm_output/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output/slurm_prokka_%A-%a.out

DB=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.faa
TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv
HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/PF00999
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea
HMMdir=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Archaea



module load HMMER/3.3.2-GCC-12.3.0
module load Python/3.11.6-foss-2023a


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/$1_HMM_e03vsArchaea.tsv ${HMMdir}/${1}.hmm ${DB}

python getting_fasta_from_hit_extra.py ${HMMsearch}/$1_HMM_e03vsArchaea.tsv HMM ${TSV} ${HMMsearch}/$1_HMM_e03vsArchaea.fa

#mkdir -p ${HMMalign}/HMMsearch

hmmalign --amino --trim -o ${HMMalign}/$1_HMM_e03vsArchaea.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/$1_HMM_e03vsArchaea.fa

python parse_stockholm_filter.py ${HMMalign}/$1_HMM_e03vsArchaea.sthk ${HMMalign}/$1_HMM_e03vsArchaea.fa 258
