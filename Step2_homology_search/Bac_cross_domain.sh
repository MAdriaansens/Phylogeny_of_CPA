#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      BAC00999
#SBATCH --time          72:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 16
#SBATCH --array         0-61
#SBATCH --error         slurm_output_cross/BACslurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_cross/BACslurm_prokka_%A-%a.out

DB=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/DB/fasta
TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/DB/tsv

HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria
HMMdir=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Bacteria

ARCHMM=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Archaea/Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned.hmm
EUKHMM=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMM/PF00999/Manual_seq_cov30_e05_seqid0.7_genafpair_aligned.hmm

module load HMMER/3.3.2-GCC-12.3.0
module load Python/3.11.6-foss-2023a


declare -a array=($(seq 0 61))

hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.tsv ${ARCHMM} ${DB}/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta

python getting_fasta_from_hit_extra.py ${HMMsearch}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.tsv HMM ${TSV}/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv ${HMMsearch}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa

#mkdir -p ${HMMalign}/HMMsearch

hmmalign --amino --trim -o ${HMMalign}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa

python parse_stockholm_filter.py ${HMMalign}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.sthk ${HMMalign}/ARC_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa 258


hmmsearch --noali --cpu 15 -E 0.001 --tblout ${HMMsearch}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.tsv ${EUKHMM} ${DB}/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta

python getting_fasta_from_hit_extra.py ${HMMsearch}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.tsv HMM ${TSV}/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv ${HMMsearch}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa

#mkdir -p ${HMMalign}/HMMsearch

hmmalign --amino --trim -o ${HMMalign}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm ${HMMsearch}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa

python parse_stockholm_filter.py ${HMMalign}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.sthk ${HMMalign}/EUK_HMM_e03vsBacteria${array[$SLURM_ARRAY_TASK_ID]}.fa 258
