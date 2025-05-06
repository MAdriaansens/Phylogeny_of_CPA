#!/bin/bash
#SBATCH --job-name=00999_Homology_Euk_direct
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=18500MB          # Memory in MB
#SBATCH --partition=milan
#SBATCH --cpus-per-task=11
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/00999_Homology_Euk_direct_output%A.out
#SBATCH --error=slurm_output/00999_Homology_Euk_direct_error%A.err

module load Python/3.11.6-foss-2023a
module load HMMER/3.3.2-GCC-12.3.0

HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign
MMseqs=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/MMseq
HMMdir=/nesi/nobackup/uc04105/results/HMM

#this code takes the mmseq output and directly runs a hmmalign and filter on it. 

hmmalign --amino --trim -o ${HMMalign}/Pfam00999_seq_cov50_e03_alignedPF00999.sthk ${HMMdir}/PF00999.hmm  
${MMseqs}/Pfam00999_seq_cov50_e03.faa

python parse_stockholm_filter.py ${HMMalign}/Pfam00999_seq_cov50_e03_alignedPF00999.sthk ${HMMalign}/Pfamm
00999_seq_cov50_e03_alignedPF00999 257
