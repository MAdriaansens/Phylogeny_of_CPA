#!/bin/bash
#SBATCH --job-name=06450_Homology_Euk_direct
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=18500MB          # Memory in MB
#SBATCH --partition=milan
#SBATCH --cpus-per-task=11
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/06450_Homology_Euk_direct_output%A.out
#SBATCH --error=slurm_output/06450_Homology_Euk_direct_error%A.err

module load Python/3.11.6-foss-2023a
module load HMMER/3.3.2-GCC-12.3.0

HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign
MMseqs=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/MMseq/Iterative_MMseqs
HMMdir=/nesi/nobackup/uc04105/results/HMM

#this code takes the mmseq output of Step 4 and directly runs a hmmalign and filter on it. 
#this is just eukarya, but can be (and has been) made iterative for Archaea and Bacteria
hmmalign --amino --trim -o ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450.sthk ${HMMdir}/PF06450.hmm  ${MMseqs}/PF06450_seq_Step4_e03vsEuk.faa

python parse_stockholm_filter.py ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450.sthk ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450 235
