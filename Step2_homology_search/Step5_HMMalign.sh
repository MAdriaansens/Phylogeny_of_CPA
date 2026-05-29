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

HMMalign=~/results/HMMalign
MMseqs=~/results/MMseqs/Iterative_MMseqs
HMMdir=/nesi/nobackup/uc04105/results/HMM
mkdir ${MMseqs}
#this code takes the mmseq output of Step 4 and directly runs a hmmalign and filter on it. 
#this is just eukarya, but can be (and has been) made iterative for Archaea and Bacteria similar to HMMsearch etc in previous scripts
#for CPA this was PF00999 and 263 (0.7*375)
hmmalign --amino --trim -o ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450.sthk ${HMMdir}/PF06450.hmm  ${MMseqs}/PF06450_seq_Step4_e03vsEuk.faa

python parse_stockholm_filter.py ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450.sthk ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450_Step5_HMMscan 235
python parse_stockholm_HMMscan.py ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450.sthk ${HMMalign}/Pfam06450_seq_Step4_vsEuk_alignedPF06450_Step5 235
