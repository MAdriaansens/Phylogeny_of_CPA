#!/bin/bash
#SBATCH --account       uc04105
#SBATCH --job-name      06450prokka
#SBATCH --time          72:00:00
#SBATCH --mem           50GB
#SBATCH --cpus-per-task 14
#SBATCH --error         slurm_output_06450/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_06450/slurm_prokka_%A-%a.out


module load Python/3.11.6-foss-2023a
HMMscan=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan
OUTDIR=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign
TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv

module load HMMER/3.3.2-GCC-12.3.0

hmmscan -E 0.001 -T 12 --tblout ${HMMscan}/Eukarya_PF03600_aligned_hmmscanned_Step6.tsv /nesi/nobackup/uc04105/results/HMM/HMM_db/CPA_IT_pfam ${OUTDIR}/Pfam06450_seq_Step4_vsEuk_alignedPF06450_Step5.fasta
