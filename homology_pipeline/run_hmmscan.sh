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

python retrieve_flEUk.py ${OUTDIR} PF03600_Eukarya_merged_alignedPF03600.fasta_exes.fasta PF03600_Eukarya_merged_alignedPF03600_full_length.faa ${TSV}
module load HMMER/3.3.2-GCC-12.3.0

hmmscan --tblout ${HMMscan}/Eukarya_PF03600_aligned_hmmscanned.tsv /nesi/nobackup/uc04105/results/HMM/HMM_db/CPA_IT_pfam ${OUTDIR}/PF03600_Eukarya_merged_alignedPF03600_full_length.faa
python parse_hmmscan.py ${HMMscan}/Eukarya_PF03600_aligned_hmmscanned.tsv PF03600 ${HMMscan}/Eukarya_PF03600_aligned_hmmscanned_onlyPF03600.json
