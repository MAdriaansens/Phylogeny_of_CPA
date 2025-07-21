#!/bin/bash -e

#SBATCH --account       uc04105
#SBATCH --job-name      Hmmscan
#SBATCH --time          72:00:00
#SBATCH --mem           15GB
#SBATCH --cpus-per-task 17
#SBATCH --error         pipeline_scan/slurm_prokka_%A-%a.err
#SBATCH --output        pipeline_scan/slurm_prokka_%A-%a.out
#SBATCH --array         1-12

declare -a array=($(seq 1 12))
HMMALGIN=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign
HMMDB=/nesi/nobackup/uc04105/results/HMM/Pfam_db/Merged_PFAMA.hmm
module load HMMER/3.3.2-GCC-12.3.0

hmmscan --cpu 10 -E 0.001 --tblout /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea/All_Arcseq_7juli_vsArchaea_subset_${array[$SLURM_ARRAY_TASK_ID]}FULLPFAMscanned.tsv ${HMMDB}  /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/All_Arcseq_7juli_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}.faa

hmmscan --cpu 10 -E 0.001 --tblout /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea/All_Bacseq_7juli_vsArchaea_subset_${array[$SLURM_ARRAY_TASK_ID]}FULLPFAMscanned.tsv ${HMMDB}  /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/All_Bacseq_7juli_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}.faa
