#!/bin/bash
#SBATCH --job-name=BacSCAN
#SBATCH --time=72:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=10GB          # Memory in MB
#SBATCH --cpus-per-task=10
#SBATCH --array=0-61
#SBATCH --account=uc04105 
#SBATCH --output=slurm_bacoutput/euk6450Run_python_output%A-%a.out
#SBATCH --error=slurm_bacoutput/euk6450Run_python_error%A-%a.err

module load Python/3.11.6-foss-2023a


HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria
declare -a array=($(seq 1 61))

HMMscan=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Bacteria

module load HMMER/3.3.2-GCC-12.3.0

HMMDB=/nesi/nobackup/uc04105/results/HMM/Pfam_db/Merged_PFAMA.hmm
module load HMMER/3.3.2-GCC-12.3.0

#python parse_stockholm_special.py ${HMMalign}/PF03553/cross_domain/PF03553_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}.sthk  ${HMMalign}/PF03553/cross_domain/PF03553_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter 212
#python parse_stockholm_special.py ${HMMalign}/PF03600/cross_domain/PF03600_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}.sthk  ${HMMalign}/PF03600/cross_domain/PF03600_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter 235
#python parse_stockholm_special.py ${HMMalign}/PF06450/cross_domain/PF06450_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}.sthk  ${HMMalign}/PF06450/cross_domain/PF06450_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter 361

hmmscan -E 0.001 --cpu 10 --tblout ${HMMscan}/PF06450/PF06450_Bacteria_merged_alignedPF06450_refilter_${array[$SLURM_ARRAY_TASK_ID]}HMMscanned.tsv ${HMMDB} ${HMMalign}/PF06450/cross_domain/PF06450_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter.fasta

#hmmscan -E 0.001 --cpu 10 --tblout ${HMMscan}/PF03553/PF03553_Bacteroa_merged_alignedPF03553_refilter_${array[$SLURM_ARRAY_TASK_ID]}HMMscanned.tsv ${HMMDB} ${HMMalign}/PF03553/cross_domain/PF03553_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter.fasta

#hmmscan -E 0.001 --cpu 10 --tblout ${HMMscan}/PF03600/PF03600_Bacteroamerged_alignedPF03600_refilter_${array[$SLURM_ARRAY_TASK_ID]}HMMscanned.tsv ${HMMDB} ${HMMalign}/PF03600/cross_domain/PF03600_Bacteria_crossdomain_searched_all_${array[$SLURM_ARRAY_TASK_ID]}_refilter.fasta
