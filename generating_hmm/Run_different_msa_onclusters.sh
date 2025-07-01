#!/bin/bash
#SBATCH --job-name=06BacMan50
#SBATCH --time=120:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=48GB          # Memory in MB
#SBATCH --cpus-per-task=27
#SBATCH --array=0-3
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/06Bac_HMM_output%A-%a.out
#SBATCH --error=slurm_output/06Bac_HMM_error%A-%a.err

module load MAFFT/7.505-gimkl-2022a-with-extensions


MMseq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/PF00999
MAFFT=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/PF00999
HMM_out=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMM/Archaea


declare -a array=("0.6" "0.7" "0.8" "0.9")

mafft --thread 25 --auto ${MMseq}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq.fasta > ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_autoaligned.faa

module load HMMER/3.3.2-GCC-12.3.0

hmmbuild ${HMM_out}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_autoaligned.hmm ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_autoaligned.faa

mafft --thread 25 --globalpair ${MMseq}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq.fasta > ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_globalpairaligned.faa

#module load HHMMER/3.3.2-GCC-12.3.0
hmmbuild ${HMM_out}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_globalpairaligned.hmm ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_globalpairaligned.faa

mafft --thread 25 --localpair ${MMseq}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq.fasta > ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_localpairaligned.faa

#module load HMMER/3.3.2-GCC-12.3.0

hmmbuild ${HMM_out}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_localpairaligned.hmm ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_localpairaligned.faa


mafft --thread 25 --genafpair ${MMseq}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq.fasta > ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_genafpairaligned.faa

#module load HMMER/3.3.2-GCC-12.3.0

hmmbuild ${HMM_out}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_genafpairaligned.hmm ${MAFFT}/Archaea_Manual_seq_e05_cov30_seqid${array[$SLURM_ARRAY_TASK_ID]}.faa_rep_seq_genafpairaligned.faa
