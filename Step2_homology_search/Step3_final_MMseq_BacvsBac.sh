#!/bin/bash -e
#SBATCH --job-name      BacA
#SBATCH --time          24:00:00
#SBATCH --mem           15GB
#SBATCH --cpus-per-task 10
#SBATCH --error         slurm_outputE/slurm_BacA_%A-%a.err
#SBATCH --output        slurm_outputE/slurm_BacA_%A-%a.out
#SBATCH --array         0-15

declare -a array=($(seq 0 127))

Bac_TSV=~/tsv/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Bac_db=~/Bac_DB/fasta/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
MMseqs=~/results/MMseq/Bacteria/PF00999/part3
Seq=~/results/MMseq/PF00999/PF00999_cross_vsBacteria_dupes_removed_sequences.fasta

#module load MMseqs2/15-6f452-gompi-2023a

#mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Bac_db} ${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp


if [ -f "${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta" ]; then
    echo "${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta"
else
    echo "${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta doesnt exist"
    module load MMseqs2/15-6f452-gompi-2023a
    mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 10 ${Seq}  ${Bac_db} ${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp
    module load Python/3.11.6-foss-2023a
    python getting_fasta_from_hit_extra_Bac.py ${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/Bacfoundseq_vsBacarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.fasta
    rm -r ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp
fi
#

