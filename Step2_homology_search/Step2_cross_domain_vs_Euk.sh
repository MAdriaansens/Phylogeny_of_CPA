#!/bin/bash -e
#SBATCH --job-name      EukA
#SBATCH --time          8:00:00
#SBATCH --mem           30GB
#SBATCH --cpus-per-task 20
#SBATCH --error         slurm_outputE/slurm_EukA_%A-%a.err
#SBATCH --output        slurm_outputE/slurm_BacA_%A-%a.out
#SBATCH --array         0-15

declare -a array=($(seq 0 15))

HMMdir=~/HMM/PF00999.hmm

#in this case we decided to split up Eukarya into 15 subsets. 
Euk_TSV=~/Euk_DB/tsv/Euk_db_May_protein_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Euk_db=~/Euk_DB/fasta/Euk_db_May_protein_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta

#else
#Euk_TSV= ~/Euk_db_7April_protein.tsv
#Euk_db=~ ~/Euk_db_7April.fasta

HMMsearch=~/results/HMMsearch/PF00999
MMseqs=~/results/MMseq/PF00999/part2
Seq=~/results/Bacteria_all_part1A_remove_dupes.fasta

module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Euk_db} ${MMseqs}/PF00999Bacseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp

Seq=~/results/Eukarya_part1A_dupes_removed_sequences.fasta
mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Euk_db} ${MMseqs}/PF00999Eukseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp


Seq=~/results/Archaea_all_part1A_remove_dupes.fasta
mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Euk_db} ${MMseqs}/PF00999Arcseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Euk_tmp

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF00999Arcseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF00999Arcseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF00999Bacseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF00999Bacseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta


python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF00999Eukseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF00999Eukseq_vsEukarya_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
