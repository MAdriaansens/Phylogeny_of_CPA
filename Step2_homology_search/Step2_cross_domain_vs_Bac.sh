#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      BacA
#SBATCH --time          44:00:00
#SBATCH --mem           40GB
#SBATCH --cpus-per-task 20
#SBATCH --error         slurm_outputB/slurm_BacA_%A-%a.err
#SBATCH --output        slurm_outputB/slurm_BacA_%A-%a.out
#SBATCH --array         0-127%40

declare -a array=($(seq 0 127))

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM/PF00999.hmm
Bac_TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bac_DB/tsv/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Bac_db=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bac_DB/fasta/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria
MMseqs=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF00999/part2
Seq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Bacteria_all_part1A_remove_dupes.fasta

module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Bac_db} ${MMseqs}/PF00999Bacseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

rm -r ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp
Seq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Archaea_all_part1A_remove_dupes.fasta
mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Bac_db} ${MMseqs}/PF00999Arcseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

rm -r ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp
Seq=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/Eukarya_part1A_dupes_removed_sequences.fasta
mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}  ${Bac_db} ${MMseqs}/PF00999Eukseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

rm -r ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

#this retrieving sequences can be run in the same script, but to save memmory we had a seperate bash script with 12 GB and 1 cpu per task which was more efficient memmory and CPU wise. 

#module load Python/3.11.6-foss-2023a
#python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Bacseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF00999Bacseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
#python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Eukseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF00999Eukseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
#python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Arcseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF00999Arcseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
