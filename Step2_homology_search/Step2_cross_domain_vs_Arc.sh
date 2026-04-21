#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      ArcA
#SBATCH --time          2:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 20
#SBATCH --error         slurm_outputA/slurm_ArcA_%A-%a.err
#SBATCH --output        slurm_outputA/slurm_ArcA_%A-%a.out
#SBATCH --array         0-33

declare -a array=($(seq 0 33))

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM/PF00999.hmm
Arc_TSV=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/GTDB_226/ARCDB/Archaea_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Arc_db=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/ARCDB/Archaea_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
HMMsearch=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea
MMseqs=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/PF00999/part2
Seq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Bacteria_all_part1A_remove_dupes.fasta


#---------------------CPA-------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq} ${Arc_db} ${MMseqs}/PF00999Bacseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Arc_tmp
Seq=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/Archaea_all_part1A_remove_dupes.fasta

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq} ${Arc_db} ${MMseqs}/PF00999Arcseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Arc_tmp
Seq=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/results/Eukarya_part1A_dupes_removed_sequences.fasta

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq} ${Arc_db} ${MMseqs}/PF00999Eukseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Arc_tmp

module load Python/3.11.6-foss-2023a
python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Arcseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF00999Arcseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Eukseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF00999Eukseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF00999Bacseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF00999Bacseq_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
