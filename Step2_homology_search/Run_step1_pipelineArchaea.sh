#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      ArcA
#SBATCH --time          102:00:00
#SBATCH --mem           100GB
#SBATCH --cpus-per-task 25
#SBATCH --error         slurm_outputA/slurm_ArcA_%A-%a.err
#SBATCH --output        slurm_outputA/slurm_ArcA_%A-%a.out
#SBATCH --array         0-12

declare -a array=($(seq 0 12))

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM
Arc_TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv
Arc_db=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/ARCDB/Archaea_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
HMMsearch=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/HMMsearch/Archaea
MMseqs=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/MMseqs/Archaea
Seq=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/sequences


#---------------------------------------------------------NhaC------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF03553_NhaC_sequences.fasta  ${Arc_db} ${MMseqs}/PF03553_NhaCvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF03553_NhaCvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF03553_NhaCvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF03553_NhaChmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF03553_NhaC.hmm ${Arc_db}

module load Python/3.11.3-gimkl-2022a
python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF03553_NhaChmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Arc_TSV} ${HMMsearch}/PF03553_NhaCHmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta

#---------------------------------------------------------NhaB-----------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF06450_NhaB_sequences.fasta  ${Arc_db} ${MMseqs}/PF06450_NhaBvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF06450_NhaBvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF06450_NhaBvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF06450_NhaBhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF06450_NhaB.hmm ${Arc_db}
module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF06450_NhaBhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Arc_TSV} ${HMMsearch}/PF06450_NhaBhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta


#---------------------------------------------------------NhaD------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF03600_NhaD_sequences.fasta  ${Arc_db} ${MMseqs}/PF03600_NhaDvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF03600_NhaDvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Arc_TSV} ${MMseqs}/PF03600_NhaDvsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF03600_NhaDhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF03600_NhaD.hmm ${Arc_db}

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF03600_NhaDhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Arc_TSV} ${HMMsearch}/PF03600_NhaDhmmvs_vsArchaea_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta

