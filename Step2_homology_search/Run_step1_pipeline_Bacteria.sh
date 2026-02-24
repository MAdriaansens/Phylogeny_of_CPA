#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      BacA
#SBATCH --time          102:00:00
#SBATCH --mem           100GB
#SBATCH --cpus-per-task 25
#SBATCH --error         slurm_outputB/slurm_BacA_%A-%a.err
#SBATCH --output        slurm_outputB/slurm_BacA_%A-%a.out
#SBATCH --array         0-61

declare -a array=($(seq 0 61))

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM
Bac_TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bac_DB/tsv/Bacteria_GTDB226_protein_May92025_chunk_${array[$SLURM_ARRAY_TASK_ID]}.tsv
Bac_db=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bac_DB/fasta/Bacteria_GTDB226_protein_May92025_subset${array[$SLURM_ARRAY_TASK_ID]}.fasta
HMMsearch=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/HMMsearch/Bacteria
MMseqs=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/MMseqs/Bacteria
Seq=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/sequences

#---------------------------------------------------------NhaB-----------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF06450_NhaB_sequences.fasta  ${Bac_db} ${MMseqs}/PF06450_NhaBvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF06450_NhaBvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF06450_NhaBvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF06450_NhaBhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF06450_NhaB.hmm ${Bac_db}
module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF06450_NhaBhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Bac_TSV} ${HMMsearch}/PF06450_NhaBhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta


#---------------------------------------------------------NhaC------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF03553_NhaC_sequences.fasta  ${Bac_db} ${MMseqs}/PF03553_NhaCvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF03553_NhaCvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF03553_NhaCvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF03553_NhaChmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF03553_NhaC.hmm ${Bac_db}

module load Python/3.11.3-gimkl-2022a
python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF03553_NhaChmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Bac_TSV} ${HMMsearch}/PF03553_NhaCHmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta



#---------------------------------------------------------NhaD------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 20 ${Seq}/PF03600_NhaD_sequences.fasta  ${Bac_db} ${MMseqs}/PF03600_NhaDvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv ${array[$SLURM_ARRAY_TASK_ID]}_Bac_tmp

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${MMseqs}/PF03600_NhaDvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.tsv MMSEQ ${Bac_TSV} ${MMseqs}/PF03600_NhaDvsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 20 -E 0.001 --tblout ${HMMsearch}/PF03600_NhaDhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv ${HMMdir}/PF03600_NhaD.hmm ${Bac_db}

module load Python/3.11.3-gimkl-2022a

python getting_fasta_from_hit_extra_Arc.py ${HMMsearch}/PF03600_NhaDhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03.tsv HMM ${Bac_TSV} ${HMMsearch}/PF03600_NhaDhmmvs_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_fl_sequence.fasta
