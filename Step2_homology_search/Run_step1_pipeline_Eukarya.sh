#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      EukA
#SBATCH --time          102:00:00
#SBATCH --mem           100GB
#SBATCH --cpus-per-task 25
#SBATCH --error         slurm_output/slurm_EukA_%A.err
#SBATCH --output        slurm_output/slurm_EukA_%A.out

HMMdir=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM
Euk_TSV=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv
Euk_db=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta
HMMsearch=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/HMMsearch/Eukarya
MMseqs=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/results/MMseqs/Eukarya
Seq=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/sequences
#---------------------------------------------------------NhaB-----------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 10 ${Seq}/PF06450_NhaB_sequences.fasta  ${Euk_db} ${MMseqs}/PF06450_NhaBvsEukarya_e03_mmseq.tsv tmp

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF06450_NhaBvsEukarya_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF06450_NhaBvsEukarya_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 10 -E 0.001 --tblout ${HMMsearch}/PF06450_NhaBhmmvs_Eukarya_e03.tsv ${HMMdir}/PF06450_NhaB.hmm ${Euk_db}

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${HMMsearch}/PF06450_NhaBhmmvs_Eukarya_e03.tsv HMM ${Euk_TSV} ${HMMsearch}/PF06450_NhaBhmmvs_Eukarya_e03_fl_sequence.fasta


#---------------------------------------------------------NhaC------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 10 ${Seq}/PF03553_NhaC_sequences.fasta  ${Euk_db} ${MMseqs}/PF03553_NhaCvsEukarya_e03_mmseq.tsv tmp

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF03553_NhaCvsEukarya_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF03553_NhaCvsEukarya_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 10 -E 0.001 --tblout ${HMMsearch}/PF03553_NhaChmmvs_Eukarya_e03.tsv ${HMMdir}/PF03553_NhaC.hmm ${Euk_db}

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${HMMsearch}/PF03553_NhaChmmvs_Eukarya_e03.tsv HMM ${Euk_TSV} ${HMMsearch}/PF03553_NhaCHmmvs_Eukarya_e03_fl_sequence.fasta



#---------------------------------------------------------NhaD------------------------------------------------------------------------------------------------
module load MMseqs2/15-6f452-gompi-2023a

mmseqs easy-search -e 1.00E-03 -c 0.0 --threads 10 ${Seq}/PF03600_NhaD_sequences.fasta  ${Euk_db} ${MMseqs}/PF03600_NhaDvsEukarya_e03_mmseq.tsv tmp

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${MMseqs}/PF03600_NhaDvsEukarya_e03_mmseq.tsv MMSEQ ${Euk_TSV} ${MMseqs}/PF03600_NhaDvsEukarya_e03_mmseq_fl_seq.fasta

module purge

module load HMMER/3.4-GCC-12.3.0

#HMMsearch
hmmsearch --noali --cpu 10 -E 0.001 --tblout ${HMMsearch}/PF03600_NhaDhmmvs_Eukarya_e03.tsv ${HMMdir}/PF03600_NhaD.hmm ${Euk_db}

module load Python/3.11.6-foss-2023a

python getting_fasta_from_hit_extra_Euk.py ${HMMsearch}/PF03600_NhaDhmmvs_Eukarya_e03.tsv HMM ${Euk_TSV} ${HMMsearch}/PF03600_NhaDhmmvs_Eukarya_e03_fl_sequence.fasta
