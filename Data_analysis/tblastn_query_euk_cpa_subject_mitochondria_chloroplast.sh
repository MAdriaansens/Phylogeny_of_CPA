#!/bin/bash
#SBATCH --account       uc04105
#SBATCH --job-name      Euk_blastn
#SBATCH --time          2:00:00
#SBATCH --mem           2GB
#SBATCH --array         0-32
#SBATCH --cpus-per-task 2
#SBATCH --error         slurm_output_00999/slurm_prokka_%A-%a.err
#SBATCH --output        slurm_output_00999/slurm_prokka_%A-%a.out
genomedir=/nesi/nobackup/uc04105/new_databases_May/chloroplast_mitochondria_genomes/database
seed_dir=/nesi/nobackup/uc04105/new_databases_May/chloroplast_mitochondria_genomes/euk_cpa_seeds
tblastnout=/nesi/nobackup/uc04105/new_databases_May/chloroplast_mitochondria_genomes/blastn_output
module load BLAST/2.16.0-GCC-12.3.0

declare -a array=("Q56XP4_NHX2" "Q68KI4_NHX1" "Q8VYR9_KEA5" "Q9ZUN3_KEA4" "Q9LKW9_NHX7" "Q84WG1_NHX3" "Q8S397_NHX4" "Q9ZUV9_CHX7" "Q58P69_CHX10" "Q9M353_CHX20" "Q8VYD4_CHX23" "Q8GX92_CHX6A" "P0CG16_CHX6B" "Q9SIT5_CHX15" "Q9LUN4_CHX19" "Q9FYB9_CHX11" "Q9FFB7_CHX9" "Q8S396_NHX5" "Q8RWU6_NHX6" "Q9SUQ7_CHX17" "Q9SA37_CHX1" "Q58P71_CHX8" "O22920_CHX13" "Q9FFR9_CHX18" "Q9SKA9_CHX21" "Q9M008_CHX26" "Q8L709_CHX28" "Q9M0Z3_KEA3" "O65272_KEA2" "Q9ZTZ7_KEA1" "B5X0N6_KEA6" "Q3yl57_NHX8")

tblastn -query ${seed_dir}/Arabidopsis_thaliana_${array[$SLURM_ARRAY_TASK_ID]}.fna -subject ${genomedir}/Arabidopsis_thaliana_chrMT.fna -out ${tblastnout}/Arabidopsis_thaliana_chrMT_vs_${array[$SLURM_ARRAY_TASK_ID]}.tsv -outfmt 6


tblastn -query ${seed_dir}/Arabidopsis_thaliana_${array[$SLURM_ARRAY_TASK_ID]}.fna -subject ${genomedir}/Arabidopsis_thaliana_chrPltd.fna -out ${tblastnout}/Arabidopsis_thaliana_chrPltd_vs_${array[$SLURM_ARRAY_TASK_ID]}.tsv -outfmt 6

