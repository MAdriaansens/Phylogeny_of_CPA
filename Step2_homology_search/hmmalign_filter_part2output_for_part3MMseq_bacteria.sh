#!/bin/bash -e
#SBATCH --account       uc04105
#SBATCH --job-name      ArcA
#SBATCH --time          4:00:00
#SBATCH --mem           8GB
#SBATCH --cpus-per-task 1
#SBATCH --error         slurm_outputA/slurm_ArcA_%A-%a.err
#SBATCH --output        slurm_outputA/slurm_ArcA_%A-%a.out
#SBATCH --array         0-123

declare -a array=($(seq 4 127))
HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/part2
HMM=/nesi/nobackup/uc04105/cross_biome_metagenome/Protein/HMM/PF00999.hmm
Seq=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF00999/part2/PF00999Arcseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
module load HMMER/3.4-GCC-12.3.0
module load Python/3.11.6-foss-2023a

hmmalign --amino --trim -o ${HMMalign}/ArcPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMM} ${Seq}

python parse_stockholm_filterFL.py ${HMMalign}/ArcPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMMalign}/ArcPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_FL_e03_mmseq_alignedPF00999.fasta 263

Seq=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF00999/part2/PF00999Eukseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta
hmmalign --amino --trim -o ${HMMalign}/EukPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMM} ${Seq}

python parse_stockholm_filterFL.py ${HMMalign}/EukPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMMalign}/EukPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_FL_e03_mmseq_alignedPF00999.fasta 263


Seq=/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/GTDB_226/results/MMseq/Bacteria/PF00999/part2/PF00999Bacseq_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq_fl_seq.fasta

hmmalign --amino --trim -o ${HMMalign}/BacPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMM} ${Seq}

python parse_stockholm_filterFL.py ${HMMalign}/BacPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_e03_mmseq.sthk ${HMMalign}/BacPF00999_vsBacteria_subset${array[$SLURM_ARRAY_TASK_ID]}_FL_e03_mmseq_alignedPF00999.fasta 263
