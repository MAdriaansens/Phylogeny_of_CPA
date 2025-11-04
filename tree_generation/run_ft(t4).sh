#!/bin/bash -e
#SBATCH --job-name=Iq-tree_TIGR1009_trimmed  # job name (shows up in the queue)
#SBATCH --time=168:00:00      # Walltime (HH:MM:SS
#SBATCH --mem=4GB          # Memory in MB
#SBATCH --cpus-per-task=3
#SBATCH --account=uc04105 
#SBATCH --output=slurm_output/Iqtree_0TRK_output%A.out
#SBATCH --error=slurm_output/Iqtree_0TRK_error%A.err


module load FastTree/2.1.11-GCC-12.3.0
cd /nesi/nobackup/uc04105/new_databases_May/final_tree_set 
 export OMP_NUM_THREADS=3
FastTree -lg -out CPA_TREE_ALGINEDMMSEQ1_PF00999_3Nov_fasttree.treefile second_set/CPA_phylogeny_3novft.fasta
