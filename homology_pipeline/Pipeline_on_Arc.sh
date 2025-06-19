#!/bin/bash -e

#SBATCH --account       uc04105
#SBATCH --job-name      new_pfam_Archaea
#SBATCH --time          72:00:00
#SBATCH --mem           50GB
#SBATCH --cpus-per-task 14
#SBATCH --error         slurm_output_00999/slurm_prokka_%J.err
#SBATCH --output        slurm_output_00999/slurm_prokka_%j.out

PFAM_fasta=/nesi/nobackup/uc04105/fasta_files/pfams/PF00999_Pfam_June_9th_size_excluded_cdhit09.faa
DB=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.faa
module load MMseqs2/15-6f452-gompi-2023a


mmseqs easy-search -e 1.00E-03 --threads 12 ${PFAM_fasta} ${DB} /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/Archaea_vsNewPFAM009999_vsArchaea.tsv  tmp


module load Python/3.11.6-foss-2023a


TSV=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv

python getting_fasta_from_hit_extra.py /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/Archaea_vsNewPFAM009999_vsArchaea.tsv MMSEQ ${TSV} /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/Archaea_vsNewPFAM009999_vsArchaea.fasta

modue load HMMER/3.3.2-GCC-12.3.0

HMMalign=/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea

hmmalign -o ${HMMalign}/Archaea_vsNewPFAM009999_vsArchaea.sthk /nesi/nobackup/uc04105/results/HMM/PF00999.hmm /nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/MMseq/Archaea/Archaea_vsNewPFAM009999_vsArchaea.fasta

parse_stockholm_filter.py ${HMMalign}/Archaea_vsNewPFAM009999_vsArchaea.sthk ${HMMalign}/Archaea_vsNewPFAM009999_vsArchaea.faa 257
