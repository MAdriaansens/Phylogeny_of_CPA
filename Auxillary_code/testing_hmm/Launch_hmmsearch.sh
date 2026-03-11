#!/bin/bash

cd /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMM/PF00999

for file in *hmm; do
        base=$(basename "$file" .hmm);

        sbatch /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/hmmsearch_cluster.sh $base
done
