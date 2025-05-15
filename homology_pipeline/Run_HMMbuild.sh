for file in *fasta; do
 base=$(basename "$file" .fasta);
 hmmbuild /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMM/${base}.hmm $file;
 done
