for file in *fasta; do
 base=$(basename "$file" .fasta);
 hmmbuild ~/results/HMM/${base}.hmm $file;
 done
