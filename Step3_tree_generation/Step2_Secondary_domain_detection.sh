HMM=/nesi/nobackup/uc04105/results/HMM/pfams
HMMalign=/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign

module load Python/3.11.6-foss-2023a

module load HMMER/3.3.2-GCC-12.3.0
HMMid=06965
HMMcutoff=258
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=02254
HMMcutoff=82
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=02080
HMMcutoff=50
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=23256
HMMcutoff=99

hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=23259
HMMcutoff=105
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=00027
HMMcutoff=63
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=08619
HMMcutoff=258
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}

HMMid=00582
HMMcutoff=99
hmmalign -o ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMM}/PF${HMMid}.hmm /nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_all_hmmscanned_aligned_fl.fasta

python parse_stockholm_filter.py ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.sthk ${HMMalign}/${HMMid}_vsAllEukarya_${HMMid}_aligned.fa ${HMMcutoff}
