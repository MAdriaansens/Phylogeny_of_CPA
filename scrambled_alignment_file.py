from Bio import SeqIO
import json
import os
import random

id_list = []
prot_seq ={}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/CPA_Phylogeny_redone_6Okt_PF00999_aligned_hmmscanned_scrambled.fasta', 'w') as X:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/CPA_Phylogeny_redone_6Okt_PF00999_aligned_hmmscanned.fasta', 'fasta'):
        id_list.append(record.id)
        prot_seq[record.id] = str(record.seq)
    random.shuffle(id_list)
    for i in id_list:
        sequence = '>' + i + '\n' + prot_seq[i] + '\n'

        X.write(sequence)
