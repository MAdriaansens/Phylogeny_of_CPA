import json
import os
from Bio import SeqIO
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain/Arc_cross_domain_hmmaligned_fulllength_pf00999.faa', 'fasta'):
    Arc_dict[record.id.split('_tax')[0]] = str(record.seq)

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta', 'w') as WA:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq.fasta', 'fasta'):
        if record.id in Arc_dict:
            Line = '>' + record.id + '\n' + Arc_dict[record.id] + '\n'
            WA.write(Line)
        

Bac_dict = {}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/cross_domain_final/All_merged_Bacteria_FL.faa', 'fasta'):
    Bac_dict[record.id.split('_tax')[0]] = str(record.seq)

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta', 'w') as WA:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq.fasta', 'fasta'):
        if record.id in Bac_dict:
            Line = '>' + record.id + '\n' + Bac_dict[record.id] + '\n'
            WA.write(Line)

Euk_dict ={}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/All_merged_PF00999.faa', 'fasta'):
    Euk_dict[record.id.split('_tax')[0]] = str(record.seq)
    
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq_fl.fasta', 'w') as WA:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/Eukarya_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_rep_seq.fasta', 'fasta'):
        if record.id in Euk_dict:
            Line = '>' + record.id + '\n' + Euk_dict[record.id] + '\n'
            WA.write(Line)
    
