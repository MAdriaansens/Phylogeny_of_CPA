from Bio import SeqIO
import sys

fasta_chunk = sys.argv[1]

#the amount of sequences detected in 
#Bacteria 341212
Bac_passed_seq = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/Bacteria_passed_all_filters_OKT3_GTDBreps_alignedPF00999.fasta'

#taxonomy_bacteria
import os
import json

bac_tax_dic = {}

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/MMseq/All_Bac_CPAseq_GTDB226_ids_tax_{}.json'.format(fasta_chunk), 'r') as test:
     bac_tax_dic.update(json.load(test))
print(len(bac_tax_dic.keys()))

bac_tax = {}
for record in SeqIO.parse(Bac_passed_seq, 'fasta'):
    if record.id in bac_tax_dic:
        tax_info = (bac_tax_dic[record.id]).replace('GTDB_id', 'GTDB_id:')
        bac_tax[str(record.id)] =(tax_info)
chunk_dic = {}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bac_DB/fasta/Bacteria_GTDB226_protein_May92025_subset{}.fasta'.format(fasta_chunk), 'fasta'):
    if record.id in bac_tax:
        taxa = bac_tax[str(record.id)]
        seq = (str(record.seq))
        chunk_dic[str(record.id) + '_' + taxa] = seq
        
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/fl_taxa/Bacteria_cpa_fl_taxa_subset{}.json'.format(fasta_chunk), 'w') as f:
    json.dump(chunk_dic, f)

