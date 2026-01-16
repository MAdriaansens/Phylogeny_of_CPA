from Bio import SeqIO
import sys

fasta_chunk = sys.argv[1]

#the amount of sequences detected in 
#Eukarya 3414
#Bacteria 341212
#Archaea 14101
Euk_passed_seq = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/Eukarya_4nov_passed_all_setcpa.fasta'

#taxonomy_bacteria
import os
import json

euk_passed = {}

for record in SeqIO.parse(Euk_passed_seq, 'fasta'):
    euk_passed[record.id] = ''

with open('/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_ids_seq_tax_tax.json', 'r') as euk_tax:
    euk_dic = json.load(euk_tax)
    for key in euk_dic.keys():
        if key in euk_passed:
            euk_passed[key] = euk_dic[key]
euk_tax_seq_dic = {}
for record in SeqIO.parse('/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/Euk_database_May/Euk_db_May_protein.fasta', 'fasta'):
    if str(record.id) in euk_passed:
        tax = euk_passed[record.id]
        seq = str(record.seq)
        euk_tax_seq_dic[record.id] = [tax, seq]
print(len(euk_tax_seq_dic.keys()))


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/fl_taxa/Eukarya/Eukarya_cpa_fl_taxa.json', 'w')  as output:
    json.dump(euk_tax_seq_dic, output)
