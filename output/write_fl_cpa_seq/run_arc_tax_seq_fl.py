from Bio import SeqIO
import sys
#the amount of sequences detected in 
#Archaea 14101
Arc_passed_seq = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Archaea_passed_all_filters_7OKT_GTDBreps_alignedPF00999.fasta'
#taxonomy_bacteria
import os
import json
arc_passed = {}
arc_tax_dic = {}
for record in SeqIO.parse(Arc_passed_seq, 'fasta'):
    arc_passed[record.id] = ''
print(len(arc_passed))
with open('/home/mad149/00_nesi_projects/uc04105_nobackup/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv', 'r') as test:
    next(test, None)
    for line in test:
        protein_id = line.split('\t')[0]
        if protein_id in arc_passed:
            gtdb_id = line.split('\t')[3]
            taxa = line.split('\t')[2]
            seq = line.split('\t')[-1].split('\n')[0]
            arc_tax_dic[protein_id] = [gtdb_id, taxa, seq]
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/second_set/fl_taxa/Archaea/Archaea_cpa_fl_taxa.json', 'w')  as output:
    json.dump(arc_tax_dic, output)
