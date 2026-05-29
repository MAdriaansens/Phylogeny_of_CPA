from Bio import SeqIO
import sys
#the amount of sequences detected in 


Arc_passed_seq = '~/tree_input/fl/Archaea_CPA_fl.fasta'
#taxonomy_bacteria
import os
import json
arc_passed = {}
arc_tax_dic = {}
for record in SeqIO.parse(Arc_passed_seq, 'fasta'):
    arc_passed[record.id] = ''
print(len(arc_passed))
with open('~/Archaea_GTDB226_protein_May92025.tsv', 'r') as test:
    next(test, None)
    for line in test:
        protein_id = line.split('\t')[0]
        if protein_id in arc_passed:
            gtdb_id = line.split('\t')[3]
            taxa = line.split('\t')[2]
            seq = line.split('\t')[-1].split('\n')[0]
            arc_tax_dic[protein_id] = [gtdb_id, taxa, seq]
with open('~Archaea_cpa_fl_taxa.json', 'w')  as output:
    json.dump(arc_tax_dic, output)
