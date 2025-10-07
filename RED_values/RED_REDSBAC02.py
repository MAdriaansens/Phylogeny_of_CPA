import json
import sys

from Bio import SeqIO

diction_in = sys.argv[1]

Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Bacteria_allhmmscanned_final_clustered_at0.7.fasta_all_seqs.fasta'
Allseq_list = []

print('parsing sequences')
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list.append(str(record.id))
Allseq_list = set(Allseq_list)


Seq_tax_dic = {}
Seq_tax = diction_in.split('.json')[0].split('_tax_')[1]

with open('{}'.format(diction_in), 'r') as F:
    Allseq_dict = json.load(F)
    for Seq in Allseq_list:
        if Seq in Allseq_dict:
            Seq_tax_dic[Seq] = Allseq_dict[Seq]
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/All_Bac_CPAseq_GTDB226_ids_tax_{}.json'.format(Seq_tax), 'w') as X:
    json.dump(Seq_tax_dic, X)
