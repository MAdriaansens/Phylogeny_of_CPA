#goal is to get representatives, the red of mrca, the last shared taxonomy of the most distant ones, the distance between them. 
#the average red distance, the number of species in there

#representatives

from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.Consensus import *

import json



#the tree contains the IDs we are interested in 
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/RED_decorated_trees/RED_A53_ALIGNEDRED_30julift.nw', 'newick')

Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_all_seqs.fasta'
Allseq_list = []
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list.append(record.id)
    
termini = tree.get_terminals()
print(len(set(Allseq_list)))
print(len(termini))


Allseq_dict = {}
    #make a dictionairy for representatives and what their GTDB_id is
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv' ,'r') as Proteins:
    next(Proteins, None)
    for protein in Proteins:
        protein_id = protein.split('\t')[0]
        GTDB_id = protein.split('\t')[3]
        protein_tax = protein.split('\t')[2].replace(' ', '_')
        if protein_id in Allseq_list:
           Allseq_dict[protein_id] = 'GTDB_id{}__GTDB_tax:{}'.format(GTDB_id, protein_tax) 
        else:
            pass
print(len(set(list(Allseq_dict.keys()))))


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/All_Arc_CPAseq_GTDB226_ids_tax.json', 'w') as Out:
    json.dump(Allseq_dict, Out)
