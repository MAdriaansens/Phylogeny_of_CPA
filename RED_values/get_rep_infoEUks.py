from Bio import SeqIO
import json
import os
Euk_meta = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv'

file_path = '/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax.json'

protein_ids_all = []

for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Euk_allhmmscanned_final_clustered_at0.7.fasta_all_seqs.fasta', 'fasta'):
    protein_ids_all.append(record.id)
if os.path.exists(file_path):
    print('yah')
    
    with open('/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax.json', 'r') as F:
       data = json.load(F)
else:
    data = {}
    #make a dictionairy for representatives and what their GTDB_id is
    with open(Euk_meta, 'r') as Euks:
        for Euk in Euks:
            protein_id = Euk.split('\t')[0]
            if protein_id in protein_ids_all:
                Tax_inf = Euk.split('\t')[1] +'_Taxa:' + Euk.split('\t')[4]
                data[protein_id] = Tax_inf 
    with open('/nesi/nobackup/uc04105/new_databases_May/Euk_ids_seq_tax_tax.json', 'w') as J:
        json.dump(data, J)



MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Euk_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
rep_list = []


with open(MMseq, 'r') as clusters:
    for cluster in clusters:
        rep = cluster.split('\t')[0]
        rep_list.append(rep)
    clusters.close()

from collections import Counter

#the rep list, is just a list of all values in the representative columns, this means that ids are present multiple times

rep_list = (rep_list)
reps = Counter(rep_list)

#reps returns is a counter of how many sequences a rep represents
items = reps.items()


singleton_list = []
multiple_list = []

for i in items:
    if i[-1] == 1:
        singleton_list.append(i[0])
    else:
        multiple_list.append(i[0])

print(len(multiple_list))
print(len(singleton_list))
replist = []
Rep_dict = {}

with open(MMseq, 'r') as clustermap:
    count = 0
    second_list = []
    element_count = 0
    for entry in clustermap:
        count = reps[entry.split('\t')[0]]
        repm = entry.split('\t')[0]

        if repm in multiple_list:
            if repm in replist:
                second = (entry.split('\t')[1].split('\n')[0])
                second_list.append(second)
                if len(second_list) == count:
                    Rep_dict[repm] = second_list
            else:

                replist.append(repm)
                second_list = []
                second_list.append(repm)
        else:
            pass
print(len(set(replist)))


