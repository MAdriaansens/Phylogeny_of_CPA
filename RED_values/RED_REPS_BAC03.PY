
print('test')

from Bio import SeqIO

from Bio.Phylo.Consensus import *
import os

import json


Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_all_seqs.fasta'
Allseq_list = {}

print('parsing sequencces')
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list[str(record.id)] = record.id




MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'
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
print(len(singleton_list))replist
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

replist = set(replist)
nummers
small_count = 0
larger_count = 1
sub_dict ={}
for key in Rep_dict.keys():
    small_count +=1
    sub_dict[key] = Rep_dict[key]
    if small_count % 1523 == 0:
        with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/dictionaries/Bacteria_red_reps_subset_{}.json'.format(larger_count), 'w') as JSON:
            json.dump(sub_dict, JSON)
            JSON.close()
            sub_dict = {}
            larger_count +=1
#import os

#for dic in os.listdir('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/dictionaries'):
#    with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/dictionaries/{}'.format(dic), 'r') as test:
#        run = json.load(test)
#        for key in run.keys():
#            print(key, dic)
#            break
