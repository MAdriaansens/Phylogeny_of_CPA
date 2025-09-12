#goal is to get representatives, the red of mrca, the last shared taxonomy of the most distant ones, the distance between them. 
#the average red distance, the number of species in there

#representatives

from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.Consensus import *

import json

from itertools import combinations

def generate_all_pairs(input_list):
  """
  Generates all unique pairs from a given list.

  Args:
    input_list: The list from which to generate pairs.

  Returns:
    A list of tuples, where each tuple represents a unique pair.
  """
  return list(combinations(input_list, 2))

#the tree contains the IDs we are interested in 
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/RED_decorated_trees/RED_A53_ALIGNEDRED_30julift.nw', 'newick')

Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Archaea_allhmmscanned_final_clustered_at0.7.fasta_all_seqs.fasta'
Allseq_list = []
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list.append(record.id)
    
termini = tree.get_terminals()
print(len(set(Allseq_list)))
print(len(termini))

GTDB_id_tax_dic = {}
import json
#make a dictionairy of all GTDB species in GTDB, with their taxonomy as the key and their GTDB id as the value
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/ar53_metadata.tsv', 'r') as Meta:
    next(Meta, None)
    for line in Meta:
        if (line.split('\t')[18]) == 't':

            GTDB_id = line.split('\t')[0]
            GTDB_tax = line.split('\t')[19].replace(' ', '_')
            GTDB_id_tax_dic[GTDB_tax] = GTDB_id

print(len(GTDB_id_tax_dic.keys()))


import os
file_path = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/All_Arc_seq_GTDB226_tax.json'
if os.path.exists(file_path):
    print('yah')
    
    with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/All_Arc_seq_GTDB226_tax.json', 'r') as F:
        Allseq_dict = json.load(F)

else:
    Allseq_dict = {}
    #make a dictionairy for representatives and what their GTDB_id is
    with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Archaea_GTDB226_protein_May92025.tsv' ,'r') as Proteins:
        next(Proteins, None)
        for protein in Proteins:
            protein_id = protein.split('\t')[0]
            protein_tax = protein.split('\t')[2].replace(' ', '_')
            if protein_id in Allseq_list:
                Allseq_dict[protein_id] = protein_tax
            else:
                pass
    print(len(set(list(Allseq_dict.keys()))))
    
    with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/All_Arc_seq_GTDB226_tax.json', 'w') as outfile:
        json.dump(Allseq_dict, outfile)
print(type(Allseq_dict))
if len(Allseq_dict.keys()) != len(set(Allseq_list)):
    raise ValueError('Dictionary is wrong length of sequences')

for key in Allseq_dict.keys():
    print(key)
    break
MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Archaea_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
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

replist = set(replist)

count = 0
for i in replist:
    count +=(len(Rep_dict[i]))
if 13829 - count != len(singleton_list):
    raise ValueError

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Archaea_clusters_12sep.tsv', 'w') as WRED:

    Header = 'Cluster_representatitve' + '\t' + 'no_proteins_in_cluster' + '\t' + 'red_of_mrca' +  '\t' + 'red_of_mrca2noot' + '\t' + 'average_RED' +  '\t'+ 'common_taxonomic_grouping' + '\t' + 'shared_taxa' + '\t' + 'Phyla_in_cluster' + '\t' + 'no_cross_phyla_pairs' + '\t'+ 'no_cross_class_pairs' + '\t'+ 'no_cross_order_pairs' + '\t' + 'no_cross_family_pairs' + '\t' +  'no_cross_genera_pairs' + '\t' +'list_of_clusters' '\n'
    WRED.write(Header)
    Shared_tax_list = [] 
    sumtax = 0
    
    RED_mrca =0
    RED_mrca_2_root = 0
    average_RED = 0
    for i in replist:
        taxP_list = []
        taxC_list = []
        taxO_list = []
        taxF_list = []
        taxG_list=[]
    
        ids_in_cluster = (Rep_dict[i])
        GTDBids_in_cluster = []
        #number cross taxonomic group in cluster
        for prot_id in ids_in_cluster:
            tax = Allseq_dict[prot_id]
            GTDB_id = GTDB_id_tax_dic[tax]
            GTDBids_in_cluster.append(GTDB_id)
            taxP_list.append(tax.split(';')[1])
            taxC_list.append(tax.split(';')[2])
            taxO_list.append(tax.split(';')[3])
            taxF_list.append(tax.split(';')[4])
            taxG_list.append(tax.split(';')[5])
        
            
        
        
        #latest shared taxonomic group
        if (len(list(set(taxP_list)))) != 1:
            shared_taxonomy = 'Archaea'
            common = 'Archaea'
        elif len(list(set(taxG_list))) == 1:
            common = 'Genera'
            shared_taxonomy = taxG_list[0]
        elif (len(list(set(taxF_list)))) == 1:
            shared_taxonomy = taxF_list[0]
            common = 'Family'
        elif len(list(set(taxO_list))) == 1:
            shared_taxonomy = taxO_list[0]
            common='Order'
        elif len(list(set(taxC_list))) == 1:
            shared_taxonomy = taxC_list[0]
            common='Class'
        elif len(list(set(taxP_list))) == 1:
            shared_taxonomy = taxP_list[0]
            common = 'Phyla'
    
        Cluster_rep = i
        no_proteins_in_cluster = len(ids_in_cluster)
        no_cross_Phyla = str(len(set(taxP_list)))
        no_cross_phyla_pairs=str(len(list(set(taxP_list))))
        no_cross_class_pairs=str(len(list(set(taxC_list))))
        no_cross_order_pairs=str(len(list(set(taxO_list))))
        no_cross_family_pairs=str(len(list(set(taxF_list))))
        no_cross_genera_pairs=str(len(list(set(taxG_list))))
    
    
        distance_list = []
        Allpairs_GTDBids_cluster = generate_all_pairs(GTDBids_in_cluster )
        distance_pairs = 1.1
        distances_list = []
        for pair in Allpairs_GTDBids_cluster:
        
            distance_pairs = (tree.distance(pair[0], pair[1]))
            distances_list.append(distance_pairs)
                
        RED_mrca = (max(distances_list))
        RED_mrca_2_root = 1-(RED_mrca/2)
        average_RED = sum(distances_list)/len(distances_list)
    
    
        Line = (Cluster_rep + '\t' + str(no_proteins_in_cluster) + '\t' + str(RED_mrca) + '\t' + str(RED_mrca_2_root) + '\t' + str(average_RED) + '\t' +  common + '\t' + shared_taxonomy + '\t' + str(set(taxP_list)) + '\t' + no_cross_phyla_pairs + '\t' + no_cross_class_pairs + '\t' + no_cross_order_pairs + '\t' + no_cross_family_pairs + '\t' + no_cross_genera_pairs + '\t' + str(ids_in_cluster) + '\n')
        WRED.write(Line)
         
