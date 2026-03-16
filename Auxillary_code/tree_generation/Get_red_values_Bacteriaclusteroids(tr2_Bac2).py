print('test')
#first run Get_red_values_Bacteriaclusteroids(tr2_Bac2).py

from Bio import Phylo
import sys
from Bio import SeqIO

from Bio.Phylo.Consensus import *
import os

import json
from itertools import combinations
print('loaded packages')

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

Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Bacteria_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_all_seqs.fasta'
Allseq_list = {}

print('parsing sequencces')
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list[str(record.id)] = record.id



#representatives part


import json
import os

Allseq_dict = {}
for diction in os.listdir('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/'):
    if diction.split('.')[1] == 'json':
        if 'All_Bac' in diction:
            with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/{}'.format(diction), 'r') as F:
                Tax_dic = json.load(F)
                for key in Tax_dic.keys():
                    Allseq_dict[key] = Tax_dic[key]

print('reading tree')
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/RED_decorated_trees/RED_Bac120_ALIGNEDRED.nw', 'newick')
print('finsihed reading tree')
count = 0
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.name = str(count)
        count = count + 1
REPD = sys.argv[1]

#first run Get_red_values_Bacteriaclusteroids(tr2_Bac2).py
with open('{}'.format(REPD), 'r') as REPRED:
    rep_diction = json.load(REPRED)
rep_dic_id = REPD.split('.json')[0].split('subset_')[1]

rep_list = list(rep_diction.keys())

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/dictionaries/RED_values_Bacteria_clusters_8okt_rep_dic_id{}.tsv'.format(rep_dic_id), 'w') as WRED:
    Header = 'Cluster_representatitve' + '\t' + 'no_proteins_in_cluster' + '\t' + 'mrca' + '\t' + 'red_of_mrca' +  '\t' + 'red_of_mrca2noot' + '\t' + 'average_RED' +  '\t'+ 'common_taxonomic_grouping' + '\t' + 'shared_taxa' + '\t' + 'Phyla_in_cluster' + '\t' + 'no_cross_phyla_pairs' + '\t'+ 'no_cross_class_pairs' + '\t'+ 'no_cross_order_pairs' + '\t' + 'no_cross_family_pairs' + '\t' +  'no_cross_genera_pairs' + '\t' +'list_of_clusters' + '\t' + 'GTDB_id_in_cluster' + '\t' + 'most_distant_pairs' + '\n'
    WRED.write(Header)
    Shared_tax_list = []
    sumtax = 0

    RED_mrca =0
    RED_mrca_2_root = 0
    average_RED = 0
    for i in rep_list:
        taxP_list = []
        taxC_list = []
        taxO_list = []
        taxF_list = []
        taxG_list=[]

        ids_in_cluster = (rep_diction[i])
        GTDBids_in_cluster = []
        #number cross taxonomic group in cluster
        for prot_id in ids_in_cluster:
            tax = Allseq_dict[prot_id].split('GTDB_tax:')[1]
            GTDB_id = Allseq_dict[prot_id].split('__GTDB_tax:')[0].split('GTDB_id')[1]
            GTDBids_in_cluster.append(GTDB_id)
            taxP_list.append(tax.split(';')[1])
            taxC_list.append(tax.split(';')[2])
            taxO_list.append(tax.split(';')[3])
            taxF_list.append(tax.split(';')[4])
            taxG_list.append(tax.split(';')[5])




        #latest shared taxonomic group
        if (len(list(set(taxP_list)))) != 1:
            shared_taxonomy = 'Bacteria'
            common = 'Bacteria'
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
            if distance_pairs > RED_mrca:
                RED_mrca =  distance_pairs
                most_distant_pairs = (pair[0], pair[1])
                mrca = tree.common_ancestor(pair[0], pair[1])
            distances_list.append(distance_pairs)
        
        RED_mrca_2_root = 1-(RED_mrca/2)
        average_RED = sum(distances_list)/len(distances_list)


        Line = (Cluster_rep + '\t' + str(no_proteins_in_cluster) + '\t' + mrca.name + '\t' + str(RED_mrca) + '\t' + str(RED_mrca_2_root) + '\t' + str(average_RED) + '\t' +  common + '\t' + shared_taxonomy + '\t' + str(set(taxP_list)) + '\t' + no_cross_phyla_pairs + '\t' + no_cross_class_pairs + '\t' + no_cross_order_pairs + '\t' + no_cross_family_pairs + '\t' + no_cross_genera_pairs + '\t' + str(ids_in_cluster) +  '\t' + str(GTDBids_in_cluster) + '\t' + str(most_distant_pairs) + '\n')
        WRED.write(Line)
