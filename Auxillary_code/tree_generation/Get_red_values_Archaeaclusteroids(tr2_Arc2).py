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
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.name = str(count)
        count = count + 1
        
Allseqs = '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_all_seqs.fasta'
Allseq_list = []
for record in SeqIO.parse(Allseqs, 'fasta'):
    Allseq_list.append(record.id)
    
termini = tree.get_terminals()
print(len(set(Allseq_list)))
print(len(termini))
#we need a seperate script (Get_red_values_Archaeaclusteroids(tr2_Arc1).py for this
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/All_Arc_CPAseq_GTDB226_ids_tax.json', 'r') as F:
    Allseq_dict = json.load(F)


print(len(Allseq_dict.keys()))
for key in Allseq_dict.keys():
    print(Allseq_dict[key])
    break




MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/prefered_sequences/MMseq/Archaea_all_hmmscanned_aligned_treeinput_clusterd_at_0.7.faa_cluster.tsv'
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


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Archaea_clusters_7OKT.tsv', 'w') as WRED:

    Header = 'Cluster_representatitve' + '\t' + 'no_proteins_in_cluster' + '\t' + 'mrca' + '\t' + 'red_of_mrca' +  '\t' + 'red_of_mrca2noot' + '\t' + 'average_RED' +  '\t'+ 'common_taxonomic_grouping' + '\t' + 'shared_taxa' + '\t' + 'Phyla_in_cluster' + '\t' + 'no_cross_phyla_pairs' + '\t'+ 'no_cross_class_pairs' + '\t'+ 'no_cross_order_pairs' + '\t' + 'no_cross_family_pairs' + '\t' +  'no_cross_genera_pairs' + '\t' +'list_of_clusters' + '\t' + 'GTDB_id_in_cluster' + '\t' + 'most_distant_pairs' + '\n'
    WRED.write(Header)
    Shared_tax_list = [] 
    sumtax = 0
    
    RED_mrca =0
    RED_mrca_2_root = 0
    average_RED = 0
    for i in replist:
        RED_mrca =0
        RED_mrca_2_root = 0
        average_RED = 0
        
        taxP_list = []
        taxC_list = []
        taxO_list = []
        taxF_list = []
        taxG_list=[]
    
        ids_in_cluster = (Rep_dict[i])
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
            if distance_pairs > RED_mrca:
                RED_mrca =  distance_pairs
                most_distant_pairs = (pair[0], pair[1])
                mrca = tree.common_ancestor(pair[0], pair[1])
            distances_list.append(distance_pairs)
        
        RED_mrca_2_root = 1-(RED_mrca/2)
        average_RED = sum(distances_list)/len(distances_list)
    
    
        Line = (Cluster_rep + '\t' + str(no_proteins_in_cluster) + '\t' + mrca.name + '\t' + str(RED_mrca) + '\t' + str(RED_mrca_2_root) + '\t' + str(average_RED) + '\t' +  common + '\t' + shared_taxonomy + '\t' + str(set(taxP_list)) + '\t' + no_cross_phyla_pairs + '\t' + no_cross_class_pairs + '\t' + no_cross_order_pairs + '\t' + no_cross_family_pairs + '\t' + no_cross_genera_pairs + '\t' + str(ids_in_cluster) +  '\t' + str(GTDBids_in_cluster) + '\t' + str(most_distant_pairs) + '\n')
        WRED.write(Line)
         
