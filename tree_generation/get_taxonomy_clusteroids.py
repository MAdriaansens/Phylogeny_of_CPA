#important in analysis, this tells us if certain clusteroids are made up of different taxonomic groupings
from Bio import SeqIO

#diction contains all the protein_ids and matching accession id + taxonomy


#get signles and multiples:


Tax_dict = {}
import os
for entry in os.listdir('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain'):
    if 'faa.fasta' in entry:
        print(entry)
        for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain/{}'.format(entry), 'fasta'):
            Tax_dict[record.id.split('_tax:')[0]] = record.id.split('tax:d__')[1]
        
MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Archaea_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
rep_list = []
with open(MMseq, 'r') as clusters:
    for cluster in clusters:
        rep = cluster.split('\t')[0]
        rep_list.append(rep)
    clusters.close()        

from collections import Counter

rep_list = (rep_list)
reps = Counter(rep_list)

singletonlist = []
multiplelist = []
items = reps.items()

for i in items:
    if i[-1] == 1:
        singletonlist.append(i[0])
    else:
        multiplelist.append(i[0])

print(len(singletonlist))
print(len(multiplelist))

replist = []
Rep_dict = {}
with open(tsv_cluster, 'r') as clustermap:
    count = 0
    second_list = []
    element_count = 0
    for entry in clustermap:
        count = blank[entry.split('\t')[0]]
        repm = entry.split('\t')[0]
        if repm in multiplelist:
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
        #latest element to check in, drop count = 1 line as in multiple list does the same
import json
def get_taxonomy_cluster(cluster_rep):
    
    entry = Rep_dict[cluster_rep]
    num_proteins = (len(entry))
    species_list = []
    genus_list = []
    family_list = []
    order_list = []
    class_list = []
    phylum_list = []
    for prot_id in entry:
        species_list.append(Tax_dict[prot_id][-1].split(';s__')[1])
        genus_list.append(Tax_dict[prot_id][-1].split(';g__')[1].split(';s')[0])
        family_list.append(Tax_dict[prot_id][-1].split(';f__')[1].split(';g')[0])
        order_list.append(Tax_dict[prot_id][-1].split(';o__')[1].split(';f')[0])
        class_list.append(Tax_dict[prot_id][-1].split(';c__')[1].split(';o')[0])
        phylum_list.append(Tax_dict[prot_id][-1].split(';p__')[1].split(';c')[0])
    species_list = list(set(species_list))
    genus_list = list(set(genus_list))
    family_list = list(set(family_list))
    order_list = list(set(order_list))
    class_list = list(set(class_list))
    phylum_list = list(set(phylum_list))

    ls = len(species_list)
    lg = len(genus_list)
    lf = len(family_list)
    lo = len(order_list)
    lc = len(class_list)
    lp = len(phylum_list)
    string = 'num proteins: {}, numb uniq s: {}, num uniq g: {},  number of uniq f: {},  num uniq o: {}, num uniq c: {}, num uniq p {}'.format(num_proteins, ls, lg, lf, lo, lc, lp)
    return(string)



nclusters = 0
nproteins =  0
phyla_list = []
for cluster_rep in Rep_dict.keys():
    clusteroid = (get_taxonomy_cluster(cluster_rep))
    if clusteroid[-5] > 1:
        if clusteroid[-4] > 1:
            pass
        else:
            
            nclusters = nclusters + 1
    
            print('\n' + 'taxa of cluster seqs with more than one order, with rep {}'.format(cluster_rep))
            nproteins = nproteins + clusteroid[0] 
            entry = Rep_dict[cluster_rep]
            for prot_id in entry:
               print(Tax_dict[prot_id])
               phyla_list.append(Tax_dict[prot_id].split('f__')[1].split(';g')[0])
                
print(nclusters)
print(nproteins)
print(Counter(phyla_list))
