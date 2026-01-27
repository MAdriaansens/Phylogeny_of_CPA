#part 1
#use this script to identify high confidence nodes were two taxa merge together 

from Bio import Phylo
from Bio.Phylo.Consensus import *
import json

#part 3
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_7Nov_ft_lg_cat_gamma.nw', 'newick')
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.comment = clade.name
        clade.bootstrap = clade.confidence
        clade.name = str(count)
        count = count + 1
    clade.combi = 'empty'
        
#

#give all the clades and preterminal clades a combi name
for clade in tree.get_terminals():
    if 'Bac' in clade.name:
        clade.combi = 'B'
    elif 'Arc' in clade.name:
        clade.combi = 'A'
    elif 'Euk' in clade.name:
        clade.combi = 'E'
    else:
        clade.combi = 'E'

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/parent_CPA_tree.json', 'r') as Out:
    Good_parents = json.load(Out)

clade_dict = {}
for clade in tree.find_clades():
    clade_dict[str(clade.name)] = clade

for clade in tree.find_clades(order='postorder'):
    origin_list = []
    if clade.combi != 'empty':
        pass
    else:
        children = (Good_parents[clade.name])
        for child_id in children:
            origin_list.append(clade_dict[child_id].combi)
            
        if origin_list.count('B') == len(origin_list):
            clade.combi = 'B'
        elif origin_list.count('A') == len(origin_list):
            clade.combi = 'A'
        elif origin_list.count('E') == len(origin_list):
            clade.combi = 'E'
        elif origin_list.count('AB') == len(origin_list):
            clade.combi = 'AB'
        else:
            origin_list.sort()
            origin_set = ''.join(origin_list)
            if 'AE' == origin_set:
                clade.comni = 'AE'
            elif 'AB' == origin_set:
                clade.combi = 'AB'
            elif 'BE' == origin_set:
                clade.combi = 'BE'
            

        if origin_list == ['B', 'BE']:
            clade.combi = 'B'
        elif origin_list == ['A', 'AB']:
            clade.combi = 'A'
        elif origin_list == ['AB', 'B']:
            clade.combi = 'B'
        elif origin_list == ['BE', 'E']:
            clade.combi = 'E'
        elif origin_list == ['BE', 'B']:
            clade.combi = 'B'
        elif origin_list == ['AB', 'BE']:
            clade.combi = 'B'
    if clade.combi == 'empty':
        print(clade.name)

Red_dict = {}
#this file contains all distances within a tree file
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/clades_CPA_phylogeny_FTtree_15sept.tsv', 'r') as clade_inf:
    next(clade_inf, None)
    for inf in clade_inf:
        RED =(inf.split('\t')[5])
        clade_id = inf.split('\t')[0]
        Red_dict[clade_id] = RED
print(len(Red_dict.keys()))
#overfamily

#match red values per tree
for clade in tree.find_clades():
    entry = str(clade.name)
    clade.RED = float(Red_dict[entry])
