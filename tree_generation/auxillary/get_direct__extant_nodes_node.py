#part 1
from Bio import Phylo
from Bio.Phylo.Consensus import *
import json

#part 3
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_12sept_ft.nw', 'newick')
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

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    if len(node_path) == 1:
        return (tree.root)
    else:
        return node_path[-2]

#returns the parent of each clade
parent_dic = {}
good_parent_dic = {}
for clade in tree.find_clades():
    if clade != tree.root:
        parent = str(get_parent(tree, clade))

        if parent in parent_dic:
            children = (str(clade.name), parent_dic[parent])
            good_parent_dic[parent]= children
        else:
            parent_dic[parent] = str(clade.name)
with open('parent_CPA_tree.json', 'w') as Out:
    json.dump(good_parent_dic, Out)
