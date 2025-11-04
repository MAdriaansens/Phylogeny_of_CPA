from Bio import Phylo
from Bio.Phylo.Consensus import *
import json

clade_id = {}
#outputs midrooted tree, a tsv file of all clade info + midrooted RED normalized tree
import os
file_path = "/nesi/nobackup/uc04105/new_databases_May/final_tree_set/edit_trees/CPA_TREE_ALGINEDMMSEQ1_PF00999_3Nov_fasttree_midrooted.treefile"

if os.path.exists(file_path):
    print('yah')
    tree = Phylo.read("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/edit_trees/CPA_TREE_ALGINEDMMSEQ1_PF00999_3Nov_fasttree_midrooted.treefile", "newick")

else:
    tree = Phylo.read("/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_TREE_ALGINEDMMSEQ1_PF00999_3Nov_fasttree.treefile", "newick")
    tree.root_at_midpoint()
    
    Phylo.write(tree, "/nesi/nobackup/uc04105/new_databases_May/final_tree_set/edit_trees/CPA_TREE_ALGINEDMMSEQ1_PF00999_3Nov_fasttree_midrooted.nw", "newick")
#this has to occur after midpoint root
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.comment = clade.name
        clade.name = str(count)
        count = count + 1

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)

    return node_path
for clade in tree.find_clades(order='postorder'):
    avg_div = 0
    if clade.is_terminal() == True:
        clade.mean_dist = 0.0
    else:
        clade.num_taxa = clade.count_terminals()
        for c in clade.clades:
            avg_div +=(float(c.count_terminals())/clade.count_terminals())*(c.mean_dist + c.branch_length)
    clade.mean_dist = avg_div

for clade in tree.find_clades(order='preorder'):
        clade.og_branch_length = clade.branch_length
        clade.bootstrap = clade.confidence

        if clade == tree.root:
            clade.red = 0.0

        elif clade.is_terminal() == True:
            clade.red = 1.0
        else:
            a = clade.branch_length
            b = clade.mean_dist
            parent = get_parent(tree, clade)
            if len(parent) == 1:
               x =  0
            else:
                #the parent in this case is the root
                x = parent[-2].red
            if (a + b) != 0:
               red = x + (a / (a + b)) * (1.0 - x)
            else:
               red = x

            clade.red = red
print('red done')


import json

#important file when generating and calling clades later on
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/info_RED_clades_CPA_phylogeny_FTtree_3Nov.tsv', 'w') as C:
    header = 'clade_id' + '\t' + 'distance_from_root' + '\t' + 'og_branch_length' + '\t' + 'bootstrap' + '\t' + 'children' + '\t' + 'RED' + '\t' + 'summed_red' + '\t' + 'RED_BL' + '\n'
    C.write(header)
    for clade in tree.find_clades(order='postorder'):
       if clade == tree.root:
           clade.branch_length = 0
       else:
           parent = get_parent(tree, clade)
           if len(parent) == 1:
                #the parent in this case is the root
               x =  0
           else:
               x = parent[-2].red

           clade.branch_length = clade.red - x
       if clade.is_terminal() == False:
           summed_red = 0
           for node in tree.get_path(clade):
               summed_red += node.red
       else:
           summed_red = 0
           for node in tree.get_path(clade):
               if node.name != clade.name:
                   summed_red += node.red
               
       line = clade.name + '\t' + str(tree.distance('{}'.format(clade.name))) + '\t' + str(clade.og_branch_length) + '\t' + str(clade.confidence) + '\t' + str(len(clade.get_terminals())) + '\t' + str(clade.red) + '\t' +  str(summed_red) + '\t' + str(clade.branch_length) + '\n'
       C.write(line)

for clade in tree.find_clades(order='postorder'):
    if clade == tree.root:
         clade.branch_length = 0
    else:
        parent = get_parent(tree, clade)
        if len(parent) == 1:
             #the parent in this case is the root
            x =  0
        else:
            x = parent[-2].red

        clade.branch_length = clade.red - x

for clade in tree.find_clades():
    if clade.is_terminal() != True:
        clade.name = clade.comment
    clade.confidence = clade.bootstrap

Phylo.write(tree, "./CPA_phylogeny_midroot_TREE_ALIGNEDRED_3Nov_ft.nw", "newick")

