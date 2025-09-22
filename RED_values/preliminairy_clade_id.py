from Bio import Phylo
from Bio.Phylo.Consensus import *



rep_dict = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/representatives_per_clades.tsv', 'r') as reps:
    next(reps, None)
    for rep in reps:
        rep_id = (rep.split('\t')[0])
        no_reps = rep.split('\t')[1].split('\n')[0]
        rep_dict[rep_id] = int(no_reps)
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/edit_trees/CPA_TREE_ALGINEDMMSEQ1_PF0099_28juli_fasttree_midrooted_REDNORM_final.nw', 'newick')
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.name = str(count)
        count = count + 1


Red_dict = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/clades_CPA_phylogeny_FTtree_15sept.tsv', 'r') as clade_inf:
    next(clade_inf, None)
    for inf in clade_inf:
        RED =(inf.split('\t')[5])
        clade_id = inf.split('\t')[0]
        Red_dict[clade_id] = RED
print(len(Red_dict.keys()))
#overfamily

for clade in tree.find_clades():
    entry = str(clade.name)
    clade.RED = float(Red_dict[entry])



def get_subset(tree,lowerlimit,upperlimit):
    parent_list = []
    for clade in tree.find_clades():
        if clade != tree.root:
            if clade.is_terminal() == False:
                if clade.RED > lowerlimit:
                    if clade.RED < upperlimit:
                        if clade.confidence > 0.94:
                            if clade.name not in parent_list:
                                if len(parent_list) == 0:
                                    parent_list.append(clade)
        
                                else:
                                    path_list = tree.get_path(clade)
                                    #this flag allows me to see if I am not redoing trees which already contain certain clades. 
                                    path_already_done = False
                    
                                    for x in parent_list:
                                        if x in path_list:
                                            path_already_done = True
                                    if path_already_done == True:
                                        pass
                                    else:
                                        parent_list.append(clade)

    return(parent_list)

paren_list =[]
parent_list = get_subset(tree, 0.042, 0.21)
zero = list(parent_list[0].get_terminals())
print(len(zero))
one = list(parent_list[1].get_terminals())
print(len(one))

print(parent_list) 
print(len(parent_list))
#first are 0.1811 and 0.0408 in RED
#size removal of clades smaller than 10 representatives
parent_list = []
parent_list = get_subset(tree, 0.21, 0.39)
print(len(parent_list))
for i in parent_list:
    if rep_dict[i.name] < 1000:
        parent_list.remove(i)

print(parent_list)
print(len(parent_list))

mice =list(parent_list[0].get_terminals())
dog =list(parent_list[1].get_terminals())
ape =list(parent_list[2].get_terminals())
pig =list(parent_list[3].get_terminals())
sheep =list(parent_list[4].get_terminals())
rat =list(parent_list[5].get_terminals())
horse =list(parent_list[6].get_terminals())
fly =list(parent_list[7].get_terminals())

print(len(fly))
print(len(horse))
print(len(pig))
print(len(rat))
print(len(sheep))
print(len(ape))
print(len(dog))
print(len(mice))

for entry in parent_list:
    termini_list = entry.get_terminals()
    for terminus in termini_list:

        node_dic[str(terminus)] = entry.name
print(len(node_dic.keys()))
parent_list = []
parent_list = get_subset(tree, 0.39, 0.56)
print(len(parent_list))
for i in parent_list:
    if rep_dict[i.name] < 1000:
        parent_list.remove(i)
print(len(parent_list))
print(parent_list)

node_dic = {}

for entry in parent_list:
    termini_list = entry.get_terminals()
    for terminus in termini_list:

        node_dic[str(terminus)] = entry.name
print(len(node_dic.keys()))

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Red_informed_clade_bigfamily.tsv', 'w') as Big:
    header = 'protein_id' + '\t' + 'big family0.042-0.21' + '\t' + 'subfamily0.21-0.39' + '\t' + 'small_family0.39-0.56' + '\n'
    Big.write(header)
    for clade in tree.find_clades():
        if clade != tree.root:
            if clade.is_terminal() == True:
                if clade in zero:
                    big_family = 'zero'
                elif clade in one:
                    big_family = 'one'
                else: 
                    big_family = 'R'
                #subfamily
                if clade in mice:
                    subfamily = 'mice'
                elif clade in dog:
                    subfamily = 'dog'

                elif clade in ape:
                    subfamily = 'ape'

                elif clade in pig:
                    subfamily = 'pig'

                elif clade in sheep:
                    subfamily = 'sheep'

                elif clade in rat:
                    subfamily = 'rat'

                elif clade in horse:
                    subfamily = 'horse'

                elif clade in fly:
                    subfamily = 'fly'


                else:
                    subfamily = 'R'

                if str(clade.name) in node_dic.keys():
                    small_family = node_dic[clade.name]
                else:
                    small_family = 'R'
                    
                protein_id = clade.name
                Line = protein_id + '\t' + big_family + '\t' + subfamily + '\t' + small_family + '\n'
                Big.write(Line)
