#this is a python script to test clade set up
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import Phylo
from Bio.Phylo.Consensus import *


#redo the clade name as it was done for the tree in previous code so we are uniform and consistent

        
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_7Nov_ft_lg_cat_gamma.nw', 'newick')
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
short_dict = {}
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        if len(clade.name) > 37:
            short = clade.name[0:37] + '...'
            short_dict[clade.name] = short
            pass
    else:
        clade.name = str(count)
        count = count + 1
        

#get the reps and their matching rep value

Red_dict = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/info_RED_clades_CPA_phylogeny_FTtree_lg_cat_gamma_7Nov.tsv', 'r') as clade_inf:
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
def get_subset_inlcuding_terminal(tree,lowerlimit,upperlimit, clade_2_exclude):
    parent_list = []
    skipped_list = []
    for clade in tree.find_clades():
        if clade != tree.root:
            if str(clade.name) != str(clade_2_exclude):
                #perform actions for terminal nodes
                if clade.is_terminal() == True:
                   #in this case the clade is terminal
    
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
    
    
                else:
                    if clade.is_preterminal() == True:
                        if len(clade.get_terminals()) == 2:
                            length_leaves = 0
                            for leaves in clade.get_terminals():
                                length_leaves += leaves.branch_length
                            if length_leaves == 0:
                                continue
                
                    if clade.RED > lowerlimit:
                        if clade.RED < upperlimit:
                                # we keep running into an issue were clade Bac226_18731726 keeps being labbeled as a non terminal node, and it has a branch_length of 0.38 while having a RED value of 1.
                                #something which
                                if type(clade.confidence) != float:
                                    continue
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
                        
                            
            else: 
                pass

    return(parent_list)

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.73, 1,'None')
print('Subgroup')
node_dic4 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic4[str(terminus)] = entry.name
print(len(node_dic4.keys()))
Subgroup_list = []
Subgroup_list =  parent_list
print(len(Subgroup_list))


parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.479, 0.729, 'None')
print('Subclade')
node_dic2 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic2[str(terminus)] = entry.name
print(len(node_dic2.keys()))
Subclade_list = []
Subclade_list =  parent_list
print(len(Subclade_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.14, 0.478, 'None')
print('Clade')
node_dic1 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic1[str(terminus)] = entry.name

print(len(node_dic1.keys()))
Clade_list = []
Clade_list =  parent_list
print(len(Clade_list))


parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.001, 0.14, 'None')
print('Subfamily')
node_dic3 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic3[str(terminus)] = entry.name
print(len(node_dic3.keys()))
Subfamily_list = []
Subfamily_list =  parent_list
print(len(Subfamily_list))


parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.001, 0.31, 27958)
print('Family')
node_dic8 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic8[str(terminus)] = entry.name
print(len(node_dic8.keys()))
Family_list = []
Family_list =  parent_list
print(len(Family_list))


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Red_informed_clade_20_lg_cat_gamma_RED_interval.tsv', 'w') as Big:
    header = 'Protein_id' + '\t' + 'Family' + '\t' + 'Subfamily' + '\t' + 'Clade' + '\t' + 'Subclade' + '\t' + 'Subgroup' + '\n'
    
    Big.write(header)

    for clade in tree.find_clades():
        if clade != tree.root:
            if clade.is_terminal() == True:
                if 'Homo' in clade.name:
                    clade.temp = clade.name
                elif 'Mus' in clade.name:
                    clade.temp = clade.name
                elif len(clade.name) > 37:
                    clade.temp = short_dict[clade.name]

                else:
                    clade.temp = clade.name
                if str(clade.temp) in node_dic1.keys():
                    Clade = node_dic1[str(clade.temp)]
                else:
                    Clade = 'RR'
                if str(clade.temp) in node_dic2.keys():
                    Subclade = node_dic2[str(clade.temp)]
                else:
                    Subclade = 'RR'
                if str(clade.temp) in node_dic3.keys():
                    Subfamily = node_dic3[str(clade.temp)]
                else:
                    Subfamily = 'RR'
                if str(clade.temp) in node_dic8.keys():
                    Family = node_dic8[str(clade.temp)]
                else:
                    Family = 'RR'
                if str(clade.temp) in node_dic4.keys():
                    Subgroup = node_dic4[str(clade.temp)]
                else:
                    Subgroup = 'RR'
               
                protein_id = clade.name

                Line = protein_id + '\t' +  Family + '\t' + Subfamily + '\t' + Clade + '\t' + Subclade + '\t' + Subgroup + '\n'
                Big.write(Line)
