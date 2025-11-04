#this is a python script to test clade set up
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import Phylo
from Bio.Phylo.Consensus import *


#redo the clade name as it was done for the tree in previous code so we are uniform and consistent

        
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_7Oktsept_ft.nw', 'newick')
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
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_info_CPA_phylogeny_FTtree_7oct.tsv', 'r') as clade_inf:
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
def get_subset_inlcuding_terminal(tree,lowerlimit,upperlimit):
    parent_list = []
    skipped_list = []
    for clade in tree.find_clades():
        if clade != tree.root:

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
                    
                        


    return(parent_list)



parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.225, 0.5)
print('genus')
node_dic7 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic7[str(terminus)] = entry.name
print(len(node_dic7.keys()))

genus_list = []
genus = parent_list

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.24,0.4275)
print('family')
node_dic6 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic6[str(terminus)] = entry.name
print(len(node_dic6.keys()))
family_list = []
family_list =  parent_list
print(len(family_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.24, 0.47)
print('order')
node_dic5 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic5[str(terminus)] = entry.name
print(len(node_dic5.keys()))
order_list = []
order_list =  parent_list
print(len(order_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.225, 0.4275)
print('class')
node_dic4 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic4[str(terminus)] = entry.name
print(len(node_dic4.keys()))
class_list = []
class_list =  parent_list
print(len(class_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.1991,0.4275)
print('Phylum1')
node_dic3 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic3[str(terminus)] = entry.name
print(len(node_dic3.keys()))
Phylum1_list = []
Phylum1_list =  parent_list
print(len(Phylum1_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0,0.225)
print('Phylum2')
node_dic2 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic2[str(terminus)] = entry.name
print(len(node_dic2.keys()))
Phylum2_list = []
Phylum2_list =  parent_list
print(len(Phylum2_list))

parent_list =[]
parent_list = get_subset_inlcuding_terminal(tree, 0.0,0.24)
print('Domain')
node_dic1 = {}
print(len(parent_list))
for entry in parent_list:
    if entry.is_terminal() == False:
        termini_list = entry.get_terminals()
        for terminus in termini_list:

            node_dic1[str(terminus)] = entry.name

print(len(node_dic1.keys()))
Domain_list = []
Domain_list =  parent_list
print(len(Domain_list))

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Red_informed_clade_bigfamily_3nov.tsv', 'w') as Big:
    header = 'protein_id' + '\t' + 'Domain' + '\t' + 'Phyla1' + '\t' + 'Phyla2' + '\t' + 'Class' + '\t' + 'Order' + '\t' + 'Family' + '\t' + 'Genus' + '\n'
    Big.write(header)

    for clade in tree.find_clades():
        if clade != tree.root:
            if clade.is_terminal() == True:
                if len(clade.name) > 37:
                    clade.temp = short_dict[clade.name]
                else:
                    clade.temp = clade.name
                if str(clade.temp) in node_dic1.keys():
                    Domain = node_dic1[str(clade.temp)]
                else:
                    Domain = 'RR'
                if str(clade.temp) in node_dic2.keys():
                    Phyla1 = node_dic2[str(clade.temp)]
                else:
                    Phyla1 = 'RR'
                if str(clade.temp) in node_dic3.keys():
                    Phyla2 = node_dic3[str(clade.temp)]
                else:
                    Phyla2 = 'RR'
                if str(clade.temp) in node_dic4.keys():
                    Class = node_dic4[str(clade.temp)]
                else:
                    Class = 'RR'
                if str(clade.temp) in node_dic5.keys():
                    Order = node_dic5[str(clade.temp)]
                else:
                    Order = 'RR'
                if str(clade.temp) in node_dic6.keys():
                    Family = node_dic6[str(clade.temp)]
                else:
                    Family = 'RR'
                if str(clade.temp) in node_dic7.keys():
                    Genus = node_dic7[str(clade.temp)]
                else:
                    Genus = 'RR'
    
                protein_id = clade.name

                Line = protein_id + '\t' + Domain + '\t' + Phyla1 + '\t' + Phyla2 + '\t' + Class + '\t' + Order + '\t' + Family + '\t' + Genus + '\n'
                Big.write(Line)
