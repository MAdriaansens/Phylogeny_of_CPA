#this is a python script to test clade set up
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import Phylo
from Bio.Phylo.Consensus import *


#redo the clade name as it was done for the tree in previous code so we are uniform and consistent

#loas red normalized tree        
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_7Nov_ft_lg_cat_gamma.nw', 'newick')
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
short_dict = {}
#short is done because biopython is wierd with long names and messess it up sometimes
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
#since I am bad at decorating trees I have a seperate tsv file with for each node/termini the number of reps and their red in the lg_cat_gamma tree
#that is what lg_cat_gamma.tsv is, I open that file, matches the clade id with the red value in a dict so it is easy to parse
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/info_RED_clades_CPA_phylogeny_FTtree_lg_cat_gamma_7Nov.tsv', 'r') as clade_inf:
    next(clade_inf, None)
    for inf in clade_inf:
        RED =(inf.split('\t')[5])
        clade_id = inf.split('\t')[0]
        Red_dict[clade_id] = RED
#overfamily

#match red values per tree by making new trait to a given node
for clade in tree.find_clades():
    entry = str(clade.name)
    clade.RED = float(Red_dict[entry])

    
def get_subset_inlcuding_terminal(tree,lowerlimit,upperlimit, clade_2_exclude):

    #the parent list is a list of clade_ids which are well supported and not within any other clades. 
    parent_list = []
    bad_parent_list = []
    skipped_list = []
    for clade in tree.find_clades():

        #root is not of use here
        if clade != tree.root:
            if str(clade.name) != str(clade_2_exclude):
                #perform actions for terminal nodes
                if clade.is_terminal() == True:
                    pass
                else:     
                    #all clades are currently non-terminal so they are split/bifurcations
                    #we check if the clades fall within RED window
                    if clade.RED < lowerlimit:
                        pass
                    else:
                        if clade.RED > upperlimit:
                            pass
                        else:
                        # we keep running into an issue were clade Bac226_18731726 keeps being labbeled as a non terminal node, and it has a branch_length of 0.38 while having a RED value of 1.
                                if type(clade.confidence) != float:
                                    continue
                                else:
                                    pass
                                #now it is time to check for anything equal or greater than 95
                                if clade.confidence > 0.94:
                                    if clade.name not in parent_list:
                                        if len(parent_list) == 0:
                                            parent_list.append(clade)
                                            path_already_done = False
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
                                    #if they fall within the red winow and have <95 bootstrap we exclude them, we also check if they are already part of a good clade
                                    path_list = tree.get_path(clade)
                                            #this flag allows me to see if I am not redoing trees which already contain certain clades. 
                                    path_already_done = False
                            
                                    for x in parent_list:
                                        if x in path_list:
                                            path_already_done = True
                                    if path_already_done == True:
                                        pass
                                    else:
                                        bad_parent_list.append(clade)
                            


    return(parent_list, bad_parent_list)
parent_list =[]
well_supported_clade_list = [] 


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def is_semipreterminal(clade):
    #"""True if any direct descendent is terminal."""
    for child in clade:
        if child.is_terminal():
            return True
    return False


#this bit of code checks if a bad parent results in a good bifurcations down the road
def parse_parents_list(tree, parent_list, bad_parent_list, lowerlimit, upperlimit, to_exclude):
    missing_parent_list = []
    bad_parent_good_child = []
    node_dic = {}
    missed = 0
    parents_beyond_upperlimit = {}
    #all children and the parents themself of the parent list are of limits, since they are good
    off_limit = []

    #we take all the clades and children of the good parent to make sure they will not be included
    for clade in tree.find_clades():
        if clade.name == to_exclude:
            pass
        if clade in parent_list:
            off_limit.append(clade)
            for child in clade.find_clades():
                off_limit.append(child)
    print(len(off_limit))
    #now i run over the list and check if any of the parents contain any extant clades which are present in 'good' parent list
    #if they have good extant taxa they are placed into bad_parent_good_child list (meaning that a split results in a well supported clade and one poorly supported one)
    #if they only have bad extant taxa within the border set they only contain poorly supported taxa and are classified as poorly supported and kept in bad_parent_list
    #because bad parents with bad children are bad

    
    for clade in tree.find_clades():
        if clade.name == to_exclude:
            pass
        if clade in bad_parent_list:
            for parent in (tree.get_path(clade)):
                #here we check if some of the children from a bad clade are well supported within the red interval
                if parent in parent_list:
                    bad_parent_list.remove(parent)
                    #these nodes result in well supported children
                    bad_parent_good_child.append(parent)

       
    for clade in tree.find_clades():
        if clade.name == to_exclude:
            pass
        #clade.clades returns the first extant clade from a clade called, which could be a bad or a good clade
        if clade in bad_parent_good_child:
            if clade  != tree.root:
                for extant_clade in clade.find_clades():
                    if extant_clade in off_limit:
                        pass
                        #this means we do nothing with the good clades in the bifurcation and only perform actions on the bad clades in the bifurcation
                    elif extant_clade.is_terminal():
                        pass
                    else:
                        bad_parent_list.append(extant_clade)
                        if extant_clade in bad_parent_good_child:
                            bad_parent_good_child.remove(extant_clade)
            
    #these nodes have no well supproted (good) bifurcations as children within the given range, the bad parent bad children will be named pb_: (poor bootstrap)
    bad_parent_bad_children =  bad_parent_list
    #I want to idenitfy the termini which will become RR in the future
    missed_out_list = []

        #this identifies things within the red interval with poor bootstrap support bifurcations
    for entry in bad_parent_bad_children:
        if entry.name == to_exclude:
            pass
        if isinstance(entry, int) == True:

            pass
        if entry in off_limit:
            pass
        else:
            if entry.is_terminal() == False:
                if is_semipreterminal(entry) == True:
                    pass
                else:
                    termini_list = entry.get_terminals()
                    #I want to seperate out which clades are well supported
                    for terminus in termini_list:
            
                        node_dic[str(terminus)] = 'pb_' +  entry.name 


    for entry in parent_list:
        if entry.name == to_exclude:
            pass
        
        if entry.is_terminal() == False:
            termini_list = entry.get_terminals()
            #I want to seperate out which clades are well supported
            for terminus in termini_list:
    
                node_dic[str(terminus)] = entry.name
    #shows me which termini are left out
    for clade in tree.find_clades():
        if clade.name == to_exclude:
            pass
        #seed sequences have issues with their naming and might freak out the system so those are skipped
        if clade != tree.root:
            if clade.is_terminal():
                #this is to check the representative sequences
                if clade.name[3] == '_':
                    pass
                else:
                    if clade.name not in node_dic.keys():
                        #now we check if the reason why we miss a termini is because its split is before the lowerlimit of the RED interval
                        #if it does then the termini will just be kept seperate since the split occured before the lowerlimit and since they entered
                        #the interval as already seperate termini
    
    
                        #rework, get red of parent directly after the split happend (why did I forget to add this??)
                        if get_parent(tree, clade).RED < lowerlimit:
                            node_dic[str(clade.name)] = 'tbl_' + clade.name
                        elif get_parent(tree, clade).RED > upperlimit:
                            #now if the split occurs after red upperlimit we see if the the clade is terminal or
                            #if the split contains multiple
                            if len(get_parent(tree, clade).get_terminals()) > 1:
                                if is_semipreterminal(get_parent(tree, clade)) == True:
                                    node_dic[str(clade.name)] = 'tal_' + clade.name
                                else:
                                    node_dic[str(clade.name)] = 'nal_' + get_parent(tree, clade).name
                            else:
                                node_dic[str(clade.name)] = 'tal_' + clade.name 
                                missing_parent_list.append(get_parent(tree, clade).name)
                        else:
                            #it falls in between both

                            if clade.is_terminal() == True:
                                #poor bootstrap terminal
                                node_dic[str(clade.name)] = 'pbt_' + clade.name
                            else:
                                node_dic[str(clade.name)] = 'in between'
                    else:
                        pass


                        
    #checks
    if len(node_dic) != 66978:
        if len(bad_parent_good_child) == 0:
            print('missing branches are already cut off by lower threshold')
        else:
            print('soo bad parents and good children')
            print(len(bad_parent_good_child))
    print(missed, 'missed')
    Counter(missing_parent_list)
    return(node_dic, bad_parent_good_child)

#generate dictionary for each protein group
#Clade
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.001, 0.14, 'None')
node_dic1, bad_parent_good_child = parse_parents_list(tree, parent_list, bad_parent_list, 0.001, 0.14, 'None')

#Subfamily
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.001, 0.31, 27958)
node_dic2, bad_parent_good_child= parse_parents_list(tree, parent_list, bad_parent_list, 0.001, 0.31, 27958)

#Subclade
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.14, 0.478, 'None')
node_dic3, bad_parent_good_child = parse_parents_list(tree, parent_list, bad_parent_list, 0.14, 0.478, 'None')

#Group
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.479, 0.729, 'None')
node_dic4, bad_parent_good_child = parse_parents_list(tree, parent_list, bad_parent_list, 0.479, 0.729,'None')

#Subgroup
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.73, 1,'None')
node_dic5, bad_parent_good_child = parse_parents_list(tree, parent_list, bad_parent_list, 0.73, 1,'None')

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Red_informed_clade_20_lg_cat_gamma_RED_interval_15Jan.tsv', 'w') as Big:
    header = 'Protein_id' + '\t' + 'Subfamily' + '\t' + 'Clade' + '\t' + 'Subclade' + '\t' + 'Group' + '\t' + 'Subgroup' + '\n'
    
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
                    Subfamily = node_dic2[str(clade.temp)]
                else:
                    Subfamily = 'RR'
                if str(clade.temp) in node_dic3.keys():
                    Subclade = node_dic3[str(clade.temp)]
                else:
                    Subclade = 'RR'
                if str(clade.temp) in node_dic4.keys():
                    Group = node_dic4[str(clade.temp)]
                else:
                    Group = 'RR'
                if str(clade.temp) in node_dic5.keys():
                    Subgroup = node_dic5[str(clade.temp)]
                else:
                    Subgroup = 'RR'
               
                protein_id = clade.name

                Line = protein_id + '\t' + Subfamily + '\t' + Clade + '\t' + Subclade + '\t' + Group + '\t' + Subgroup + '\n'
                Big.write(Line)
