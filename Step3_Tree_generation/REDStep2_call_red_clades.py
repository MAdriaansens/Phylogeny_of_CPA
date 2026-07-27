#this is a python script to test clade set up
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import Phylo
from Bio.Phylo.Consensus import *

from collections import Counter

#redo the clade name as it was done for the tree in previous code so we are uniform and consistent

#loas red normalized tree        
tree = Phylo.read('~/RED/CPA_TREE_ALGINEDMMSEQ1_PF00999_28April_fasttree_midrooted.treefile', 'newick')
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
with open('~/RED/info_RED_clades_CPA_phylogeny_FTtree_28April.tsv', 'r') as clade_inf:
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

#def parse_parents_list(tree, parent_list, bad_parent_list, lowerlimit, upperlimit, to_exclude):

def assign_clades(tree, parent_list, bad_parent_list, lowerlimit, upperlimit, to_exclude):
    presence_good_step = 'no'
    node_dic = {}
    missed = 0
        #all children and the parents themself of the parent list are of limits, since they are good
    off_limit = []
    print(len(tree.get_terminals()))
        #we take all the clades and children of the good parent to make sure they will not be included
    for clade in tree.find_clades():
        if clade.name == to_exclude:
           pass
        if clade in parent_list:
           off_limit.append(clade)
           for child in clade.find_clades():
               off_limit.append(child)
    #this gives a value to all the clades and their termini if their latest common ancestor is well supported
    for entry in parent_list:
        if entry.name == to_exclude:
            pass
            
        if entry.is_terminal() == False:
           termini_list = entry.get_terminals()
                #I want to seperate out which clades are well supported
           for terminus in termini_list:
        
               node_dic[str(terminus)] = entry.name
    print(len(node_dic), 'node dic')
    
    
    to_sort_out = []
    
    #the issue is that this code assumes that all termini have nodes (be they well supported or not) within the red-interval,
    #if that is not the case, i.e it branches of before the red-interval it is not detected
    for bad_parent in bad_parent_list:
       for clade in bad_parent.find_clades():
           if clade in off_limit:
               pass
           else:
               if clade.is_terminal():
                   path = bad_parent.get_path(clade)
                   presence_good_step = 'no'
                   for step in path:
                       if step.is_terminal():
                               pass
                       elif step.RED < lowerlimit:
                           pass
                       elif type(step.confidence) != float:
                            if step.RED == 1.0:
                                pass
                       elif step.confidence > 0.94:
                           if presence_good_step == 'yes':
                               pass
                           else:
                               presence_good_step = 'yes'
                               #we know that it not part of the good parents,
                               node_dic[clade.name] = "lgs_{}".format(step.name)
                   if presence_good_step == 'yes':
                       pass
                   else:
                       node_dic[clade.name] = "t_{}".format(clade.name)
               else:
                   pass
    
    
    
    missing_termini_list=[]
    count=0
    print(len(node_dic), 'node dic 2')
    if len(node_dic) != len(tree.get_terminals()):
        print('missing some')
        for terminus in tree.get_terminals():
            if terminus.name not in node_dic:
                if terminus.name not in short_dict:
                    presence_good_step = 'no'
                    path = tree.root.get_path(terminus)
                    for step in path:
                        if step.is_terminal():
                            pass
                        elif step.RED < lowerlimit:
                            pass
                        elif step.confidence != float:
                            if step.RED==1.0:
                                pass
                        elif step.confidence > 0.94:
                            if presence_good_step == 'yes':
                                pass
                            else:
                                presence_good_step = 'yes'
                                node_dic[terminus.name] = "lgs_{}".format(step.name)
                    if presence_good_step == 'yes':
                        pass
                    else:
                        node_dic[terminus.name] = "t_{}".format(clade.name)
    
    print(count)
    print(len(node_dic))
    
    if len(node_dic) != len(tree.get_terminals()):
        print('still missing some')
    else:
        print('complete')
    return(node_dic)


parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree,  0.001, 0.11, 'None')
node_dic1 = assign_clades(tree, parent_list, bad_parent_list, 0.001, 0.11, 'None')

#Subfamily

parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.11, 0.478, 'None')
node_dic2= assign_clades(tree, parent_list, bad_parent_list, 0.11, 0.478, 'None')

#subclade

parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.478, 0.78, 'None')
node_dic3 = assign_clades(tree, parent_list, bad_parent_list, 0.478, 0.78, 'None')


#Group
parent_list =[]
parent_list, bad_parent_list = get_subset_inlcuding_terminal(tree, 0.78, 1, 'None')
node_dic4= assign_clades(tree, parent_list, bad_parent_list,  0.78, 1,'None')


with open('Red_informed_clade_20_lg_cat_gamma_RED_interval_28Juli.tsv', 'w') as Big:
    header = 'Protein_id' + '\t' + 'Family' + '\t' + 'Subfamily' + '\t' + 'Subclade' + '\t' + 'Group' + '\n'
    
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
                    Family = node_dic1[str(clade.temp)]
                else:
                    Family = 'RR'
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

               
                protein_id = clade.name

                Line = protein_id + '\t' + Family + '\t' + Subfamily + '\t' + Subclade + '\t' + Group + '\n'
                Big.write(Line)
