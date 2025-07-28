from Bio import Phylo
from Bio.Phylo.Consensus import *
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/ITOL_TREE_ALGINEDMMSEQ70_PF0099_28juli_fasttree.treefile', 'newick')

#first generate a midpoint rooted tree.
tree.root_at_midpoint()
#save midpoint rooted tree so this time expensive step does not need to be repeated
Phylo.write(tree, "./midrooted_ITOL_TREE_ALIGNED70_28JULI.nw", "newick")
count = 0


#decorate and identify all clades/nodes in a tree
for clade in tree.find_clades():
   target_clade = clade
   target_clade.name = str(count)
   count = count + 1
#if a a node has only terminal children then it will return false on is_semipreterminal(). 
#i.e. if a tree clade has two termini as children then it returns false
clade_depths = tree.depths()
def is_semipreterminal(clade):
    """True if any direct descendent is terminal."""
    for child in clade:
        if child.is_terminal():
            return True
    return False

#passed_thres_count =  0

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/clades_CPA_phylogeny.tsv', 'a') as C:
    header = 'clade_id' + '\t' + 'distance_from_root' + '\t' + 'branch_length' + '\t' + 'bootstrap' + '\t' + 'children' + '\t' + 'bootstrap_value' +  '\n'
    C.write(header)
    bootstrap_value = ''
    distance_from_root = ''
    for clade in tree.find_clades():
        #We are not interested in termini
        if clade.is_terminal() == True:
            pass
        else:
            #this is to ensure that only nodes with confidence/bootstrap and branchlength are taken
            if clade.confidence is None:
                pass
            elif clade.branch_length is None:
                pass
            else:
                #set a confidence threshold
                if clade.confidence > 0.949:
                    #what is_semipreterminal checks is if it is a clade/node for 2 termini
                    if (is_semipreterminal(clade)) == True:
                        pass
                    else:
                       #we take the distance by taking the tree.root
                        distance_from_root = tree.distance(tree.root, '{}'.format(clade.name))
                        #print(clade, clade.confidence, clade.branch_length, distance_from_root)
                        children = (len(clade.get_terminals()))
                        passed_thres_count = passed_thres_count + 1
                        bootstrap_value  = 'good'
                else:
                    #if you want we can include all clades which have a below 0.95 bootstrap
                    pass
                     #if (is_semipreterminal(clade)) == True:
                     #   pass
                    # else:
                    #    distance_from_root = tree.distance(tree.root, '{}'.format(clade.name))
                        #print(clade, clade.confidence, clade.branch_length, distance_from_root)
                    #    children = (len(clade.get_terminals()))
                        #Phylo.draw_ascii(clade)
                     #   passed_thres_count = passed_thres_count + 1
                  #      bootstrap_value  = 'bad'
                line =  clade.name + '\t' + str(tree.distance(tree.root, '{}'.format(clade.name))) + '\t' + str(clade.branch_length) + '\t' + str(clade.confidence) + '\t' + str(len(clade.get_terminals())) + '\t' + bootstrap_value +  '\n'
                C.write(line)
