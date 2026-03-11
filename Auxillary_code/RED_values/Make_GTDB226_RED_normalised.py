from Bio import Phylo
from Bio.Phylo.Consensus import *

clade_id = {}

tree = Phylo.read("/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/226_tree/ar53_r226.tree", "newick")


#I still midroot the tree
tree.root_at_midpoint()
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.name = str(count)
        count = count + 1
#in Phylo there is not a clear and easy way to get the parent of a node, atleast as far as I could tell. 
#Because GTDBTK is written with a different phylogentic python software I had to create the function get_parent

#it basically terurns me  a list of the path it takes to get to the child clade 
#the last element in the list is information about the child clade itself, and the -2 is the information about the parent

#so for node 874 this would return node 871 as the parent

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)

    return node_path

#add clade.mean_dist attribute
#this part is copied from gtdbtk/relative_disance.py, it is important to note that GTDBTK does not use BioPython.Phylo but rather DendroPy
for clade in tree.find_clades(order='postorder'):
    avg_div = 0
    #check if clade is not terminal and has children
    if clade.is_terminal() == True:
        clade.mean_dist = 0.0
    else:
        #if it has children
        clade.num_taxa = clade.count_terminals()
        for c in clade.clades:
            avg_div +=(float(c.count_terminals())/clade.count_terminals())*(c.mean_dist + c.branch_length)
    clade.mean_dist = avg_div
#add clade.red attribute

for clade in tree.find_clades(order='preorder'):
  #make an attribute which keeps the original branchlength for when, in the future, I replace branchlength with RED
        clade.og_branch_length = clade.branch_length
  #RED of root is always 0
        if clade == tree.root:
            clade.red = 0.0
     
        elif clade.is_terminal() == True:
          #RED of a leaf or terminal is always 1
            clade.red = 1.0
        else:
            a = clade.branch_length
            b = clade.mean_dist
            parent = get_parent(tree, clade)
            if len(parent) == 1:
              #if len(parent) == 1 it is the root
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

#originally I had the branch length just be the RED of the tree but that leads to issues and does not make sense when you think about it. 
#So with the knowledge that the root RED was always 0 and the termini/leaf was laways 1 I could now start replacing the branch length with the normalized RED length.
#calculation of RED is done postoder, i.e. first the left most root completely, for more details see: https://www.geeksforgeeks.org/dsa/postorder-traversal-of-binary-tree/

#for further analysis purposes I also create a .tsv file of the tree
with open("/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/RED_decorated_trees/Ar53_GTDB226_phylo_midroot.tsv", 'w') as C:
    header = 'clade_id' + '\t' + 'distance_from_root' + '\t' + 'og_branch_length' + '\t' + 'bootstrap' + '\t' + 'children' + '\t' + 'RED' + '\t' + 'RED_BL' + '\n'
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
               #retrieves the red of the parent node. 
               x = parent[-2].red
            #we are working or starting at the root. So the first node after the root has a red of 0.012629 and the root has a red of 0
            #so 0.012629 - 0 equals 0.012629 and therefore the branchlength is 0.012629. You can check this by running "tree.common_ancestor('1')"
           clade.branch_length = clade.red - x
       line = clade.name + '\t' + str(tree.distance('{}'.format(clade.name))) + '\t' + str(clade.og_branch_length) + '\t' + str(clade.confidence) + '\t' + str(len(clade.get_terminals())) + '\t' + str(clade.red) + '\t' + str(clade.branch_length) + '\n' 
       C.write(line)
  
#the old root is now node 4686: (GB_GCA_024860865.1|GB_GCA_036482035.1), which in the new RED has an red of 0.2089
#i.e p__JACRDV01 g__JACRDV01_sp024860865 and p__Hydrothermarchaeota g__JAHJSX01_sp036482035

#I do not decorate the tree or save it in such a way that it retains all the attributes, this is something I still will have to do. 
Phylo.write(tree, "/nesi/nobackup/uc04105/new_databases_May/GTDB_226/GTDBK/RED_decorated_trees/RED_A53_ALIGNEDRED_30julift.nw", "newick")
