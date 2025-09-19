
from Bio import Phylo
from Bio.Phylo.Consensus import *

def Get_representatives(MMseq):
    rep_list = []
    with open(MMseq, 'r') as clusters:
        for cluster in clusters:
            rep = cluster.split('\t')[0]
            rep_list.append(rep)
        clusters.close()
    
    from collections import Counter
    
    #the rep list, is just a list of all values in the representative columns, this means that ids are present multiple times
    
    rep_list = (rep_list)
    reps = Counter(rep_list)
    
    #reps returns is a counter of how many sequences a rep represents
    items = reps.items()
    
    
    singleton_list = []
    multiple_list = []
    
    for i in items:
        if i[-1] == 1:

            singleton_list.append(i[0])
        else:
            multiple_list.append(i[0])
    
    print(len(multiple_list))
    print(len(singleton_list))
    
    replist = []
    Rep_dict = {}
    
    with open(MMseq, 'r') as clustermap:
        count = 0
        second_list = []
        element_count = 0
        for entry in clustermap:
            count = reps[entry.split('\t')[0]]
            repm = entry.split('\t')[0]
    
            if repm in multiple_list:
                if repm in replist:
                    second = (entry.split('\t')[1].split('\n')[0])
                    second_list.append(second)
                    if len(second_list) == count:
                        Rep_dict[repm] = second_list
                else:
    
                    replist.append(repm)
                    second_list = []
                    second_list.append(repm)
            else:
                pass
    return(Rep_dict)

MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Euk_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
Erepdic = Get_representatives(MMseq)
MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Archaea_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
Arepdic = Get_representatives(MMseq)
MMseq= '/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/MMseq/Bacteria_allhmmscanned_final_clustered_at0.7.fasta_cluster.tsv'
Brepdic = Get_representatives(MMseq)
Seq_taxa = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_full_tree_aligned_annotation_13sept.tsv', 'r') as A:
    next(A, None)
    for line in A:
        Seq_taxa[(line.split('\t')[0])] =(line.split('\t')[5]) 
        
def get_clusteroids_per_clade(list_termini):
    total_clusteroids = 0

    for entry in list_termini:
        entry = str(entry)
        origin = entry.split('_')[0]
        if origin == 'Bac226':
            if entry not in Brepdic.keys():
                reps=1
            else:
                reps = (len(Brepdic[entry]))
        elif origin == 'Arc226':
            if entry not in Arepdic.keys():
                reps=1
            else:
                reps = (len(Arepdic[entry]))
        elif origin == 'EukM6':
            if entry not in Erepdic.keys():
                reps=1
            else:
                reps = (len(Erepdic[entry]))
        else:
            #representatives
            reps = 1
        
        total_clusteroids += reps
    return(total_clusteroids)
def enough_proteins_in_clade(i):
    non_euk_count = 0
    enough = 'no'
    if i.reps > 67254:
        enough = 'too_much'
    else:
        for term in i.get_terminals():
            if str(term).split('_')[0] != 'EukM6':
                                    #make sure it is not a seed
                if len(str(term).split('_')) == 2:
                    non_euk_count +=1
    return(non_euk_count)
    return(total_clusteroids)

tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/edit_trees/CPA_TREE_ALGINEDMMSEQ1_PF0099_28juli_fasttree_midrooted.nw', 'newick')
count = 0
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.name = str(count)
        count = count + 1
        
for clade in tree.find_clades(order='postorder'):
    list_termini = []
    representatives_in_clade = 0
    if clade == clade.root:
        list_termini = clade.get_terminals()
        representatives_in_clade = get_clusteroids_per_clade(list_termini)
        clade.reps = representatives_in_clade
#generate mini trees

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/representatives_per_clades.tsv', 'w') as Rep_clades:
    Header = 'clade_id' + '\t' + 'no_reps' + '\n'
    Rep_clades.write(Header)
    for clade in tree.find_clades(order='postorder'):
        clade_id = clade.name
        no_reps = clade.reps
        Line = clade_id + '\t' + str(no_reps) + '\n'
        Rep_clades.write(Line)

  import os
skip_list = []
seq_skip_list = []

for clade in tree.find_clades(order='postorder'):
    if clade.is_terminal() != True:
        pass
    else:
        check = False
        origin = str(clade.name)
        if origin.split('_')[0] == 'EukM6':
            if origin in seq_skip_list:
                pass
            else:
                if clade.reps < 5:
                    pass
                else:
                    #does not draw terminals
                    path_list = tree.get_path(clade)
                    #this flag allows me to see if I am not redoing trees which already contain certain clades. 
                    path_already_done = False
    
                    for x in skip_list:
                        if x in path_list:
                            path_already_done = True
                    if path_already_done == True:
                        pass
                    else:
                        check = False
                        for i in path_list:
                            if origin in seq_skip_list:
                                pass
                            else:
                                if i.reps > 67254:
                                    pass
                                else:
                                    no_reps_in_clade = int(enough_proteins_in_clade(i))
                                    if no_reps_in_clade < 2000:
                                        if no_reps_in_clade > 70 :
                                                if os.path.exists('terminal_euk_prot_{}_no_proteins_total_{}.nwk'.format(origin, no_reps_in_clade)):
                                                    break
                                                else:
                                                    print(i, origin, no_reps_in_clade)
        
                                                    Phylo.write(i, 'terminal_euk_prot_{}_no_proteins_total_{}.nwk'.format(origin, no_reps_in_clade), 'newick')
                                                    skip_list.append(i)
                                                    seq_skip_list.append(origin)
