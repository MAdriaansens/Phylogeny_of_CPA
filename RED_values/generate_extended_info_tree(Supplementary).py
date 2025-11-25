from Bio import Phylo
header = 'clade_id' + '\t' + 'distance_from_root' + '\t' + 'og_branch_length' + '\t' + 'bootstrap' + '\t' + 'children' + '\t' + 'RED' + '\t' + 'summed_red' + '\t' + 'RED_BL' +  '\t' + 'Termini_list' + '\t' + 'daughter_nodes' + '\n'





Info_dict = {}
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/info_RED_clades_CPA_phylogeny_FTtree_lg_cat_gamma_7Nov.tsv', 'r') as info:
    next(info, None)
    for line in info:
        if (int(line.split('\t')[4])) < 2:
            pass
        else:
            Info_dict[line.split('\t')[0]] = line.split('\n')[0] 
tree = Phylo.read('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/CPA_phylogeny_midroot_TREE_ALIGNEDRED_7Nov_ft_lg_cat_gamma.nw', 'newick')

count = 0
for clade in tree.find_clades():
    if clade.is_terminal() == True:
        pass
    else:
        clade.comment = clade.name
        clade.name = str(count)
        count = count + 1
#only give clade a number/name if it is not terminal
# otherwise let it keep the same name, this prevents future back paddeling and renaming. 
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/info_RED_clades_CPA_phylogeny_FTtree_lg_cat_gamma_7Nov_plus_extended_info.tsv', 'w') as extra_info:
    extra_info.write(header)
    for clade in tree.find_clades():
        if clade.is_terminal() == True:
            newline = oldline + '\t' + 'termini' + '\t' + 'termini' + '\n'
            extra_info.write(newline)
        else:
            if clade.name in Info_dict:
                if clade.name != '0':
                    oldline =  (Info_dict[clade.name])
                    termini_list = []
                    for termini in clade.get_terminals():
                        termini_list.append(termini.name)
                    daughter_node_list =[]
                    for i in clade.find_clades():
                        daughter_node_list.append(i)
                    newline = oldline + '\t' + str(termini_list) + '\t' + str(daughter_node_list) + '\n'
                    extra_info.write(newline)
