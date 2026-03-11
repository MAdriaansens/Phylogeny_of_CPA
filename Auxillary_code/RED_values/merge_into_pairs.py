onzin = '/nesi/nobackup/uc04105/redvals/All_distances_all_pairs_Archaea.tsv'
output = '/nesi/nobackup/uc04105/redvals/All_distances_all_pairs_Archaea_edit.tsv'

with open(output, 'a') as R:
    header = 'pair' + '\t' + 'RED' + '\t' + 'Phylo' + '\n'
    R.write(header)
    with open(onzin, 'r') as O:
        next(O, None)
        for line in O:

            pair = line.split('\t')[0] + '*' + line.split('\t')[1] 
            RED = line.split('\t')[2]
            Phylo = line.split('\t')[3].split('\n')[0]
            line = pair + '\t' + RED + '\t' + Phylo + '\n'
            R.write(line)
