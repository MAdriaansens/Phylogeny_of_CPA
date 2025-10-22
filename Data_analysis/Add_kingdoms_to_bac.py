Kingdom_dict = {}
with open('/nesi/nobackup/uc04105/new_databases_May/Phyla_Bac226kingdoms - Phyla.tsv', 'r') as x:
    next(x, None)
    for line in x:
        Phyla = line.split('\t')[0]
        Kingdom = line.split('\t')[1]
        Kingdom_dict[Phyla] = Kingdom


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Bacteria_clusters_12sep.tsv', 'r') as clusters:
    for line in clusters:
        header = line
        break
print(header)


with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Bacteria_clusters_22sep.tsv', 'w') as Main:
    Title = header.split('\n')[0] + '\t' + 'Multi_kingdom' + '\t' + 'Which_kingdoms' + '\n'
    Main.write(Title)
    with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/RED_values_Bacteria_clusters_12sep.tsv', 'r') as clusters:
        next(clusters, None)
        for cluster in clusters:
            Single_Kingdoms = 'yes'
            if len(cluster.split('\t')[8].split(",")) > 1:
                Which_Kingdoms = []

                phylas = (cluster.split('\t')[8][2:-2].replace("'", "").replace("p__", "").replace(' ', '').split(","))
                for phyla in phylas:
                    Kingdom = Kingdom_dict[phyla]
                    Which_Kingdoms.append(Kingdom)
            else:
                Which_Kingdoms = []
                Phyla = str(cluster.split('\t')[8][5:-2])
                Kingdom = Kingdom_dict[Phyla]
                Which_Kingdoms.append(Kingdom)
            if len(set(Which_Kingdoms)) != 1:
                Multi_kingdom = 'yes'
                Which_Kingdoms = set(Which_Kingdoms)
            else:
                Multi_kingdom = 'no'
                Which_Kingdoms = set(Which_Kingdoms)
            Line = cluster.split('\n')[0] + '\t' + str(Multi_kingdom) + '\t' + str(Which_Kingdoms) + '\n'
            Main.write(Line)
