directory = '/nesi/nobackup/uc04105/fasta_files/pipeline/CPA/tree_16082024/input_files_tree/cat/cluster_mmseqs'

diction = '/nesi/nobackup/uc04105/Bacteria_hits_full_mar17.json'
tsv_cluster = '{}/BacteriaclusterResSI7C8_cluster.tsv'.format(directory)
reps = '{}/BacteriaclusterResSI7C8_rep_seq.fasta'.format(directory)

#diction contains all the protein_ids and matching accession id + taxonomy


#get signles and multiples:

blank_list = []
with open(tsv_cluster, 'r') as clustermap:
    for entry in clustermap:
        representative = entry.split('\t')[0]
        blank_list.append(representative)

clustermap.close()

from collections import Counter
import json

singletonlist = []
multiplelist = []
blank = Counter(blank_list)
items = blank.items()
for i in items:
    if i[-1] == 1:
        singletonlist.append(i[0])
    else:
        multiplelist.append(i[0])

print(len(singletonlist))
print(len(multiplelist))

replist = []
Rep_dict = {}
with open(tsv_cluster, 'r') as clustermap:
    count = 0
    second_list = []
    element_count = 0
    for entry in clustermap:
        count = blank[entry.split('\t')[0]]
        repm = entry.split('\t')[0]
        if repm in multiplelist:
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
        #latest element to check in, drop count = 1 line as in multiple list does the same
import json
def get_taxonomy_cluster(cluster_rep):
    with open(diction) as json_file:
        data = json.load(json_file)


    entry = Rep_dict[cluster_rep]
    num_proteins = (len(entry))
    species_list = []
    genus_list = []
    family_list = []
    order_list = []
    class_list = []
    phylum_list = []
    for prot_id in entry:
        species_list.append(data[prot_id][-1].split(';s__')[1])
        genus_list.append(data[prot_id][-1].split(';g__')[1].split(';s')[0])
        family_list.append(data[prot_id][-1].split(';f__')[1].split(';g')[0])
        order_list.append(data[prot_id][-1].split(';o__')[1].split(';f')[0])
        class_list.append(data[prot_id][-1].split(';c__')[1].split(';o')[0])
        phylum_list.append(data[prot_id][-1].split(';p__')[1].split(';c')[0])
    species_list = list(set(species_list))
    genus_list = list(set(genus_list))
    family_list = list(set(family_list))
    order_list = list(set(order_list))
    class_list = list(set(class_list))
    phylum_list = list(set(phylum_list))

    ls = len(species_list)
    lg = len(genus_list)
    lf = len(family_list)
    lo = len(order_list)
    lc = len(class_list)
    lp = len(phylum_list)
    string = 'num proteins: {}, numb uniq s: {}, num uniq g: {},  number of uniq f: {},  num uniq o: {}, num uniq c: {}, num uniq p {}'.format(num_proteins, ls, lg, lf, lo, lc, lp)
    if lf > 1:
        return(num_proteins)
    else:
        return(0)
