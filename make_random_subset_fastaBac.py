from Bio import SeqIO
import json
import os
import random

direction='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Bacteria/PF00999'
file_list = os.listdir(direction)

ID_list = []

for i in file_list:
    if 'json' not in i:
        pass
    else:

        f = open('{}/{}'.format(direction, i))
        data = json.load(f)
        for j in data.keys():
            main_hmm = data[j][0]
            if 'PF00999' != main_hmm:
                pass
            elif 'PF00999' == main_hmm:
                ID_list.append(j)

CPA_id_list = list(set(ID_list))
                
            
infile= '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/MERGED/Bac_merged_all_hmmscanned.faa'




hit_dict = {}
for record in SeqIO.parse('{}'.format(infile), 'fasta'):
    hit_dict[record.id]=record.seq

def return_fasta_from_list(hit_dict, random_list2, count, outdir):
    written_list = []
    with open('{}/Bacteria_HMMscanned_PF00999_aligned_randomset_{}.fasta'.format(outdir, count), 'a') as OUT:
        for i in random_list2:
            sequence = hit_dict[i]print(len(list(hit_dict.keys())))
            outline = '>' + i + '\n' + str(sequence)  + '\n'
            OUT.write(outline)

outdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/MERGED'
for count in range(1,101):
    random_list = random.choices(CPA_id_list, k=100)
    if len(list(set(random_list))) != 100:
        for i in random.sample(CPA_id_list, k=len(CPA_id_list)):
            if len(list(set(random_list))) < 100:
                random_list.append(i)
            elif len(list(set(random_list))) == 100:
                break
    random_list2 = list(set(random_list))
    print(len(random_list2))
    random_list = []
    return_fasta_from_list(hit_dict, random_list2, count, outdir)
