def best_hit_dict(HMMscan, scan_dict):
    count = 0
    with open(HMMscan, 'r') as H:
        for line in H:
            if line[0] != '#':
                e_count =0
                best_match_hmm = line.split(' ')[0].split(' ')[-1]
                protein_id = line.split('_tax:')[0].split(' ')[-1]
                Taxonomy = (line.split('_tax:')[-1].split(' - ')[0])
                full_id = protein_id + '_tax:' + Taxonomy
    
                if line.count(' - ') == 2:
                    prelim_evalue = line.split(' - ')[2]
                elif line.count(' - ') == 1:
                    prelim_evalue = line.split(' - ')[1]
                else:
                    print(line)
                    break
                for i in (list(set(prelim_evalue.split(' ')))):
                    if e_count == 0:
                        if 'e-' in i:
                            evalue = i
                            e_count =  1
                        else:
                            pass
                if full_id not in scan_dict.keys():
                    entry_list = (best_match_hmm, evalue)
                    scan_dict[full_id] = entry_list
    
                elif scan_dict[full_id][1].split('-')[1] > evalue.split('-')[1]:

                    pass
                else:
                    entry_list = (best_match_hmm, evalue)
                    scan_dict[full_id] =entry_list
    return(scan_dict)

scanned_dict = {}

print(len(scanned_dict.keys()))

for key in scanned_dict.keys():
    print(key)
    break
from Bio import SeqIO

import os

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain'
Seq_dict = {}
for j in os.listdir(HMMalign_dir):
    if 'fa.fasta' in j:
        for record in SeqIO.parse('{}/{}'.format(HMMalign_dir, j), 'fasta'):
            Seq_dict[record.id] = record.seq

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/Euk_hmmscanned_hmmaligned.fasta', 'a') as O:

    for key in scanned_dict.keys():
        sequence_id = key.split('_tax')[0]
        sequence = Seq_dict[key]
        Line = '>{}'.format(sequence_id) + '\n' + str(sequence) + '\n'
        O.write(Line)
