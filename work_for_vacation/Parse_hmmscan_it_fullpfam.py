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
                else:
                    prelim_evalue = line.split(' - ')[1]
                
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

scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan'
HMMscan = '{}/Bacteria_03600_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'CitMHS':
        hit_list.append(key)
print(len(hit_list))

from Bio import SeqIO
passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/Archaea_merged_PF06450_unique_aligned_PF06450.fasta.faa.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id)
print(len(passed_list))
