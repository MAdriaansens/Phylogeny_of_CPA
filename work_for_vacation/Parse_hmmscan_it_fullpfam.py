#make sure this gets incorperated into the full HMM of all GTDB

def best_hit_dict(HMMscan, scan_dict):
    #open a HMMscan file provided
    with open(HMMscan, 'r') as H:
        for line in H:
            if line[0] != '#':
                #skip lines with # as this contains information not needed
    
                e_count =0
                #marker for a 

                best_match_hmm = line.split(' ')[0].split(' ')[-1]
                protein_id = line.split('_tax:')[0].split(' ')[-1]
                Taxonomy = (line.split('_tax:')[-1].split(' - ')[0])
                full_id = protein_id + '_tax:' + Taxonomy
    
                if line.count(' - ') == 2:
                    prelim_evalue = line.split(' - ')[2]

                else:
                    prelim_evalue = line.split(' - ')[1]
                    if 'Fe' in prelim_evalue:
                        prelim_evalue = line.split(' - ')[1].split('Fe')[0]
                    else:
                        prelim_evalue = line.split(' - ')[1]
                e_list = []
                for i in prelim_evalue.split(' '):
                    if i != '':
                        e_list.append(i)
                j = e_list[0]
                if 'e' in j:
                    evalue = pow(10,int(j.split('e')[1]))*float(j.split('e')[0])
                
                else:
                    if any(x.isalpha() for x in j) == True:
                        print(j)
                        break
                    else:
                        evalue = float(j)
                    
                if full_id not in scan_dict.keys():
                    entry_list = (best_match_hmm, evalue)
                    scan_dict[full_id] = entry_list
    
                elif scan_dict[full_id][1] < evalue:
                    pass
                else:
                    entry_list = (best_match_hmm, evalue)
                    scan_dict[full_id] =entry_list
    return(scan_dict)

scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan'
HMMscan = '{}/Archaea_03553_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'Na_H_antiporter':
        hit_list.append(key)
print(len(hit_list))

from Bio import SeqIO
passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/Archaea_merged_PF03553_unique_aligned_PF03553.fasta.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id.split('tax:')[1])
print(len(set(passed_list)))

scan_dict = best_hit_dict(HMMscan, scan_dict)
