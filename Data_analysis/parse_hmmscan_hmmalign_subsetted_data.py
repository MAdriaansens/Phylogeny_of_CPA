from Bio import SeqIO
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
CPA_hit_list = ['Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Na_H_antiport_1', 'CHX17_2nd', 'CHX_17', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned']
scanned_list = []

import os
Scan_dict1 = {}
#for Bacteria HMMscan_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Bacteria/PF00999/'

HMMscan_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea'
not_cpa_list = []

for entry in os.listdir(HMMscan_dir):
    if 'All_euk' in entry:
        scan_dict = {}
        HMMscan = '{}/{}'.format(HMMscan_dir, entry)
        scan_dict = best_hit_dict(HMMscan, scan_dict)
        Scan_dict1.update(scan_dict)
    
        
for key in Scan_dict1.keys():
    if Scan_dict1[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(Scan_dict1[key][0])

Scan_dict1 = {}

for entry in os.listdir(HMMscan_dir):
    if 'All_bac' in entry:
        scan_dict = {}
        HMMscan = '{}/{}'.format(HMMscan_dir, entry)
        scan_dict = best_hit_dict(HMMscan, scan_dict)
        Scan_dict1.update(scan_dict)
        
for key in Scan_dict1.keys():
    if Scan_dict1[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(Scan_dict1[key][0])

Scan_dict1 = {}

for entry in os.listdir(HMMscan_dir):
    if 'All_Arc' in entry:
        #for Bacteria it should be All_arc
        scan_dict = {}
        HMMscan = '{}/{}'.format(HMMscan_dir, entry)
        scan_dict = best_hit_dict(HMMscan, scan_dict)
        Scan_dict1.update(scan_dict)
        
for key in Scan_dict1.keys():
    if Scan_dict1[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(Scan_dict1[key][0])

print(len(set(scanned_list)))

#if you want to parse the not CPA
from collections import Counter
print(len(not_cpa_list))
print(len(scanned_list))
print(Counter(not_cpa_list))

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain'
#For Bacteria HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/cross_domain_final'
Seq_dict = {}
for j in os.listdir(HMMalign_dir):
    if 'faa.fasta' in j:
        for record in SeqIO.parse('{}/{}'.format(HMMalign_dir, j), 'fasta'):
            Seq_dict[record.id] = record.seq
print(len(Seq_dict.keys()))

with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/ARC_hmmscanned_hmmaligned.fasta', 'a') as O:
#For Bacteria with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/BAC_hmmscanned_hmmaligned.fasta', 'a') as O:

    for entry in set(scanned_list):
        sequence_id = entry.split('_tax')[0]
        sequence = Seq_dict[entry]
        Line = '>{}'.format(sequence_id) + '\n' + str(sequence) + '\n'
        O.write(Line)
