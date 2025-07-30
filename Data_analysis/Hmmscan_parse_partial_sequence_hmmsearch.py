#parse subset HMMscan

from Bio import SeqIO
not_cpa_list = []

def best_hit_dict(HMMscan, scan_dict):
    count = 0
    with open(HMMscan, 'r') as H:
        for line in H:
            if line[0] != '#':
                e_count =0

                
                protein_id = line.split('_subset')[0].split(' ')[-1]
                best_match_hmm = line.split(' ')[0].split(' ')[0]
                full_id = protein_id
                if line.count(' - ') == 2:
                    prelim_evalue = line.split(' - ')[2]
                elif line.count(' - ') == 1:
                    prelim_evalue = line.split(' - ')[1]
                else:
                    print(line)
                    
                #print(prelim_evalue)

                for i in (list(set(prelim_evalue.split(' ')))):
                    if e_count == 0:
                        if 'e-' in i:
                            if i[0].isalpha() == False:
                                evalue = i
                                e_count =  1
                            else:
                                pass
                        else:
                            pass
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
scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/subset_sequence'
import os
from collections import Counter
CPA_list = ['Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned']
scanned_list = []
not_cpa_list = []

for hmmscan in os.listdir(scandir):
    scanned_list = [] 
    not_cpa_list = []
    HMMscan = '{}/{}'.format(scandir, hmmscan)
    
    scan_dict = {}
    scan_dict = best_hit_dict(HMMscan, scan_dict)
    for key in scan_dict.keys():

        if scan_dict[key][0] in CPA_list:
            scanned_list.append(scan_dict[key][0])
        else:
            not_cpa_list.append(scan_dict[key][0])
    print(HMMscan)
    print(len(not_cpa_list))
    print(len(scanned_list))
    print(Counter(not_cpa_list))
