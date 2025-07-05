import sys
import json

HMMscan = sys.argv[1] #'/nesi/nobackup/uc04105/results/HMM/HMM_db/PF03600Archaea_HMMscanned.tsv'
HMM = sys.argv[2]
out_json = sys.argv[3]

scan_dict ={}
count = 0
with open(HMMscan, 'r') as H:
    for line in H:
        if line[0] != '#':
            e_count =0
            best_match_hmm = line.split('.')[0].split(' ')[-1]
            protein_id = line.split('_tax:')[0].split(' ')[-1]
            Taxonomy = (line.split('_tax:')[-1].split(' - ')[0])
            full_id = protein_id + '_tax:' + Taxonomy

            #chunk added, if HMMs are made inhouse vs Pfam the '-' split becomes a problem, thus an if else statement is implemented
            if line.count(' - ') ==1:
                prelim_evalue = line.split(' - ')[1]
            else: 
                prelim_evalue = line.split(' - ')[2]
            
            
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
                
            elif scan_dict[full_id][1].split('-') > evalue.split('-'):
                pass
            else:
                 scan_dict[full_id] =entry_list

hit_list = []
for j in list(scan_dict.keys()):
    if '{}'.format(HMM) == scan_dict[j][0]:
        hit_list.append(j)
with open('{}'.format(out_json), 'w') as f:
    json.dump(hit_list, f)
