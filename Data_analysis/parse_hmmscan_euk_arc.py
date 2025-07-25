#this works for both Eukarya and Archaea since they are one singular hmmscan .tsv file per domain, Bacteria has been split up and the code is a bit different
def best_hit_dict(HMMscan):
    scan_dict ={}
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
    
                elif scan_dict[full_id][1].split('-') > evalue.split('-'):
                    pass
                else:
                     scan_dict[full_id] =entry_list
    return(scan_dict)
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/Eukarya_allmerged_PF00999_INH_HMMSCANNED.TSV' #PFAM



scan_dict = {}                    
scan_dict = best_hit_dict(HMMscan)

scan_dict = best_hit_dict(HMMscan)
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea/BAC_HMM_e03vsArchaea_hmmsearchFULLPFAMscanned.tsv' #PFAM
scan_dict = {}


scan_dict = best_hit_dict(HMMscan, scan_dict)

Pfam_count = 0
Eukcount =0
Arc_count =0
Bac_count = 0
count = 0
Pfam_hitlist = []
for key in scan_dict.keys():
    count =  count + 1
    Pfam_hitlist.append(scan_dict[key][0])
print(count)
print(Counter(Pfam_hitlist))
