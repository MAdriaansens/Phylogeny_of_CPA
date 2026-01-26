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

Pfam_count = 0
Eukcount =0
Arc_count =0
Bac_count = 0
count = 0
for key in scan_dict.keys():
    count =  count + 1
    if 'Na_H_Exchanger' == (scan_dict[key][0]):
        Pfam_count = Pfam_count + 1
    if 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned' == (scan_dict[key][0]):
        Eukcount = Eukcount + 1
    if 'Archaea_Manual_seq_e05_cov30_seqid0.6.faa_rep_seq_globalpairaligned' == (scan_dict[key][0]):
        Arc_count = Arc_count + 1
    elif 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned' == (scan_dict[key][0]):
        Bac_count = Bac_count + 1
print((Pfam_count/count)*100, (Bac_count/count)*100, (Eukcount/count)*100, (Arc_count/count)*100)
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/BAC_HMM_e03vsEukarya_fl_INH_scanned.tsv'#BAC HMM



scan_dict = {}                    
scan_dict = best_hit_dict(HMMscan)

Pfam_count = 0
Eukcount =0
Arc_count =0
Bac_count = 0
count = 0
for key in scan_dict.keys():
    count =  count + 1
    if 'Na_H_Exchanger' == (scan_dict[key][0]):
        Pfam_count = Pfam_count + 1
    if 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned' == (scan_dict[key][0]):
        Eukcount = Eukcount + 1
    if 'Archaea_Manual_seq_e05_cov30_seqid0.6.faa_rep_seq_globalpairaligned' == (scan_dict[key][0]):
        Arc_count = Arc_count + 1
    elif 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned' == (scan_dict[key][0]):
        Bac_count = Bac_count + 1
print((Pfam_count/count)*100, (Bac_count/count)*100, (Eukcount/count)*100, (Arc_count/count)*100)

#HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/Eukarya_allmerged_PF00999_INH_HMMSCANNED.TSV' #EUK TO be done



scan_dict = {}                    
scan_dict = best_hit_dict(HMMscan)

Pfam_count = 0
Eukcount =0
Arc_count =0
Bac_count = 0
count = 0
for key in scan_dict.keys():
    count =  count + 1
    if 'Na_H_Exchanger' == (scan_dict[key][0]):
        Pfam_count = Pfam_count + 1
    if 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned' == (scan_dict[key][0]):
        Eukcount = Eukcount + 1
    if 'Archaea_Manual_seq_e05_cov30_seqid0.6.faa_rep_seq_globalpairaligned' == (scan_dict[key][0]):
        Arc_count = Arc_count + 1
    elif 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned' == (scan_dict[key][0]):
        Bac_count = Bac_count + 1
print((Pfam_count/count)*100, (Bac_count/count)*100, (Eukcount/count)*100, (Arc_count/count)*100)

#HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/Eukarya_allmerged_PF00999_INH_HMMSCANNED.TSV' #EUK TO be done
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/ARC_HMM_e03vsEukarya_fl_INH_scanned.tsv' #ARC HMM



scan_dict = {}                    
scan_dict = best_hit_dict(HMMscan)

Pfam_count = 0
Eukcount =0
Arc_count =0
Bac_count = 0
count = 0
for key in scan_dict.keys():
    count =  count + 1
    if 'Na_H_Exchanger' == (scan_dict[key][0]):
        Pfam_count = Pfam_count + 1
    if 'Manual_seq_cov30_e05_seqid0.7_genafpair_aligned' == (scan_dict[key][0]):
        Eukcount = Eukcount + 1
    if 'Archaea_Manual_seq_e05_cov30_seqid0.6.faa_rep_seq_globalpairaligned' == (scan_dict[key][0]):
        Arc_count = Arc_count + 1
    elif 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned' == (scan_dict[key][0]):
        Bac_count = Bac_count + 1
print((Pfam_count/count)*100, (Bac_count/count)*100, (Eukcount/count)*100, (Arc_count/count)*100)
