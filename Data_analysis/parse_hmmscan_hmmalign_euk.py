from Bio import SeqIO
not_cpa_list = []

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
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/cross_domain/Self_6juli_merged_MMSEQ_vsEukarya_PF00999hmmaligned_fl_FULLPFAMHMMscanned.tsv' #PFAM
scan_dict = {}


scan_dict = best_hit_dict(HMMscan, scan_dict)

for key in scan_dict.keys():
    if scan_dict[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(scan_dict[key][0])

HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/cross_domain/7JULI_ARCHAEA_vsEukarya_PF00999hmmaligned_fl_FULLPFAMHMMscanned.tsv' #PFAM
scan_dict = {}


scan_dict = best_hit_dict(HMMscan, scan_dict)

for key in scan_dict.keys():
    if scan_dict[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(scan_dict[key][0])
        
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/cross_domain/7JULI_BACTERIA_vsEukarya_PF00999hmmaligned_fl_FULLPFAMHMMscanned.tsv' #PFAM
scan_dict = {}


scan_dict = best_hit_dict(HMMscan, scan_dict)

for key in scan_dict.keys():
    if scan_dict[key][0] in CPA_list:
        scanned_list.append(key)
    else:
        not_cpa_list.append(scan_dict[key][0])

from collections import Counter
print(len(not_cpa_list))
print(len(scanned_list))
print(Counter(not_cpa_list))

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain'
Seq_dict = {}
for j in os.listdir(HMMalign_dir):
    if 'fa.fasta' in j:
        for record in SeqIO.parse('{}/{}'.format(HMMalign_dir, j), 'fasta'):
            Seq_dict[record.id] = record.seq
with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/sequences/Euk_hmmscanned_hmmaligned.fasta', 'a') as O:

    for entry in set(scanned_list):
        sequence_id = entry.split('_tax')[0]
        sequence = Seq_dict[entry]
        Line = '>{}'.format(sequence_id) + '\n' + str(sequence) + '\n'
        O.write(Line)

