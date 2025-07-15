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
    
                elif scan_dict[full_id][1].split('-') > evalue.split('-'):
                    pass
                else:
                    entry_list = (best_match_hmm, evalue)
                    scan_dict[full_id] =entry_list
    return(scan_dict)
import os

import os
from collections import Counter

scan_dict = {}                    


count = 0

for hmmscan in os.listdir( '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/cross_domain'):
    if hmmscan.split('.')[-1] == 'tsv':
        HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/cross_domain/{}'.format(hmmscan)
        print(HMMscan)
        scan_dict = best_hit_dict(HMMscan, scan_dict)

tax_list = []    
Pfam_hitlist = []
CPA_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
scanned_dict = {}
for key in scan_dict.keys():
    value = scan_dict[key]
    if value[0] in CPA_list:
        scanned_dict[key] = value
print(len(list(scanned_dict.keys())))
tax_list = []
for key in scanned_dict.keys():
    tax_list.append(key.split('tax:')[1])
print(len(set(tax_list)))

from Bio import SeqIO

Euk_fldict = {}

for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain/Euk_cross_domain_fulllength_pf00999.faa', 'fasta'):
    Seq = str(record.seq)
    Euk_fldict[record.id] = Seq

Euk_aligneddict = {}
Seq = ''
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain/Euk_cross_domain_hmmaligned_pf00999.faa', 'fasta'):
    Seq = str(record.seq)
    Euk_aligneddict[record.id] = Seq


CPA_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
scanned_dict = {}
for key in scan_dict.keys():
    value = scan_dict[key]
    if value[0] in CPA_list:
        scanned_dict[key] = value


with open('Eukarya_all_hmmscanned_aligned.tsv', 'a') as E:
    header = 'protein_id' + '\t' + 'tax' + '\t' + 'full_length' + '\t' + 'aligned' + '\n'
    E.write(header)
    with open('Eukarya_all_hmmscanned_aligned_treeinput.fasta', 'a') as O:
        for key in scanned_dict.keys():
            aligned = (Euk_aligneddict[key])
            full_length = Euk_fldict[key]
            protein_id = key.split('_tax:')[0]
            tax = key.split('tax:')[1]
            Line_tsv = protein_id + '\t' + tax + '\t' + full_length + '\t' + aligned + '\n'
            fasta = '>' + '{}'.format(key) + '\n' + aligned + '\n'
            O.write(fasta)
            E.write(Line_tsv)


with open('Eukarya_all_hmmscanned_hmmaligned_fl.fasta', 'a') as O:
        for key in scanned_dict.keys():
  
            full_length = Euk_fldict[key]

            fasta = '>' + '{}'.format(key) + '\n' + full_length + '\n'
            O.write(fasta)
