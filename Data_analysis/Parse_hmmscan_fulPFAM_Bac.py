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
scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/subset_sequence/Bacteria'
from Bio import SeqIO

scan_dict = {}


import os
for hmmscan in os.listdir(scandir):
    if hmmscan.split('.')[-1] == 'tsv':
        HMMscan = '{}/{}'.format(scandir,hmmscan)
        scan_dict = best_hit_dict(HMMscan, scan_dict)
        print(HMMscan)
tax_list = []
Pfam_hitlist = []
CPA_HMM_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
scanned_dict = {}
for key in scan_dict.keys():
    value = scan_dict[key]
    if value[0] in CPA_HMM_list:
        if key not in list(scanned_dict.keys()):
            scanned_dict[key] = value
        else:
            pass
print(len(list(scanned_dict.keys())))
CPA_list = []
for key in scanned_dict.keys():
    CPA_list.append(key)
print(len(set(CPA_list)))
print(CPA_list[-1])

#returns all protein ids whom match with CPA PFAM
Passed_all_list = {}
HMMalign = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/cross_domain_final/All_merged_Bacteria_pf00999aligned.faa'
for record in SeqIO.parse(HMMalign, 'fasta'):
    if record.id in CPA_list:
        Passed_all_list[record.id] = str(record.seq)
print(len(set(list(Passed_all_list.keys()))))
#NhaB

scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan'

HMMscan = '{}/Bacteria_06450_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'NhaB':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF06450/MERGED/PF06450_Bacteria_all_merged_aligned.fasta','fasta'):
    if record.id.split('|')[1] in hit_list:
        passed_list.append(record.id.split('tax:')[1])

NhaB_list = passed_list
print(len(set(NhaB_list)))

#NhaC
HMMscan = '{}/Bacteria_03553_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'Na_H_antiporter':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03553/MERGED/PF03553_Bacteria_all_merged_aligned.fasta', 'fasta'):
    if record.id.split('|')[1] in hit_list:
        passed_list.append(record.id.split('tax:')[1])

NhaC_list = passed_list
print(len(set(NhaC_list)))


#NhaD
HMMscan = '{}/Bacteria_03600_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'CitMHS':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03600/PF03600_merged_Bac_fl.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id.split('tax:')[1])
NhaD_list = passed_list
print(len(set(NhaD_list)))
scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan'
count = 0
Tax_dict={}
representatives_list = []
#add sample site
CPA_list_forIT = []
for key in Passed_all_list.keys():
    CPA_list_forIT.append(key.split('tax:')[1])

with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/bac120_metadata.tsv', 'r') as Meta:
    with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_Bacteria_3OKT.tsv', 'w') as Out:
        header = 'GTDB_id' + '\t' + 'GTDB_tax' + '\t' + 'Completeness' + '\t' + 'Contamination' + '\t' + 'sample' + '\t'+ 'CPA_count' + '\t' + 'CPA_binary' + '\t' + 'NhaB_count' + '\t' +  'NhaB_binary' + '\t' + 'NhaC_count' + '\t' + 'NhaC_binary'+'\t' + 'NhaD_count' + '\t' + 'NhaD_binary' + '\n'
        Out.write(header)
        for line in Meta:
            if line.split('\t')[18] == 't':
                GTDB_id =line.split('\t')[0]
                completeness = line.split('\t')[2]
                contamination = line.split('\t')[3]
                sample = line.split('\t')[-52]
                GTDB_tax = (line.split('\t')[19].replace(' ', '_'))
                representatives_list.append(GTDB_tax)

                NhaB_count = NhaB_binary = NhaC_count = NhaC_binary = NhaD_count = NhaD_binary = CPA_count = CPA_binary= 0
                if GTDB_tax in NhaD_list:
                    NhaD_count = NhaD_list.count(GTDB_tax)
                    NhaD_binary = 1
                if GTDB_tax in NhaC_list:
                    NhaC_count = NhaC_list.count(GTDB_tax)
                    NhaC_binary = 1
                if GTDB_tax in NhaB_list:
                    NhaB_count = NhaB_list.count(GTDB_tax)
                    NhaB_binary = 1
                if GTDB_tax in CPA_list_forIT:
                    CPA_count = CPA_list_forIT.count(GTDB_tax)
                    CPA_binary = 1
                Tax_dict[GTDB_id] = GTDB_tax
                Wline = GTDB_id + '\t' + GTDB_tax + '\t' + str(completeness) + '\t' + str(contamination) + '\t' + sample + '\t' + str(CPA_count) + '\t' + str(CPA_binary) + '\t' + str(NhaB_count) + '\t' +  str(NhaB_binary) + '\t' + str(NhaC_count) +'\t' + str(NhaC_binary) +'\t' + str(NhaD_count) + '\t' + str(NhaD_binary) + '\n'
                Out.write(Wline)
            else:
                pass
print(count)
print(len(representatives_list))
print(len(set(CPA_list_forIT)))
print(len(CPA_list_forIT))
tax = []
import json
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_tax_3OKT.json', 'w') as outfile:
    json.dump(Tax_dict, outfile)

with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_CPA_idseq_3OKT.json', 'w') as outfile:
    json.dump(Passed_all_list, outfile)

with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_passed_all_filters_OKT3_GTDBreps_alignedPF00999.fasta', 'w') as Passed:
    for key in Passed_all_list.keys():
        if key.split('tax:')[1] in representatives_list:
            header = key.split('_tax')[0]
            tax.append(key.split('tax:')[1])
            sequence = Passed_all_list[key]
            line = '>' + header + '\n' + str(sequence) + '\n'
            Passed.write(line)
        else:
            print(key)
print(len(set(tax)))
