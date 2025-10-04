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
from Bio import SeqIO

#NhaB
HMMscan = '{}/Archaea_06450_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'NhaB':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/Archaea_merged_PF06450_unique_aligned_PF06450.fasta.faa.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id.split('tax:')[1])

NhaB_list = passed_list
print(len(NhaB_list))



#NhaC
HMMscan = '{}/Archaea_03553_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'Na_H_antiporter':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/Archaea_merged_PF03553_unique_aligned_PF03553.fasta.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id.split('tax:')[1])

NhaC_list = passed_list
print(len(NhaC_list))

#NhaD
HMMscan = '{}/Archaea_03600_fullhmmscanned.tsv'.format(scandir)
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'CitMHS':
        hit_list.append(key)

passed_list = []
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/Archaea_merged_PF03600_unique_aligned_PF03600.fasta.fasta', 'fasta'):
    if record.id in hit_list:
        passed_list.append(record.id.split('tax:')[1])
NhaD_list = passed_list
print(len(NhaD_list))



scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/subset_sequence/Archaea'
import os
for hmmscan in os.listdir(scandir):
    if hmmscan.split('.')[-1] == 'tsv':
        HMMscan = '{}/{}'.format(scandir,hmmscan)
        scan_dict = best_hit_dict(HMMscan, scan_dict)
        print(HMMscan)

tax_list = []    
Pfam_hitlist = []

scanned_dict = {}
CPA_HMM_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Na_H_antiporter','Nha1_C','Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
for key in scan_dict.keys():
    value = scan_dict[key]
    if value[0] in CPA_HMM_list:
        scanned_dict[key] = value
print(len(list(scanned_dict.keys())))
CPA_list = []
for key in scanned_dict.keys():
    CPA_list.append(key.split('_subset')[0])

print(len(CPA_list))
#returns all protein ids whom match with CPA PFAM


print(len(CPA_list))
#returns all protein ids whom match with CPA PFAM
Passed_all_list={}
HMMalign = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain/Arc_cross_domain_hmmaligned_pf00999.faa'
for record in SeqIO.parse(HMMalign, 'fasta'):
    if record.id in scanned_dict:
        Passed_all_list[record.id] = record.seq
GTDB_ids_passed = {}

CPA_list_forIT = []
for key in Passed_all_list.keys():
    CPA_list_forIT.append(key.split('tax:')[1])
print(len(CPA_list_forIT))
print(len(set(CPA_list_forIT)))


                          
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/ar53_metadata.tsv', 'r') as Meta:
    with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_Archaea_2OKTSEPT.tsv', 'w') as Out:
        header = 'GTDB_id' + '\t' + 'GTDB_tax' + '\t' + 'Completeness' + '\t' + 'Contamination' + '\t' + 'Sample' + '\t'+ 'CPA_count' + '\t' + 'CPA_binary' + '\t' + 'NhaB_count' + '\t' +  'NhaB_binary' + '\t' + 'NhaC_count' + '\t' + 'NhaC_binary'+'\t' + 'NhaD_count' + '\t' + 'NhaD_binary' + '\n'
        Out.write(header)
        representatives_list = []
        for line in Meta:
            if line.split('\t')[18] != 't':
                pass
            else:
                GTDB_id =line.split('\t')[0]
                Sample = line.split('\t')[-52]
                completeness = line.split('\t')[2]
                contamination = line.split('\t')[3]
                GTDB_tax = (line.split('\t')[19].replace(' ', '_'))

                NhaB_count = NhaB_binary = NhaC_count = NhaC_binary = NhaD_count = NhaD_binary = CPA_count = CPA_binary= 0
                GTDB_ids_passed[GTDB_tax] = GTDB_id
                representatives_list.append(GTDB_tax)

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
                Wline = GTDB_id + '\t' + GTDB_tax + '\t' + str(completeness) + '\t' + str(contamination) + '\t' + Sample + '\t' + str(CPA_count) + '\t' + str(CPA_binary) + '\t' + str(NhaB_count) + '\t' +  str(NhaB_binary) + '\t' + str(NhaC_count) +'\t' + str(NhaC_binary) +'\t' + str(NhaD_count) + '\t' + str(NhaD_binary) + '\n'
                Out.write(Wline)
tax = []




with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Archaea_passed_all_filters_2OKT_GTDBreps_alignedPF00999.fasta', 'w') as Passed:
    for key in Passed_all_list.keys():
        if key.split('tax:')[1] in representatives_list:
            header = key.split('_tax')[0]
            tax.append(key.split('tax:')[1])
            sequence = Passed_all_list[key]
            line = '>' + header + '\n' + str(sequence) + '\n'
            Passed.write(line)

print(len(set(tax)))
