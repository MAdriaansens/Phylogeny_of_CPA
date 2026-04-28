from Bio import SeqIO
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
                    if j[0].isalpha() == True:
                        pass
                    else:
                        evalue = pow(10,int(j.split('e')[1]))*float(j.split('e')[0])
               



                else:
                    if any(x.isalpha() for x in j) == True:
                        pass
                        
                    elif 'diol' in j:
                        pass
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
scandir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Bacteria/PF00999/part3'
from Bio import SeqIO

scan_dict = {}


import os
for hmmscan in os.listdir(scandir):
    if hmmscan.split('.')[-1] == 'tsv':
        HMMscan = '{}/{}'.format(scandir,hmmscan)
        scan_dict = best_hit_dict(HMMscan, scan_dict)

tax_list = []
Pfam_hitlist = []
CPA_HMM_list = ('Manual_seq_cov30_e05_seqid0.7_genafpair_aligned', 'Na_H_antiporter','Nha1_C','Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned', 'Na_H_Exchanger', 'Manual_vsBacteria_merged_e5_cov30_seqid0.6.faa_rep_seq_autoaligned','Na_H_antiport_1', 'CHX17_2nd', 'CHX17_C')
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

CPA_tax_list = []

for CPA in CPA_list:
    CPA_tax_list.append(CPA.split('_tax:')[1])
print(len(set(CPA_tax_list)))
print(len(CPA_tax_list))

#returns all protein ids whom match with CPA PFAM > 70%
records_in ={}
for aln in os.listdir('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/part3'):
    if 'FL' in aln:
        for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/part3/{}'.format(aln), 'fasta'):
            records_in[record.id] = record.seq

Passed_all_list = []
with open('/nesi/nobackup/uc04105/new_databases_May/final_april28/Bacteria_CPA_fl.fasta', 'w') as out:
    for CPA in CPA_list:
        Passed_all_list.append(CPA)
        out.write('>' + CPA + '\n' + str(records_in[CPA]) + '\n')
##NhaB
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/results_HMMscan/PF06450_aligned_Bacteria_PF06450_retrieved_preQC_full_length_hmmaligned.fasta_aligned_hmmscanned.tsv'
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'NhaB':
        hit_list.append(key)

passed_list = []
fulllength_NhaC_file='/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/original_files/PF06450_merged_Bac_fl.fasta'
seq_dict = {}
for record in SeqIO.parse(fulllength_NhaC_file, 'fasta'):
    seq_dict[record.id] = record.seq
    
with open('/nesi/nobackup/uc04105/new_databases_May/final_april28/Bacteria_NhaB_fl.fasta', 'w') as B_out:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/original_files/PF03600_merged_Bac_fl.fasta', 'fasta'):
        if record.id in hit_list:
            passed_list.append(str(record.id.split('tax:')[1]))
            outline = '>' + record.id + '\n' + str(record.seq) + '\n'
            B_out.write(outline)
NhaB_list = passed_list
print(len(NhaB_list))
print(len(set(NhaB_list)))

#NhaC
#NhaC
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/results_HMMscan/PF03553_aligned_Bacteria_PF03553_retrieved_preQC_full_length_hmmaligned.fasta_aligned_hmmscanned.tsv'
scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'Na_H_antiporter':
        hit_list.append(key)

passed_list = []
fulllength_NhaC_file='/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/original_files/PF03553_merged_Bac_fl.fasta'
seq_dict = {}
for record in SeqIO.parse(fulllength_NhaC_file, 'fasta'):
    seq_dict[record.id] = record.seq
    
with open('/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/Bacteria_NhaC_fl.fasta', 'w') as C_out:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/original_files/PF03553_merged_Bac_fl.fasta', 'fasta'):
        if record.id in hit_list:
            passed_list.append(str(record.id.split('tax:')[1]))
            outline = '>' + record.id + '\n' + str(record.seq) + '\n'
            C_out.write(outline)
NhaC_list = passed_list
print(len(NhaC_list))
print(len(set(NhaC_list)))




    
#NhaD
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/results_HMMscan/PF03600_aligned_Bacteria_PF03600_retrieved_preQC_full_length_hmmaligned.fasta_aligned_hmmscanned.tsv'

scan_dict = best_hit_dict(HMMscan, scan_dict)
from collections import Counter
hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'CitMHS':
        hit_list.append(key)
print(len(hit_list))


passed_list = []
with open('/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/Bacteria_NhaD_fl.fasta', 'w') as D_out:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Chapter_2_data/original_files/PF03600_merged_Bac_fl.fasta', 'fasta'):
        if record.id in hit_list:
            
            passed_list.append(str(record.id.split('tax:')[1]))
            outline = '>' + record.id + '\n' + str(record.seq) + '\n'
            D_out.write(outline)
NhaD_list = passed_list
print(len(NhaD_list))

GTDB_ids_passed = {}
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/bac120_metadata_r226.tsv', 'r') as Meta:
    with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Antiporter_Bacteria_25April.tsv', 'w') as Out:
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
                if GTDB_tax in CPA_tax_list:
                    CPA_count =  CPA_tax_list.count(GTDB_tax)
                    CPA_binary = 1
                Wline = GTDB_id + '\t' + GTDB_tax + '\t' + str(completeness) + '\t' + str(contamination) + '\t' + Sample + '\t' + str(CPA_count) + '\t' + str(CPA_binary) + '\t' + str(NhaB_count) + '\t' +  str(NhaB_binary) + '\t' + str(NhaC_count) +'\t' + str(NhaC_binary) +'\t' + str(NhaD_count) + '\t' + str(NhaD_binary) + '\n'
                Out.write(Wline)
