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
scandir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/subset'
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

#returns all protein ids whom match with CPA PFAM
Passed_all_list = {}
HMMalign = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/cross_domain/Eukarya_merged_PF00999hmmaligned.fasta'

for record in SeqIO.parse(HMMalign, 'fasta'):
    if record.id in scanned_dict:
        Passed_all_list[record.id] = record.seq
GTDB_ids_passed = {}
            
print(len(set(list(Passed_all_list.keys()))))



#NhaB
#HMMscan was performed on the part of the sequence whom matched with the HMMalign, then it was ran agianst the full pfam 
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/PF06450_Eukarya_merged_alignedPF06450_refilter_HMMscanned.tsv'


scan_dict = {}

scan_dict = best_hit_dict(HMMscan, scan_dict)
#note that HMMalign sometimes wants to add an '|' to the front of a protein id
scan_cleaned = {}
for key in scan_dict.keys():
    scan_cleaned[key] = scan_dict[key]
scan_dict = scan_cleaned

hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'NhaB':
        hit_list.append(key)


passed_list = []
records_in ={}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF06450_Eukarya_merged_alignedPF06450.fasta.fasta','fasta'):
    if record.id in list(set(hit_list)):
         passed_list.append(record.id.split('tax:')[1])

NhaB_list = passed_list
print(len(set(NhaB_list)))
print('finished NhaB')
#NhaC
scan_dict = {}
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/PF03553_Eukarya_merged_alignedPF03553_refilter_HMMscanned.tsv'

scan_dict = best_hit_dict(HMMscan, scan_dict)

scan_cleaned = {}
for key in scan_dict.keys():
    scan_cleaned[key.split('|')[1]] = scan_dict[key]
scan_dict = scan_cleaned



hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'Na_H_antiporter':
        hit_list.append(key)
print(len(hit_list))
print(len(set(hit_list)))


passed_list = []
records_in ={}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF03553_Eukarya_merged_alignedPF03553.fasta.fasta','fasta'):
    if record.id.split('|')[1] in list(set(hit_list)):
         passed_list.append(record.id.split('tax:')[1])

NhaC_list = passed_list
print(len(set(NhaC_list)))
print('finished NhaC')

#NhaD
scan_dict = {}
HMMscan = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMscan/PF03600_Eukarya_merged_alignedPF03600_refilter_HMMscanned.tsv'


scan_dict = best_hit_dict(HMMscan, scan_dict)

scan_cleaned = {}
for key in scan_dict.keys():
    scan_cleaned[key.split('|')[1]] = scan_dict[key]
scan_dict = scan_cleaned

hit_list = []
for key in scan_dict.keys():
    domain = scan_dict[key][0]
    if domain == 'CitMHS':
        hit_list.append(key)
print(len(hit_list))
print(hit_list[1])
print(len(set(hit_list)))

passed_list = []
records_in ={}
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF03600_Eukarya_merged_alignedPF03600.fasta.fasta','fasta'):

    if record.id.split('|')[1] in list(set(hit_list)):
         passed_list.append(record.id.split('tax:')[1])

NhaD_list = passed_list
print(len(set(NhaD_list)))


print('finished NhaD')


count = 0
Tax_dict={}
representatives_list = []
#add sample site
CPA_list_forIT = []
for key in Passed_all_list.keys():
    CPA_list_forIT.append(key.split('tax:')[1])
    
with open('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Eukarya_metadata.tsv', 'r') as Meta:
    next(Meta, None)
    with open('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/IT_Eukarya_4nov.tsv', 'w') as Out:
        header=  'species_name'  + '\t' + 'Major_tax' +  '\t' + 'Tax' + '\t' + 'CPA_count' + '\t' + 'CPA_binary' + '\t' + 'NhaB_count' + '\t' +  'NhaB_binary' + '\t' + 'NhaC_count' +'\t' + 'NhaC_binary' +'\t' + 'NhaD_count' + '\t' + 'NhaD_binary' + '\n'
        Out.write(header)
        for line in Meta:
            species_name = (line.split('\t')[1])
            if species_name == 'Agaricus_bisporus_ var. bisporus H97':
                species_name = 'Agaricus_bisporus_'
            if species_name == 'Allomyces_macrogynus_ATCC 38327':
                species_name = 'Allomyces_macrogynus_ATCC'
            if species_name == 'Pecoramyces ruminantium':
                species_name = 'Pecoramyces'
            Large_grouping = line.split('\t')[11]
            Tax = (line.split('\t')[12])
            NhaB_count = NhaB_binary = NhaC_count = NhaC_binary = NhaD_count = NhaD_binary = CPA_count = CPA_binary= 0
            if species_name  in NhaD_list:
                NhaD_count = NhaD_list.count(species_name)
                NhaD_binary = 1
            if species_name  in NhaC_list:
                NhaC_count = NhaC_list.count(species_name)
                NhaC_binary = 1
            if species_name  in NhaB_list:
                NhaB_count = NhaB_list.count(species_name)
                NhaB_binary = 1
            if species_name  in CPA_list_forIT:
                CPA_count = CPA_list_forIT.count(species_name)
                CPA_binary = 1
            if species_name == 'Agaricus_bisporus_':
                species_name = 'Agaricus_bisporus_ var. bisporus H97'
            if species_name == 'Allomyces_macrogynus_ATCC':
                species_name = 'Allomyces_macrogynus_ATCC 38327'
            if species_name == 'Pecoramyces':
                species_name = 'Pecoramyces ruminantium'
            Wline = species_name  + '\t' + Large_grouping + '\t' + Tax + '\t' + str(CPA_count) + '\t' + str(CPA_binary) + '\t' + str(NhaB_count) + '\t' +  str(NhaB_binary) + '\t' + str(NhaC_count) +'\t' + str(NhaC_binary) +'\t' + str(NhaD_count) + '\t' + str(NhaD_binary) + '\n'

            Out.write(Wline)
print(len(set(CPA_list_forIT)))
print(len(CPA_list_forIT))
tax = []
import json



with open('/nesi/nobackup/uc04105/new_databases_May/final_tree_set/Eukarya_4nov_passed_all_setcpa.fsta', 'w') as Passed:
    for key in Passed_all_list.keys():
        header = key.split('_tax')[0]
        tax.append(key.split('tax:')[1])
        sequence = Passed_all_list[key]
        line = '>' + header + '\n' + str(sequence) + '\n'
        Passed.write(line)

print(len(set(tax)))
