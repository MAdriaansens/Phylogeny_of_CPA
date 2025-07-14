
import json
from Bio import SeqIO

HMS03553='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea_PF03553_aligned_hmmscanned_onlyPF03553.json'
HMS03600='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea_PF03600_aligned_hmmscanned_onlyPF03600.json'
HMS06450='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMscan/Archaea_PF06450_aligned_hmmscanned_onlyPF06450.json'

with open(HMS03600) as json_file:
    NhaD_data=json.load(json_file)
NhaD_tax_list = []
for key in NhaD_data.keys():
    output = NhaD_data[key]
    if output[0] == 'PF03600':
        NhaD_tax_list.append(key.split('_tax:')[-1])
print(len(NhaD_tax_list))

with open(HMS03553) as json_file:
    NhaC_data=json.load(json_file)
NhaC_tax_list = []
for key in NhaC_data.keys():
    output = NhaC_data[key]
    if output[0] == 'PF03553':
        NhaC_tax_list.append(key.split('_tax:')[-1])
print(len(list(set(NhaC_tax_list))))

with open(HMS06450) as json_file:
    NhaB_data=json.load(json_file)

NhaB_tax_list = []
for key in NhaB_data.keys():
    output = NhaB_data[key]
    if output[0] == 'PF06450':
        NhaB_tax_list.append(key.split('_tax:')[-1])
print(len(list(set(NhaB_tax_list))))

CPA_hit_list=[]
for record in SeqIO.parse('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/merged_PFAM_HMMsearches_before_self_cross/PF00999_Archaea_MERGED_alignedPF00999_parsed_aligned258.faa_exes.fasta'', 'fasta'):
    hit_id = record.id.split('|')[1]
    CPA_hit_list.append(hit_id)

CPA_tax_list=[]
for i in CPA_hit_list:
    uniq_id = i.split('tax:')[1]
    CPA_tax_list.append(uniq_id)


homedir ='/nesi/nobackup/uc04105/new_databases_May/GTDB_226'
with open('{}/IT_bacteria_28may.tsv'.format(homedir), 'a') as IT:
    IT.write('Acession' + '\t' + 'completeness' + '\t' + 'contamination' + '\t' + 'GTDB_taxonomy' + '\t' + 'CPA_count' + '\t' + 'NhaB_count' + '\t' + 'NhaC_count' + '\t' + 'NhaD_count' + '\n')
    with open('{}/bac120_metadata.tsv'.format(homedir), 'r') as B:
        for line in B:
            if (line.split('\t')[18]) == 't':

                completeness = line.split('\t')[2]
                contamination = line.split('\t')[3]
                GTDB_taxonomy =  line.split('\t')[19].replace(' ', '_')
                Accession = line.split('\t')[0]
                NhaD_count = str(NhaD_tax_list.count(GTDB_taxonomy))
                NhaC_count = str(NhaC_tax_list.count(GTDB_taxonomy))
                NhaB_count = str(NhaB_tax_list.count(GTDB_taxonomy))
                CPA_count = str(CPA_tax_list.count(GTDB_taxonomy))
                outline = Accession + '\t' + completeness + '\t' + contamination + '\t' + GTDB_taxonomy + '\t' + CPA_count + '\t' + NhaB_count + '\t' + NhaC_count + '\t' + NhaD_count + '\n'
                IT.write(outline)
IT.close()
B.close()
