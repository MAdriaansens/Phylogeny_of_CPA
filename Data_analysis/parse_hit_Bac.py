#final 
from Bio import SeqIO
def get_hits(direction, final_file):
    hit_list = []

         
    for record in SeqIO.parse('{}/{}'.format(direction, final_file), 'fasta'):
            hit_list.append(record.id)
    return(hit_list)
   #returns all hits from HMMsearch and MMseq and does not filter them

def parse_tax(in_list):
    uniq_list = []
    for i in (list(set(in_list))):
        uniq_list.append(i.split('tax:')[1])
    return(uniq_list)
    #reduces the list to unique values and then takes the taxonomic component, this means that for some taxa multiple entries are present
    #we use this to count the occurance of a taxa in the list

def parse_id(in_list):
    uniq_id_list = []
    for i in (list(set(in_list))):
        uniq_id_list.append(i.split('_tax:')[0])
    return(uniq_id_list)
    #this def returns the unique ids

direction ='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999/merged'
final_file='PF00999_merged_Bacteria.fasta'

CPAhit_list = get_hits(direction,final_file)
CPAtax_list = parse_tax(CPAhit_list)
CPA_uniq_ids = parse_id(CPAhit_list)
print(len(CPA_uniq_ids))
print(len(list(set(CPAtax_list))))

import os
#final 
from Bio import SeqIO
def get_hits(direction, final_file):
    hit_list = []

         
    for record in SeqIO.parse('{}/{}'.format(direction, final_file), 'fasta'):
        hit_id = record.id.split('|')[1]
        hit_list.append(hit_id)
    return(hit_list)
   #returns all hits from HMMsearch and MMseq and does not filter them



direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF06450/MERGED'
final_file='PF06450_Bacteria_all_merged_aligned.fasta'
NhaBhit_list = get_hits(direction,final_file)
NhaBtax_list = parse_tax(NhaBhit_list)
NhaB_uniq_ids = parse_id(NhaBhit_list)
print(len(NhaB_uniq_ids))
print(len(NhaBtax_list))


direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03600/cross_domain'
final_file='PF03600_Bacteria_crossdomain_merged.fasta'
NhaDhit_list = get_hits(direction,final_file)
NhaDtax_list = parse_tax(NhaDhit_list)
NhaD_uniq_ids = parse_id(NhaDhit_list)
print(NhaD_uniq_ids[55])
print(len(NhaDtax_list))

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03553/MERGED'
final_file='PF03553_Bacteria_all_merged_aligned.fasta'
NhaChit_list = get_hits(direction,final_file)
NhaCtax_list = parse_tax(NhaChit_list)
NhaC_uniq_ids = parse_id(NhaChit_list)




homedir ='/nesi/nobackup/uc04105/new_databases_May/GTDB_226'
with open('{}/IT_bacteria_9June.tsv'.format(homedir), 'a') as IT:
    IT.write('Acession' + '\t' + 'completeness' + '\t' + 'contamination' + '\t' + 'GTDB_taxonomy' + '\t' + 'CPA_count' + '\t' + 'NhaB_count' + '\t' + 'NhaC_count' + '\t' + 'NhaD_count' + '\n')
    with open('{}/bac120_metadata.tsv'.format(homedir), 'r') as B:
        for line in B:
            if (line.split('\t')[18]) == 't':

                completeness = line.split('\t')[2]
                contamination = line.split('\t')[3]
                GTDB_taxonomy =  line.split('\t')[19].replace(' ', '_')
                Accession = line.split('\t')[0]
                NhaD_count = str(NhaDtax_list.count(GTDB_taxonomy))
                NhaC_count = str(NhaCtax_list.count(GTDB_taxonomy))
                NhaB_count = str(NhaBtax_list.count(GTDB_taxonomy))
                CPA_count = str(CPAtax_list.count(GTDB_taxonomy))
                outline = Accession + '\t' + completeness + '\t' + contamination + '\t' + GTDB_taxonomy + '\t' + CPA_count + '\t' + NhaB_count + '\t' + NhaC_count + '\t' + NhaD_count + '\n'
                IT.write(outline)
IT.close()
B.close()
