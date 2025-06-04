import os
from Bio import SeqIO
#note that Euk differs, so check out code for that
def get_hits(direction, mmseq_dir, hmm_dir):
    hit_list = []


    MMseqdir_list = os.listdir('{}/{}'.format(direction, mmseq_dir))
    for file in MMseqdir_list:
        if 'faa.fasta' not in file:
            pass
        else:
            for record in SeqIO.parse('{}/{}/{}'.format(direction, mmseq_dir, file), 'fasta'):
                Cpahit_list.append(record.id)
        HMMsearch_list = os.listdir('{}/{}'.format(direction, hmm_dir))
    for file in HMMsearch_list:
        if 'faa.fasta' not in file:
            pass
        else:
            for record in SeqIO.parse('{}/{}/{}'.format(direction, hmm_dir, file), 'fasta'):
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

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999'
mmseq_dir = 'MMseq/Manual_e03'
hmm_dir ='HMMsearch'

Cpahit_list = get_hits(direction, mmseq_dir, hmm_dir)
CPA_tax_list = parse_tax(Cpahit_list)
CPA_uniq_ids = parse_id(Cpahit_list)
#if you want the number of uniq taxa's  run (len(list(set(CPA_tax_list))))

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03600'
mmseq_dir = 'MMseq'
hmm_dir ='HMMsearch'

NhaDhit_list = get_hits(direction, mmseq_dir, hmm_dir)
NhaD_tax_list = parse_tax(NhaDhit_list)
NhaD_uniq_ids = parse_id(NhaDhit_list)

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03553'
mmseq_dir = 'MMseq'
hmm_dir ='HMMsearch'

NhaChit_list = get_hits(direction, mmseq_dir, hmm_dir)
NhaC_tax_list = parse_tax(NhaChit_list)
NhaC_uniq_ids = parse_id(NhaChit_list)

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF06450'
mmseq_dir = 'MMseq'
hmm_dir ='HMMsearch'

NhaBhit_list = get_hits(direction, mmseq_dir, hmm_dir)
NhaB_tax_list = parse_tax(NhaBhit_list)
NhaB_uniq_ids = parse_id(NhaBhit_list)

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
