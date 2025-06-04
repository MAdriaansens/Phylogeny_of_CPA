from Bio import SeqIO
direction = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May'
hits_cpa = 'EukPF00999_hits_19may.tsv' #previous hits

def get_hits(direction, MMseqfile, HMMfile):
    hit_list = []

    for record in SeqIO.parse('{}/{}'.format(direction, MMseqfile), 'fasta'):
        hit_list.append(record.id)

            
    for record in SeqIO.parse('{}/{}'.format(direction, HMMfile), 'fasta'):
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

#NhaB 

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF06450'
HMMfile='PF06450HMM_e03vsEukarya_PF06450_alignedPF06450.fasta.fasta'
MMseqfile='PF06450MMseq_vsEukarya_aligned_PF06450.fasta.fasta'
in_list = get_hits(HMMalign_dir, MMseqfile, HMMfile)
NhaB_tax_list = parse_tax(in_list)

#NhaC 

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF03553'
HMMfile='PF03553HMM_e03vsEukarya_PF03553_alignedPF03553.fasta.fasta'
MMseqfile='PF03553MMseq_vsEukarya_aligned_PF03553.sthk.fasta'

in_list = get_hits(HMMalign_dir, MMseqfile, HMMfile)
NhaC_tax_list = parse_tax(in_list)
print(len(list(set(in_list))))

#NhaD

HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF03600'
HMMfile='PF03600HMM_e03vsEukarya_PF03600_alignedPF03600.fasta.fasta'
MMseqfile='PF03600MMseq_vsEukarya_aligned_PF03600.sthk.fasta'
in_list = get_hits(HMMalign_dir, MMseqfile, HMMfile)
NhaD_tax_list = parse_tax(in_list)






with open('{}/complete_euk_antiporter_4June.tsv'.format(direction), 'w') as IT_tsv:
    IT_tsv.write('species' + '\t' + 'taxonomy' + '\t' + 'group' + '\t' + 'CPA_count' + '\t' + 'NhaB_count' + '\t' + 'NhaC_count' + '\t' + 'NhaD_count' + '\n')
    with open('{}/{}'.format(direction, hits_cpa), 'r') as CPA_tsv:
        next(CPA_tsv, None)
        for line in CPA_tsv:
            line_list = (line.split('\n')[0].split('\t'))
            species = line_list[0]
            taxonomy = line_list[1]
            group = line_list[2]
            gene_count = line_list[3]
            CPA_count = line_list[4]
            CPA_binary = line_list[5]
            NhaD_count = str(NhaD_tax_list.count(species))
            NhaC_count = str(NhaC_tax_list.count(species))
            NhaB_count = str(NhaB_tax_list.count(species))
            outline = species + '\t'+ taxonomy + '\t' + group + '\t' + CPA_count + '\t' + NhaB_count + '\t' + NhaC_count + '\t' + NhaD_count + '\n'
            IT_tsv.write(outline)


#this part allows for the retrieval of full length sequences
tsv = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/Euk_db_May_protein.tsv'


HMMalign_dir = '/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/results/HMMalign/PF03553'
HMMfile='PF03553HMM_e03vsEukarya_PF03553_alignedPF03553.fasta.fasta'
MMseqfile='PF03553MMseq_vsEukarya_aligned_PF03553.sthk.fasta'

in_list = get_hits(HMMalign_dir, MMseqfile, HMMfile)
NhaC_ids = list(set(parse_id(in_list)))
print(len(NhaC_ids))

with open('/nesi/nobackup/uc04105/new_databases_May/Euk_database_May/all_hits_fl/NhaC_Euk.fasta', 'w') as NhaCout:
    with open(tsv, 'r') as T:
        next(T, None)
        for line in T:
            if line.split('\t')[0] in NhaC_ids:
                out_fasta = '>' + line.split('\t')[0] + '_tax:_' + line.split('\t')[1] + '\n' + line.split('\t')[-1]
                NhaCout.write(out_fasta)
