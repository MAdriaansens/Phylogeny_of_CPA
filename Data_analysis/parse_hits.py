import os
from Bio import SeqIO

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF00999'

MMseqdir_list = os.listdir('{}/MMseq/Manual_e03'.format(direction))
Cpahit_list = []
for file in MMseqdir_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/MMseq/Manual_e03{}'.format(direction, file), 'fasta'):
            Cpahit_list.append(record.id)
HMMsearch_list = os.listdir('{}/HMMsearch'.format(direction))
for file in HMMsearch_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/HMMsearch/{}'.format(direction, file), 'fasta'):
            Cpahit_list.append(record.id)

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03553'

MMseqdir_list = os.listdir('{}/MMseq'.format(direction))
NhaChit_list = []
for file in MMseqdir_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/MMseq/{}'.format(direction, file), 'fasta'):
            NhaChit_list.append(record.id)
HMMsearch_list = os.listdir('{}/HMMsearch'.format(direction))
for file in HMMsearch_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/HMMsearch/{}'.format(direction, file), 'fasta'):
            NhaChit_list.append(record.id)

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF06450'

MMseqdir_list = os.listdir('{}/MMseq'.format(direction))
NhaBhit_list = []
for file in MMseqdir_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/MMseq/{}'.format(direction, file), 'fasta'):
            NhaBhit_list.append(record.id)
HMMsearch_list = os.listdir('{}/HMMsearch'.format(direction))
for file in HMMsearch_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/HMMsearch/{}'.format(direction, file), 'fasta'):
            NhaBhit_list.append(record.id)

direction = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03600'

MMseqdir_list = os.listdir('{}/MMseq'.format(direction))
NhaDhit_list = []
for file in MMseqdir_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/MMseq/{}'.format(direction, file), 'fasta'):
            NhaDhit_list.append(record.id)
HMMsearch_list = os.listdir('{}/HMMsearch'.format(direction))
for file in HMMsearch_list:
    if 'faa.fasta' not in file:
        pass
    else:
        for record in SeqIO.parse('{}/HMMsearch/{}'.format(direction, file), 'fasta'):
            NhaDhit_list.append(record.id)

NhaD_uniq_list = []
NhaB_uniq_list = []
NhaC_uniq_list = []
CPA_uniq_list = []
def parse_tax(in_list):
    uniq_list = []
    for i in (list(set(in_list))):
        uniq_list.append(i.split('tax:')[1])
    return(uniq_list)
NhaD_uniq_list = parse_tax(NhaDhit_list)
NhaC_uniq_list = parse_tax(NhaChit_list)
NhaB_uniq_list = parse_tax(NhaBhit_list)
CPA_uniq_list = parse_tax(Cpahit_list)

print(NhaD_uniq_list[7])
print(len(list(set(CPA_uniq_list))))
