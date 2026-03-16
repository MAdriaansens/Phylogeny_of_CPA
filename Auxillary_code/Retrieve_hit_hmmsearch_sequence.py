

import sys
import os
from Bio import SeqIO


#Query_HMM = #sys.argv[1]  #EUK, Bacteria, ARC or PF00999
#Bac_hmmsearch_seqdir= #sys.argv[2] #'/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/Pfam'

Seqs = {}

hmmsearch_domtbl=sys.argv[2] #'/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMdombtlout/Archaea/Archaea_vsArchaea_subset_1_fl_domtblout_PF009999.tsv'
seq_file = sys.argv[1] #'/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Archaea/cross_domain/All_arc_vsArchaea_subset1_fulllength.fasta'

outfile = sys.argv[3]

for record in SeqIO.parse(seq_file, 'fasta'):
    Seqs[record.id] = str(record.seq)

#now we have a directory containing all the proteins hit by the HMMsearch

#now I generate a dictionairy of each protein id and at which residues the alignments to the query HMM starts and ends (Bac_borders)

Borders = {}

with open(hmmsearch_domtbl, 'r') as lines:
    end = ''
    start = ''
    for line in lines:
        if line[0] != '#':
            start = line.split('  ')[-5]
            end =line.split('  ')[-4]
            if len(end) < 2:
                start = line.split('  ')[-7]
                end = line.split('  ')[-5]

            if len(start) <1:
                start = line.split('  ')[-6]
                end = line.split('  ')[-4]

            diff = int(end) - int(start)
            if diff < 1:
                start = line.split('  ')[-4]
                end = line.split('  ')[-3]
                diff = int(end) - int(start)

            start_list = [start, end]
            Borders[line.split(' ')[0]] = start_list

#this print is ran to make sure there are no empty dicts
print(len(Borders.keys()))

#here we write the fasta files
with open('{}'.format(outfile), 'w') as A:
    for entry in Borders.keys():
        if entry not in list(Seqs.keys()):
            pass
        else:
            start = int(Borders[entry][0])-1
            end = int(Borders[entry][1])-1
            sequence = Seqs[entry]
            header = '>'  + '{}'.format(entry) +'_subset_{}_until_{}'.format(start,end) + '\n' + sequence + '\n'
            A.write(header)
