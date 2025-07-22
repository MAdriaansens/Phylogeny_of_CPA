import os 
from Bio import SeqIO

HMMsearch_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/hmmsearch_for_scan'

Bac_seqs = {}
#for PFAM
Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/Pfam'

#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/INH_HMM/main' #for Bacteria INHMM

#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/cross_domain/ARC' #for Archaea INHMM

#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/cross_domain/EUK' #for Eukarya INHMM
for search_file in os.listdir(Bac_hmmsearch_seqdir):
    if search_file[-1] == 'a':
        for record in SeqIO.parse('{}/{}'.format(Bac_hmmsearch_seqdir, search_file), 'fasta'):
            Bac_seqs[record.id.split('_tax')[0]] = str(record.seq)

Bac_borders = {}
for search_tsv in os.listdir(HMMsearch_dir):
    if search_tsv.split('7juli')[0] == 'All_Bac':
       if search_tsv.split('.tsv')[0].split('vs')[-1] == 'PF00999':
           #print(search_tsv.split('.tsv')[0].split('vs')[-1])
           with open('{}/{}'.format(HMMsearch_dir,search_tsv), 'r') as lines:
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
                       Bac_borders[line.split(' ')[0]] = start_list


print(len(Bac_seqs.keys()))
print(len(Bac_borders))

with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/subset_forhmmscan/PF00999HMMsearch_matchedseq_vsBacteria.fasta', 'w') as A:
    for entry in Arc_borders.keys():
        if entry not in list(Arc_seqs.keys()):
            pass
        else:
            start = int(Arc_borders[entry][0])-1
            end = int(Arc_borders[entry][1])-1
            sequence = Arc_seqs[entry] 
            header = '>'  + '{}'.format(entry) +'_subset_{}_until_{}'.format(start,end) + '\n' + sequence + '\n'
            A.write(header)
