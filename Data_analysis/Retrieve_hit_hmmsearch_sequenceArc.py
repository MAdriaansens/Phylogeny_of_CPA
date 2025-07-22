import os 
from Bio import SeqIO
#Arc_hmmsearch_seq='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea/INH_HMM/Archaea_Manual_e5_cov50_AlignedPF00999_clustered_0.6_rep_seq_Ginsialigned_HMM_e03vsArchaea.fa for Archaea

#non detailed refers to the original HMMsearch not using the domtblout option and thus missing some details, from this hmmsearch the fl sequences were retreived
#Arc_hmmsearch_seq='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea/non_detailed/Archaea_PF00999HMM_e03.faa' for PF00999 search
#Arc_hmmsearch_seq='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea/non_detailed/BAC_HMM_e03vsArchaea.fa' for Bacteria
#Arc_hmmsearch_seq='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea/non_detailed/EUK_HMM_e03vsArchaea.fa' for Eukarya

Arc_seqs = {}

Arc_hmmsearch_seq='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea/non_detailed/BAC_HMM_e03vsArchaea.fa'

for record in SeqIO.parse(Arc_hmmsearch_seq, 'fasta'):
    Arc_seqs[record.id.split('_tax')[0]] = str(record.seq)



HMMsearch_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Archaea' #for 


Arc_borders = {}
for search_tsv in os.listdir(HMMsearch_dir):
    if search_tsv.split('7juli')[0] == 'All_Arc':
       if search_tsv.split('.tsv')[0].split('vs')[-1] == 'ARC':
           print(search_tsv.split('.tsv')[0].split('vs')[-1])
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
                       Arc_borders[line.split(' ')[0]] = start_list


print(len(Arc_seqs.keys()))
print(len(Arc_borders))
#write a new fasta file using the part of the sequence hit
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/subset_forhmmscan/ArchaeaHMMsearch_matchedseq_vsArchaea.fasta', 'w') as A:
    for entry in Arc_borders.keys():
        if entry not in list(Arc_seqs.keys()):
            pass
        else:
            start = int(Arc_borders[entry][0])-1
            end = int(Arc_borders[entry][1])-1
            sequence = Arc_seqs[entry] 
            header = '>'  + '{}'.format(entry) +'_subset_{}_until_{}'.format(start,end) + '\n' + sequence + '\n'
            A.write(header)
