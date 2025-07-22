#explanation
#The goal of this code is to extract, from a list of HMMsearch output files, the residues which match with a given HMM.
#In total 4 different queries were used against Bacteria, Archaea and Eukarya. In Bacteria the files are subsetted into 1/61th of all GTDB Bacteria proteoemes

#Three files are of importance:
#1 the fasta file fullength containing sequences hit with the HMMsearch without the --domtblout output with the same evalue threshold
#2 the hmmsearch .tsv output files, run with the -E 0.001, --cpu 10 and --domtblout as output. 
    #All_Bac means that it All bacteria sequences were the subject of the initial HMMsearch, there exists a set of files called All_Bac.., All_Euk.. and All_Arc
    #what comes after the vs denotes the query, this includes ARC (Archaea), PF00999, Bacteria and Euks (Eukarya

#3 the output fasta file of this script containing only the part of the sequences hit by the HMM, regardless of size. 


import os 
from Bio import SeqIO


#this sets the directory for the HMMsearch done with --domtblout, the --domtblout out contains the start and end of the HMMalignment to the given sequences

HMMsearch_dir = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/hmmsearch_for_scan'


#set up a dictionairy of Bacteria protein ids as Key and sequences as Values. Hence the name Bac_seq
Bac_seqs = {}

#then in a previous HMMsearch we extracted the sequences which exceeded the e-value threshold of 1*10^-3. 

#these are the fasta files containing sequences hit by HMMsearch, I only retrieve one query at a time. 
#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/INH_HMM/main' #for Bacteria INHMM
#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/cross_domain/ARC' #for Archaea INHMM
#Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/cross_domain/EUK' #for Eukarya INHMM


#In this case we use PFAM as query as an example
Bac_hmmsearch_seqdir='/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/Bacteria/PF00999/Pfam'

#since I have subsetted the GTDB into 61 and the ouput is put into a directory, we first need to pool all the sequences together.

for search_file in os.listdir(Bac_hmmsearch_seqdir):
    if search_file[-1] == 'a': #in this directory both the .tsv file and .fasta files are present, this line selects for the fasta file
        for record in SeqIO.parse('{}/{}'.format(Bac_hmmsearch_seqdir, search_file), 'fasta'):
            Bac_seqs[record.id.split('_tax')[0]] = str(record.seq)


#now we have a directory containing all the proteins hit by the HMMsearch

#now I generate a dictionairy of each protein id and at which residues the alignments to the query HMM starts and ends (Bac_borders)

Bac_borders = {}
for search_tsv in os.listdir(HMMsearch_dir):

    #All_Bac means that it All bacteria sequences were the subject of the initial HMMsearch 
    if search_tsv.split('7juli')[0] == 'All_Bac':
        #what comes after the vs denotes the query
       if search_tsv.split('.tsv')[0].split('vs')[-1] == 'PF00999':
           #print(search_tsv.split('.tsv')[0].split('vs')[-1])
           with open('{}/{}'.format(HMMsearch_dir,search_tsv), 'r') as lines:
               for line in lines:
                   #since HMMsearch.tsv is not tab delimited retrieving the start and end is akward but this solves it.
                   #I have manually parsed the start and end date of the smalles Archaea set, line for line in HMMsearch to confirm that this works.
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

#this print is ran to make sure there are no empty dicts
print(len(Bac_seqs.keys()))
print(len(Bac_borders))

#here we write the fasta files
with open('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMsearch/subset_forhmmscan/PF00999HMMsearch_matchedseq_vsBacteria.fasta', 'w') as A:
    for entry in Bac_borders.keys():
        if entry not in list(Bac_seqs.keys()):
            pass
        else:
            start = int(Bac_borders[entry][0])-1
            end = int(Bac_borders[entry][1])-1
            sequence = Bac_seqs[entry] 
            header = '>'  + '{}'.format(entry) +'_subset_{}_until_{}'.format(start,end) + '\n' + sequence + '\n'
            A.write(header)
