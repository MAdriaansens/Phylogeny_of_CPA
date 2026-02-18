import os
from Bio import SeqIO
import sys

#this requires that you got the metadata as well as the whole GTDB proteome
#the tsv can then be run to generate the db in a fasta format

tsv_file = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/bac120_metadata.tsv'
dir_files = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/protein_faa_reps/bacteria'
files = os.listdir('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/protein_faa_reps/bacteria')

outfile = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_GTDB226_protein_May92025.tsv'
fasta_out = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_GTDB226_protein_May92025.faa'
dictionary = {}

with open('{}'.format(tsv_file)) as dataset:
    for i in dataset:

        if i.split('\t')[17] != 'gtdb_genome_representative':
            dictionary[(i.split('\t')[17])] = (i.split('\t')[19]).replace(' ', '_')


header = 'Protein_id_annotation' + '\t' + 'Protein_id_GTDB' + '\t' + 'taxonomy' + '\t' + 'GTDB_ID' + '\t' + 'Sequence' + '\n'

with open(outfile, 'a') as f:
   f.write(header)

    with open(fasta_out, 'a') as fo:
        count = 0
        for file in files:
    
            GTDB_ID = file.split('_protein.faa')[0]
            taxonomy = dictionary[GTDB_ID]
            for record in SeqIO.parse('{}/{}'.format(dir_files, file), "fasta"):
                count = count + 1
                Protein_id_GTDB = record.id
                Sequence = str(record.seq)
                Protein_id_annotation = 'Bac226_{}'.format(count)
                line = Protein_id_annotation + '\t' + Protein_id_GTDB + '\t' + taxonomy + '\t' + GTDB_ID + '\t' + Sequence + '\n'
                fasta_line = '>' + Protein_id_annotation + '\n' + Sequence + '\n'
                #fo.write(fasta_line)
                #f.write(line)
f.close()
