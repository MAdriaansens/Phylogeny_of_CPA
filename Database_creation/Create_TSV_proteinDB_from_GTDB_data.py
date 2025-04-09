import os
from Bio import SeqIO


#this requires that you got the metadata as well as the whole GTDB proteome
#the tsv can then be run to generate the db in a fasta format

tsv_file = '/nesi/nobackup/uc04105/database/GTDB_220/ar53_metadata_r220.tsv'
dir_files = '/nesi/nobackup/uc04105/database/GTDB_220/protein_faa_reps/archaea'
files = os.listdir('/nesi/nobackup/uc04105/database/GTDB_220/protein_faa_reps/archaea')

outfile = '/nesi/nobackup/uc04105/database/GTDB_220/Archaea_GTDBprotein_0608024.tsv'

dictionary = {}

with open('{}'.format(tsv_file)) as dataset:
    for i in dataset:

        if i.split('\t')[17] != 'gtdb_genome_representative':
            dictionary[(i.split('\t')[17])] = (i.split('\t')[19]).replace(' ', '_')


header = 'GTDB_id' + '\t' + 'GTDB_taxonomy' + '\t' + 'Protein_id_GTDB' + '\t' + 'Sequence' + '\t' + 'Protein_id_annotation' + '\n'


with open(outfile, 'a') as f:
    f.write(header)
f.close()


with open(outfile, 'a') as f:
    count = 0
    for file in files:

        GTDB_ID = file.split('_protein.faa')[0]
        taxonomy = dictionary[GTDB_ID]
        for record in SeqIO.parse('{}/{}'.format(dir_files, file), "fasta"):
            count = count + 1
            Protein_id_GTDB = record.id
            Sequence = str(record.seq)
            Protein_id_annotation = 'Arc_{}'.format(count)
            line = Protein_id_annotation + '\t' + Protein_id_GTDB + '\t' + taxonomy + '\t' + GTDB_ID + '\t' + Sequence + '\n'
            f.write(line)
f.close()
