import sys

in_filterd = sys.argv[1] #'PF03600MMseq_e03vsBacteria_10_alignedPF03600.faa.fasta'
out_fl = sys.argv[2] #'PF03600MMseq_e03vsBacteria_10_alignedPF03600_fl.fasta'


outdir =  sys.argv[3]#'/nesi/nobackup/uc04105/new_databases_May/GTDB_226/results/HMMalign/Bacteria/PF03600/MMseq'

tsv =  sys.argv[4] #'/nesi/nobackup/uc04105/new_databases_May/GTDB_226/aDB/tsv/Bacteria_GTDB226_protein_May92025_chunk_10.tsv'

from Bio import SeqIO

hit_list = []
for record in SeqIO.parse('{}/{}'.format(outdir, in_filterd), 'fasta'):
#be aware if '|' is present, in the future maybe add a element which checks for it
    og_id = (record.id.split('_tax:')[0]) #.split('|')[1])

    hit_list.append(og_id)

with open('{}/{}'.format(outdir, out_fl), 'a') as O:
    with open(tsv, 'r') as T:
        next(T, None)
        for line in T:
            if line.split('\t')[0] in hit_list:
                out_fasta = '>' + line.split('\t')[0] + '_tax:' + line.split('\t')[2] + '\n' + line.split('\t')[-1]
                O.write(out_fasta)
