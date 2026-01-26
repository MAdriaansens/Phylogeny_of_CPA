from Bio import SeqIO
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open('/nesi/nobackup/uc04105/fasta_files/pfams/{}'.format(infile), 'a') as S:
    for record in SeqIO.parse('/nesi/nobackup/uc04105/fasta_files/pfams/{}'.format(outfile)', 'fasta'):
        if len(record.seq) < 257: #257 in case of PF00999
            pass
        else:
            outline = '>' + record.id + '\n' + str(record.seq) + '\n'
            S.write(outline)
