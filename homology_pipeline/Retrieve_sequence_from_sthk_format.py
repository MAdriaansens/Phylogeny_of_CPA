import re
import sys

#Craig herbold has been instrumental in generating the 're' part, I do not like running re

hit_file = sys.argv[1] # '/nesi/nobackup/uc04105/results/hmmalign/pipeline/CPA/Archaea_nhaA_TRKACpfam_align_CPAfold.sthk'  #sys.argv[1]  
outfile = sys.argv[2]
# define alignment
length = int(sys.argv[3])
alignment={}



# read in alignment

with open(hit_file) as hmmalignment:
    for line in hmmalignment:
        if(line[0]!="#"):
            splitLine=line.split()

            if(len(splitLine)==2):
                stringFilter = lambda text: re.sub('[a-z.*]', '', splitLine[1])
                filteredString=stringFilter(splitLine[1])

                if(splitLine[0] in alignment.keys()):
                    alignment[splitLine[0]]=alignment[splitLine[0]]+filteredString
                else:
                    alignment[splitLine[0]]=filteredString



# check for identical length

lengthDict={}

for header in alignment.keys():
    lengthDict[len(alignment[header])]=1



if len(lengthDict) > 1:
    sys.exit("ERROR! Aligned sequences differ by length")


#make a file containing all sequences above a certain length
#one file contains all sequences gapped and one replaces the gaps with X's
with open('{}.fasta'.format(outfile), 'a') as file1:
    with open('{}_exes.fasta'.format(outfile), 'a') as Xfile1:
        for key in alignment.keys():
            fasta = alignment[key]
            if len(fasta.replace('-', '')) > length:
                sequence = '>' + key + '\n' + fasta + '\n'
                x_sequence = '>' + key + '\n' + fasta.replace('-', 'X') + '\n'
                file1.write(sequence)
                Xfile1.write(x_sequence)
        Xfile1.close()
        file1.close()
