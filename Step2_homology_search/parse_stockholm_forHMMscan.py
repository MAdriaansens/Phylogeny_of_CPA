import re

import sys

#Craig herbold has been instrumental in generating the 're' part, I do not like running re

hit_file = sys.argv[1] # '/nesi/nobackup/uc04105/results/hmmalign/pipeline/CPA/Archaea_nhaA_TRKACpfam_align_CPAfold.sthk'  #sys.argv[1]  
outfile = sys.argv[2]
length = int(sys.argv[3])

# define alignment
alignment={}



# read in alignment

with open(hit_file) as hmmalignment:
    for line in hmmalignment:
        if(line[0]!="#"):
            splitLine=line.split()

            if(len(splitLine)==2):
                stringFilter = lambda text: re.sub('[.*]', '', splitLine[1])
                filteredString=stringFilter(splitLine[1])

                if(splitLine[0] in alignment.keys()):
                    alignment[splitLine[0]]=alignment[splitLine[0]]+filteredString
                else:
                    alignment[splitLine[0]]=filteredString



# check for identical length



#make a file containing all sequences above a certain length
#one file contains all sequences gapped and one replaces the gaps with X's



#make a file containing all sequences above a certain length
#one file contains all sequences gapped and one replaces the gaps with X's
import re 
#make a file containing all sequences above a certain length
#one file contains all sequences gapped and one replaces the gaps with X's
with open('{}.fasta'.format(outfile), 'w') as file1:
    for key in alignment.keys():
        fasta = alignment[key]
        if sum(1 for c in fasta if c.isupper()) < length:
            pass
        else:   
            sequence_test = alignment[key]
            p1 =re.sub("^[a-z]*", '', sequence_test)
            sequence=re.sub("[a-z]*$", '', p1).replace('-', '').upper()
            sequence_and_header = '>' + key + '\n' + str(sequence) + '\n'
            file1.write(sequence_and_header)
file1.close()

