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


file = open("{}.fa".format(outfile), "a")
x_file = open('{}_x.faa'.format(outfile), 'a')
for header in alignment.keys():
    if (len(alignment[header].replace('-', ''))) > length: 
        
    #    file = open("{}.fa".format(outfile), "a")
        sequence = ">" + header + '\n' + alignment[header] + '\n'
        x_sequence = '>' + header + '\n' + str(alignment[header]).replace('-', 'X') + '\n'
        file.write(sequence)
        x_file.write(x_sequence)
x_file.close()    
file.close()
