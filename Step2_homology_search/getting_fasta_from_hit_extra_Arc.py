import csv
import sys
import os
from Bio import SeqIO
hit_file = sys.argv[1]
hit_file_type = sys.argv[2]
fasta_tsv_file = sys.argv[3]
output = sys.argv[4]
#check if eukarya makes snese
#this just filters only unqiue hits

def unique(hits):
    # initialize a null list
    unique_list = []

    # traverse for all elements
    for x in hits:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return(unique_list)


with open("{}".format(hit_file)) as infile:
    next(infile, None)
    hits = []
    for line in infile:
        if hit_file_type == "MMSEQ":
            i = line.split("\t")
            hits.append(i[1])
        if hit_file_type == "HMM":
            i = line.split(" ")
            hits.append(i[0])

unique_list = unique(hits)
print('length list: {}'.format(len(unique_list)))

Precentage_list = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
with open("{}".format(fasta_tsv_file)) as infile:
    file1 = open("{}".format(output), "a")  # append mode
    count = 0
    for line in infile:
        entry = line.split("\t")
        if entry[0] not in unique_list:
            pass
        else:
            if len(entry) < 5:
                pass
            else:
                sequence = ('>'+ entry[0] + '_tax:' +  entry[2] + '\n' + entry[-1])
                file1.write(sequence)
                count = count + 1
                if (count/len(unique_list)) in Precentage_list:
                    print((count/len(unique_list))*10)
                else:
                    pass
    file1.close()
