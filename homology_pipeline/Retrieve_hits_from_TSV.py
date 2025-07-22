import csv
import sys

hit_file = sys.argv[1]
hit_file_type = sys.argv[2]
fasta_tsv_file = sys.argv[3]
output = sys.argv[4]
domain = sys.argv[5]
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
           i = line.split("\t")[1]
           number = int(i.split('_')[-1])
           hits.append(number)
       if hit_file_type == "HMM":
           i = line.split(" ")
           if '#' not in i[0]:
               number = int(i[0].split('_')[-1])
               hits.append(number)

unique_list = unique(hits)
unique_list=sorted(unique_list,reverse=True)
id_list = []

for entry in unique_list:
    edit_id = '{}_'.format(domain) + str(entry)
    id_list.append(edit_id)

with open("{}".format(fasta_tsv_file)) as infile:
    file1 = open("{}".format(output), "a")  # append mode
    count = 0
    for line in infile:
        entry = line.split("\t")
        if entry[0] not in id_list:
            pass
        else:
           header =  str(entry[0])
           header = header.replace(" ", "_")
           sequence = ('>'+ entry[0] + '\n' + entry[-1]+'\n')
           file1.write(sequence)

    #            print(sequence)
file1.close()
